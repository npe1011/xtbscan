import copy
import os
import shutil
import csv
import gc
from pathlib import Path
import tempfile
import traceback
from typing import Any, List, Union, Optional, Literal

import numpy as np

try:
    import warnings
    warnings.simplefilter('ignore')
    from ase import Atoms
    from ase.units import Hartree
    from ase.io import read, write
    from ase.optimize import LBFGS
    from ase.constraints import FixAtoms, FixInternals
    import torch
    from fairchem.core import FAIRChemCalculator
    from fairchem.core.units.mlip_unit import load_predict_unit
except ImportError:
    IS_UMA_VALID = False
    warnings = None
    Atoms = None
    Hartree = None
    read = None
    write = None
    LBFGS = None
    FixAtoms = None
    FixInternals = None
    torch = None
    FAIRChemCalculator = None
    load_predict_unit = None
else:
    IS_UMA_VALID = True

import config
from xtbscan import xyzutils, saddle

# Global variable to check whether setenv run or not
CHECK_SETENV_UMA = False


class UMAParams:
    def __init__(self,
                 charge: int = 0,
                 mult: int = 1,
                 fmax: float = 0.02,
                 max_cycles: int = 1000):
        # check parameters
        try:
            assert mult >= 1
            assert fmax > 0.0
            assert max_cycles > 2

        except AssertionError as e:
            raise ValueError('Given UMA parameters are not valid. ' + str(e.args))

        self.charge: int = charge
        self.mult: int = mult
        self.fmax: float = fmax
        self.max_cycles: int = max_cycles


def get_uma_calculator(model_path: Path, device: Literal['cuda', 'cpu']) -> FAIRChemCalculator:
    """
    Return new UMA calculator; the instance should be managed by the main module.
    """
    model_path = str(model_path)
    print(f'Loading {model_path} and setup uma calculator...')
    uma_predictor = load_predict_unit(path=model_path, device=device)
    calculator = FAIRChemCalculator(uma_predictor, task_name='omol')
    print('Calculator is ready.')
    return calculator


class UMATerminationError(Exception):
    pass


class UMAConstrain:

    def __init__(self, constrain_type: str, atom_indices: Union[np.ndarray, List[int]], value: Any = None):
        assert constrain_type in ['atoms', 'distance', 'angle', 'dihedral']
        if constrain_type == 'distance':
            assert len(atom_indices) == 2
        if constrain_type == 'angle':
            assert len(atom_indices) == 3
        if constrain_type == 'dihedral':
            assert len(atom_indices) == 4

        self.constrain_type: str = constrain_type
        self.atom_indices: List[int] = list(atom_indices)
        if value is None:
            self.value = None
        elif str(value).lower().strip() in ['', 'a', 'auto']:
            self.value = None
        else:
            self.value = float(value)

    def get_internal(self):
        """
        For distance/angle/dihedral, return [value, [indices]]
        For atoms, return None
        """
        if self.constrain_type == 'atoms':
            return None
        return [self.value, self.atom_indices]


    def get_print_name(self, atoms=None) -> Optional[str]:
        pass


class UMAScan:
    def __init__(self, scan_type: str, atom_indices: Union[np.ndarray, List[int]],
                 start: float, end: float, num_step: int):

        assert int(num_step) >= 2
        assert scan_type in ['distance', 'angle', 'dihedral']

        if scan_type == 'distance':
            assert len(atom_indices) == 2
        if scan_type == 'angle':
            assert len(atom_indices) == 3
        if scan_type == 'dihedral':
            assert len(atom_indices) == 4

        self.scan_type = scan_type
        self.atom_indices = list(atom_indices)
        self.start = start
        self.end = end
        self.num_step = int(num_step)
        self._values = None

    def get_values(self):
        if self._values is None:
            self._values = np.linspace(float(self.start), float(self.end), int(self.num_step), endpoint=True, dtype=float)
        return self._values

    def calc_real_value(self, coordinates: np.ndarray) -> float:
        if self.scan_type == 'distance':
            return xyzutils.calc_distance(coordinates, self.atom_indices)
        elif self.scan_type == 'angle':
            return xyzutils.calc_angle(coordinates, self.atom_indices)
        elif self.scan_type == 'dihedral':
            return xyzutils.calc_dihedral(coordinates, self.atom_indices)
        else:
            raise RuntimeError('unknown type of scan.')

    def get_print_name(self, atom_symbols=None):
        if atom_symbols is None:
            name = self.scan_type + ' unk.'
        else:
            name = self.scan_type + ' ' + '-'.join([atom_symbols[i].capitalize() + str(i + 1) for i in self.atom_indices])
        unit = ' (ang.)' if self.scan_type == 'distance' else ' (deg.)'
        return name + unit


def prepare_ase_constrain_list(uma_constrains: List[UMAConstrain]) -> List[Union[FixAtoms, FixInternals]]:
    ase_constrains = []

    fix_atom_indices = []
    bonds = []
    angles_deg = []
    dihedrals_deg = []

    for constrain in uma_constrains:
        if constrain.constrain_type == 'atoms':
            fix_atom_indices.extend(constrain.atom_indices)
        elif constrain.constrain_type == 'distance':
            bonds.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'angle':
            angles_deg.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'dihedral':
            dihedrals_deg.append([constrain.value, constrain.atom_indices])
        else:
            raise RuntimeError(f'Unknown constrain type: {constrain.constrain_type}')
    fix_atom_indices = sorted(list(set(fix_atom_indices)))
    if len(bonds) == 0:
        bonds = None
    if len(angles_deg) == 0:
        angles_deg = None
    if len(dihedrals_deg) == 0:
        dihedrals_deg = None
    if len(fix_atom_indices) > 0:
        ase_constrains.append(FixAtoms(indices=fix_atom_indices))
    if bonds or angles_deg or dihedrals_deg:
        ase_constrains.append(FixInternals(bonds=bonds, angles_deg=angles_deg, dihedrals_deg=dihedrals_deg))
    return ase_constrains


def prepare_ase_constrain_list_scan1d(scan: UMAScan, uma_constrains: List[UMAConstrain])\
        -> List[List[Union[FixAtoms, FixInternals]]]:
    ase_constrains_list = []

    fix_atom_indices = []
    bonds = []
    angles_deg = []
    dihedrals_deg = []

    for constrain in uma_constrains:
        if constrain.constrain_type == 'atoms':
            fix_atom_indices.extend(constrain.atom_indices)
        elif constrain.constrain_type == 'distance':
            bonds.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'angle':
            angles_deg.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'dihedral':
            dihedrals_deg.append([constrain.value, constrain.atom_indices])
        else:
            raise RuntimeError(f'Unknown constrain type: {constrain.constrain_type}')
    fix_atom_indices = sorted(list(set(fix_atom_indices)))

    # Iterate with scan parameter
    for value in scan.get_values():
        ase_constrains_scan = []
        bonds_scan = copy.deepcopy(bonds)
        angles_deg_scan = copy.deepcopy(angles_deg)
        dihedrals_deg_scan = copy.deepcopy(dihedrals_deg)
        if scan.scan_type == 'distance':
            bonds_scan.append([value, scan.atom_indices])
        elif scan.scan_type == 'angle':
            angles_deg_scan.append([value, scan.atom_indices])
        elif scan.scan_type == 'dihedral':
            dihedrals_deg_scan.append([value, scan.atom_indices])
        else:
            raise RuntimeError(f'Unknown scan type: {scan.scan_type}')
        if len(bonds_scan) == 0:
            bonds_scan = None
        if len(angles_deg_scan) == 0:
            angles_deg_scan = None
        if len(dihedrals_deg_scan) == 0:
            dihedrals_deg_scan = None
        # At least one parameter must exist for FixInternal
        if len(fix_atom_indices) > 0:
            ase_constrains_scan.append(FixAtoms(indices=fix_atom_indices))
        if bonds_scan or angles_deg_scan or dihedrals_deg_scan:
            ase_constrains_scan.append(FixInternals(bonds=bonds_scan, angles_deg=angles_deg_scan, dihedrals_deg=dihedrals_deg_scan))
        ase_constrains_list.append(ase_constrains_scan)

    return ase_constrains_list


def prepare_ase_constrain_list_scan2d(scan1: UMAScan, scan2: UMAScan, uma_constrains: List[UMAConstrain])\
        -> List[List[List[Union[FixAtoms, FixInternals]]]]:

    ase_constrains_list_list = []

    fix_atom_indices = []
    bonds = []
    angles_deg = []
    dihedrals_deg = []

    for constrain in uma_constrains:
        if constrain.constrain_type == 'atoms':
            fix_atom_indices.extend(constrain.atom_indices)
        elif constrain.constrain_type == 'distance':
            bonds.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'angle':
            angles_deg.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'dihedral':
            dihedrals_deg.append([constrain.value, constrain.atom_indices])
        else:
            raise RuntimeError(f'Unknown constrain type: {constrain.constrain_type}')
    fix_atom_indices = sorted(list(set(fix_atom_indices)))

    # Iterate with 2 scan parameters => list dim(scan1) x dim(scan2)
    for value1 in scan1.get_values():
        ase_constrains_list = []
        for value2 in scan2.get_values():
            ase_constrains_scan = []
            bonds_scan = copy.deepcopy(bonds)
            angles_deg_scan = copy.deepcopy(angles_deg)
            dihedrals_deg_scan = copy.deepcopy(dihedrals_deg)
            for scan, value in [(scan1, value1), (scan2, value2)] :
                if scan.scan_type == 'distance':
                    bonds_scan.append([value, scan.atom_indices])
                elif scan.scan_type == 'angle':
                    angles_deg_scan.append([value, scan.atom_indices])
                elif scan.scan_type == 'dihedral':
                    dihedrals_deg_scan.append([value, scan.atom_indices])
                else:
                    raise RuntimeError(f'Unknown scan type: {scan.scan_type}')
            if len(bonds_scan) == 0:
                bonds_scan = None
            if len(angles_deg_scan) == 0:
                angles_deg_scan = None
            if len(dihedrals_deg_scan) == 0:
                dihedrals_deg_scan = None
            # At least one parameter must exist for FixInternal
            if len(fix_atom_indices) > 0:
                ase_constrains_scan.append(FixAtoms(indices=fix_atom_indices))
            if bonds_scan or angles_deg_scan or dihedrals_deg_scan:
                ase_constrains_scan.append(FixInternals(bonds=bonds_scan,
                                                        angles_deg=angles_deg_scan,
                                                        dihedrals_deg=dihedrals_deg_scan))
            ase_constrains_list.append(ase_constrains_scan)
        ase_constrains_list_list.append(ase_constrains_list)

    return ase_constrains_list_list


def prepare_ase_constrain_list_scan_concerted(scans: List[UMAScan], uma_constrains: List[UMAConstrain])\
        -> List[List[Union[FixAtoms, FixInternals]]]:
    ase_constrains_list = []

    fix_atom_indices = []
    bonds = []
    angles_deg = []
    dihedrals_deg = []

    for constrain in uma_constrains:
        if constrain.constrain_type == 'atoms':
            fix_atom_indices.extend(constrain.atom_indices)
        elif constrain.constrain_type == 'distance':
            bonds.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'angle':
            angles_deg.append([constrain.value, constrain.atom_indices])
        elif constrain.constrain_type == 'dihedral':
            dihedrals_deg.append([constrain.value, constrain.atom_indices])
        else:
            raise RuntimeError(f'Unknown constrain type: {constrain.constrain_type}')
    fix_atom_indices = sorted(list(set(fix_atom_indices)))

    # Iterate with scan parameters in 1d concerted manner.
    # Scan length must be the same for all scan objects (should be checked before here).
    scan_step = scans[0].num_step
    for i in range(scan_step):
        ase_constrains_scan = []
        bonds_scan = copy.deepcopy(bonds)
        angles_deg_scan = copy.deepcopy(angles_deg)
        dihedrals_deg_scan = copy.deepcopy(dihedrals_deg)
        for scan in scans:
            value = scan.get_values()[i]
            if scan.scan_type == 'distance':
                bonds_scan.append([value, scan.atom_indices])
            elif scan.scan_type == 'angle':
                angles_deg_scan.append([value, scan.atom_indices])
            elif scan.scan_type == 'dihedral':
                dihedrals_deg_scan.append([value, scan.atom_indices])
            else:
                raise RuntimeError(f'Unknown scan type: {scan.scan_type}')
        if len(bonds_scan) == 0:
            bonds_scan = None
        if len(angles_deg_scan) == 0:
            angles_deg_scan = None
        if len(dihedrals_deg_scan) == 0:
            dihedrals_deg_scan = None
        # At least one parameter must exist for FixInternal
        if len(fix_atom_indices) > 0:
            ase_constrains_scan.append(FixAtoms(indices=fix_atom_indices))
        if bonds_scan or angles_deg_scan or dihedrals_deg_scan:
            ase_constrains_scan.append(FixInternals(bonds=bonds_scan, angles_deg=angles_deg_scan, dihedrals_deg=dihedrals_deg_scan))
        ase_constrains_list.append(ase_constrains_scan)

    return ase_constrains_list


def umascan(input_xyz_file: Union[str, Path],
            job_name: str,
            calculator: FAIRChemCalculator,
            uma_params: UMAParams,
            scans: List[UMAScan],
            constrains: List[UMAConstrain],
            concerted: bool = False,
            keep_log: int = 0) -> None:

    if not IS_UMA_VALID:
        raise RuntimeError('Import of ASE/torch/fairchem-core failed.')

    input_xyz_file = Path(input_xyz_file).absolute()
    if not input_xyz_file.exists():
        raise FileNotFoundError(str(input_xyz_file) + ' not found.')

    # Run
    if not scans:
        _opt(input_xyz_file, job_name, calculator, uma_params, constrains, keep_log)
    elif len(scans) == 1:
        _scan1d(input_xyz_file, job_name, calculator, uma_params, scans[0], constrains, keep_log)
    elif concerted:
        _scan_concerted(input_xyz_file, job_name, calculator, uma_params, scans, constrains, keep_log)
    elif len(scans) == 2:
        _scan2d(input_xyz_file, job_name, calculator, uma_params, scans[0], scans[1], constrains, keep_log)
    else:
        raise ValueError('Only Opt, 1D scan, 2D scan, or concerted scan is available.')


def _opt(input_xyz_file: Path,
         job_name: str,
         calculator: FAIRChemCalculator,
         params: UMAParams,
         constrains: List[UMAConstrain],
         keep_log: int = 0) -> None:

    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (job_name + '.xyz')
    stop_file: Path = input_xyz_file.parent / (job_name + config.STOP_FILE_SUFFIX)
    workdir = Path(tempfile.mkdtemp(dir=input_xyz_file.parent, prefix=job_name + '_')).absolute()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / config.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    def _check_stop():
        if stop_file.exists():
            raise UMATerminationError('Stopped by user')

    trajfile = 'opt.traj'
    trafxyzfile = 'opt.xyz'
    logfile = 'opt.log'

    # main process
    try:
        atoms = read(config.INIT_XYZ_FILE, index=0)
        atoms.info = {'charge': params.charge, 'spin': params.mult}
        atoms.calc = calculator
        if constrains:
            atoms.set_constraint(prepare_ase_constrain_list(constrains))
            atoms.set_positions(atoms.get_positions(), apply_constraint=True)  # Omajinai
        opt = LBFGS(atoms, trajectory=trajfile, logfile=logfile)
        opt.attach(_check_stop, interval=1)
        opt.run(fmax=params.fmax, steps=params.max_cycles)
        final_energy = atoms.get_potential_energy() / Hartree
        calc_success = True
        del opt
        gc.collect()
        if os.path.exists(trajfile):
            try:
                traj = read(f'{trajfile}@0:')
                write(trafxyzfile, traj, format='xyz')
            except Exception:
                print('Trajectory conversion failed.')
                traceback.print_exc()
        save_xyz_from_atoms(result_xyz_file, atoms, title=f'{final_energy:>20.12f}')
    except UMATerminationError:  # When user stop
        try:
            if os.path.exists(trajfile):
                traj = read(f'{trajfile}@0:')
                write(trafxyzfile, traj, format='xyz')
            save_xyz_from_atoms('last.xyz', atoms, title='USER_STOP')
        finally:
            calc_success = False
        raise UMATerminationError
    except Exception as e:
        raise RuntimeError(f'UMA optimization failed: {workdir}') from e

    # common post-process
    finally:
        try:
            atoms.calc = None
        except:
            pass
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def _scan1d(input_xyz_file: Path,
            job_name: str,
            calculator: FAIRChemCalculator,
            params: UMAParams,
            scan: UMAScan,
            constrains: List[UMAConstrain],
            keep_log: int = 0) -> None:

    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (job_name + '.xyz')
    stop_file: Path = input_xyz_file.parent / (job_name + config.STOP_FILE_SUFFIX)
    result_csv_file: Path = input_xyz_file.parent / (job_name + '.csv')
    workdir = Path(tempfile.mkdtemp(dir=input_xyz_file.parent, prefix=job_name + '_')).absolute()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / config.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    def _check_stop():
        if stop_file.exists():
            raise UMATerminationError('Stopped by user')

    # main process
    try:
        atoms = read(config.INIT_XYZ_FILE, index=0)
        atoms.info = {'charge': params.charge, 'spin': params.mult}
        atoms.calc = calculator
        scan_condition_list = prepare_ase_constrain_list_scan1d(scan, constrains)
        energy_list = []
        coordinates_list = []
        for n, scan_condition in enumerate(scan_condition_list):
            try:
                atoms.set_constraint(scan_condition)
                atoms.set_positions(atoms.get_positions(), apply_constraint=True)  # Needed to apply new constrain
                scan_dir = Path(f'./scan_{n+1}')
                scan_dir.mkdir(parents=True)
                trajfile = str(scan_dir / 'opt.traj')
                trafxyzfile = str(scan_dir / 'opt.xyz')
                logfile = str(scan_dir / 'opt.log')

                opt = LBFGS(atoms, trajectory=trajfile, logfile=logfile)
                opt.attach(_check_stop, interval=1)
                opt.run(fmax=params.fmax, steps=params.max_cycles)
                final_energy = atoms.get_potential_energy() / Hartree

                del opt
                gc.collect()
                if os.path.exists(trajfile):
                    try:
                        traj = read(f'{trajfile}@0:')
                        write(trafxyzfile, traj, format='xyz')
                    except Exception:
                        print('Trajectory conversion failed.')
                        traceback.print_exc()
                coordinates_list.append(np.copy(atoms.positions))
                energy_list.append(final_energy)
            except UMATerminationError:  # When user stop
                try:
                    if os.path.exists(trajfile):
                        traj = read(f'{trajfile}@0:')
                        write(trafxyzfile, traj, format='xyz')
                    save_xyz_from_atoms(str(scan_dir / 'last.xyz'), atoms, title='USER_STOP')
                finally:
                    calc_success = False
                raise UMATerminationError('Stopped by user')
            except Exception as e:
                raise RuntimeError(f'UMA optimization failed: {workdir}, scan {n+1}') from e

        energy_list = np.array(energy_list, dtype=float)
        coordinates_list = np.array(coordinates_list, dtype=float)
        atom_symbols = atoms.get_chemical_symbols()

        relative_energy_list = (energy_list - np.min(energy_list)) * config.HARTREE_TO_KCAL
        assigned_value_list = scan.get_values()
        real_value_list = [scan.calc_real_value(coord) for coord in coordinates_list]
        saddle_check_list = saddle.check_saddle_1d(energy_list)

        # output csv data
        csv_data = [['1d', str(scan.num_step)],
                    ['#', '[scan] ' + scan.get_print_name(atom_symbols), '[real] ' + scan.get_print_name(atom_symbols),
                     'E (au)', 'rel. E (kcal/mol)', 'saddle check']]

        for i in range(len(coordinates_list)):
            csv_data.append([
                '{:}'.format(i + 1),
                '{:.4f}'.format(assigned_value_list[i]),
                '{:.4f}'.format(real_value_list[i]),
                '{:.8f}'.format(energy_list[i]),
                '{:.4f}'.format(relative_energy_list[i]),
                '{:}'.format(saddle_check_list[i])])

        if os.name == 'nt':
            encoding = 'utf-8-sig'
            newline = ''
        else:
            encoding = 'utf-8'
            newline = '\n'
        with result_csv_file.open(mode='w', encoding=encoding, newline=newline) as f:
            writer = csv.writer(f)
            writer.writerows(csv_data)

        # output xyz data
        title_list = []
        for i in range(len(coordinates_list)):
            title = ', '.join(['# {:}'.format(i + 1), '{:.8f}'.format(energy_list[i])])
            title_list.append(title)
        xyzutils.save_sequential_xyz_file(result_xyz_file, atom_symbols, coordinates_list, title_list)

        calc_success = True

    # common post-process
    finally:
        try:
            atoms.calc = None
        except:
            pass
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def _scan2d(input_xyz_file: Path,
            job_name: str,
            calculator: FAIRChemCalculator,
            params: UMAParams,
            scan1: UMAScan,
            scan2: UMAScan,
            constrains: List[UMAConstrain],
            keep_log: int = 0) -> None:


    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (job_name + '.xyz')
    stop_file: Path = input_xyz_file.parent / (job_name + config.STOP_FILE_SUFFIX)
    result_csv_file: Path = input_xyz_file.parent / (job_name + '.csv')
    workdir = Path(tempfile.mkdtemp(dir=input_xyz_file.parent, prefix=job_name + '_')).absolute()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / config.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    def _check_stop():
        if stop_file.exists():
            raise UMATerminationError('Stopped by user')

    # main process
    try:
        atoms = read(config.INIT_XYZ_FILE, index=0)
        atoms.info = {'charge': params.charge, 'spin': params.mult}
        atoms.calc = calculator
        scan_condition_list_2d = prepare_ase_constrain_list_scan2d(scan1, scan2, constrains)

        energy_list_2d = []  # List[array]
        coordinates_list_2d = []  # List[array]

        for i1, clist in enumerate(scan_condition_list_2d):
            energy_list = []
            coordinates_list = []
            for i2, scan_condition in enumerate(clist):
                # Set initial structure when going new line
                if i1 > 0 and i2 == 0:
                    atoms.set_positions(coordinates_list_2d[i1-1][0])
                try:
                    scan_dir = Path(f'./scan_{i1+1}_{i2+1}')
                    scan_dir.mkdir(parents=True)
                    trajfile = str(scan_dir / 'opt.traj')
                    trafxyzfile = str(scan_dir / 'opt.xyz')
                    logfile = str(scan_dir / 'opt.log')
                    atoms.set_constraint(scan_condition)
                    atoms.set_positions(atoms.get_positions(), apply_constraint=True)
                    opt = LBFGS(atoms, trajectory=trajfile, logfile=logfile)
                    opt.attach(_check_stop, interval=1)
                    opt.run(fmax=params.fmax, steps=params.max_cycles)
                    final_energy = atoms.get_potential_energy() / Hartree
                    if os.path.exists(trajfile):
                        try:
                            traj = read(f'{trajfile}@0:')
                            write(trafxyzfile, traj, format='xyz')
                        except Exception:
                            print('Trajectory conversion failed.')
                            traceback.print_exc()
                    coordinates_list.append(np.copy(atoms.positions))
                    energy_list.append(final_energy)
                except UMATerminationError:  # When user stop
                    try:
                        if os.path.exists(trajfile):
                            traj = read(f'{trajfile}@0:')
                            write(trafxyzfile, traj, format='xyz')
                        save_xyz_from_atoms(str(scan_dir / 'last.xyz'), atoms, title='USER_STOP')
                    finally:
                        calc_success = False
                    raise UMATerminationError('Stopped by user')
                except Exception as e:
                    raise RuntimeError(f'UMA optimization failed: {workdir}, scan {i1+1}, {i2+1}') from e

            energy_list_2d.append(np.array(energy_list, dtype=float))
            coordinates_list_2d.append(np.array(coordinates_list, dtype=float))

        atom_symbols = atoms.get_chemical_symbols()

        # summarized into one-dimensional data
        coordinates_list = []
        energy_list = []
        assigned_value_list1 = []
        assigned_value_list2 = []
        real_value_list1 = []
        real_value_list2 = []
        scan1_value_list = scan1.get_values()
        scan2_value_list = scan2.get_values()
        for k1 in range(scan1.num_step):
            for k2 in range(scan2.num_step):
                coordinates_list.append(coordinates_list_2d[k1][k2])
                energy_list.append(energy_list_2d[k1][k2])
                assigned_value_list1.append(scan1_value_list[k1])
                assigned_value_list2.append(scan2_value_list[k2])
                real_value_list1.append(scan1.calc_real_value(coordinates_list_2d[k1][k2]))
                real_value_list2.append(scan2.calc_real_value(coordinates_list_2d[k1][k2]))
        energy_list = np.array(energy_list, dtype=float)
        relative_energy_list = (energy_list - np.min(energy_list)) * config.HARTREE_TO_KCAL

        # check saddle. grad_tol > default value in config
        omp_num_threads = os.environ.get('OMP_NUM_THREADS', '')
        if not omp_num_threads:
            num_procs = None
        else:
            num_procs = int(omp_num_threads.split(',')[0])
        saddle_check_list = saddle.check_saddle_2d(relative_energy_list, scan1.num_step, scan2.num_step,
                                                   grad_tol=None, num_procs=num_procs)

        # output csv data
        csv_data = [['2d', str(scan1.num_step), str(scan2.num_step)],
                    ['#',
                     '[scan] ' + scan1.get_print_name(atom_symbols), '[real] ' + scan1.get_print_name(atom_symbols),
                     '[scan] ' + scan2.get_print_name(atom_symbols), '[real] ' + scan2.get_print_name(atom_symbols),
                     'E (au)', 'rel. E (kcal/mol)', 'saddle check']]
        for i in range(len(coordinates_list)):
            csv_data.append([
                '{:}'.format(i + 1),
                '{:.4f}'.format(assigned_value_list1[i]),
                '{:.4f}'.format(real_value_list1[i]),
                '{:.4f}'.format(assigned_value_list2[i]),
                '{:.4f}'.format(real_value_list2[i]),
                '{:.8f}'.format(energy_list[i]),
                '{:.4f}'.format(relative_energy_list[i]),
                '{:}'.format(saddle_check_list[i])])
        if os.name == 'nt':
            encoding = 'utf-8-sig'
            newline = ''
        else:
            encoding = 'utf-8'
            newline = '\n'
        with result_csv_file.open(mode='w', encoding=encoding, newline=newline) as f:
            writer = csv.writer(f)
            writer.writerows(csv_data)

        # output xyz data
        title_list = []
        for i in range(len(coordinates_list)):
            title = ', '.join(['# {:}'.format(i + 1), '{:.8f}'.format(energy_list[i])])
            title_list.append(title)
        xyzutils.save_sequential_xyz_file(result_xyz_file, atom_symbols, coordinates_list, title_list)

        calc_success = True

    # common post-process
    finally:
        try:
            atoms.calc = None
        except:
            pass
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def _scan_concerted(input_xyz_file: Path,
                    job_name: str,
                    calculator: FAIRChemCalculator,
                    params: UMAParams,
                    scans: List[UMAScan],
                    constrains: List[UMAConstrain],
                    keep_log: int = 0) -> None:
    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (job_name + '.xyz')
    stop_file: Path = input_xyz_file.parent / (job_name + config.STOP_FILE_SUFFIX)
    result_csv_file: Path = input_xyz_file.parent / (job_name + '.csv')
    workdir = Path(tempfile.mkdtemp(dir=input_xyz_file.parent, prefix=job_name + '_')).absolute()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / config.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    def _check_stop():
        if stop_file.exists():
            raise UMATerminationError('Stopped by user')

    # main process
    try:
        atoms = read(config.INIT_XYZ_FILE, index=0)
        atoms.info = {'charge': params.charge, 'spin': params.mult}
        atoms.calc = calculator
        scan_condition_list = prepare_ase_constrain_list_scan_concerted(scans, constrains)
        energy_list = []
        coordinates_list = []
        for n, scan_condition in enumerate(scan_condition_list):
            try:
                scan_dir = Path(f'./scan_{n+1}')
                scan_dir.mkdir(parents=True)
                trajfile = str(scan_dir / 'opt.traj')
                trafxyzfile = str(scan_dir / 'opt.xyz')
                logfile = str(scan_dir / 'opt.log')
                atoms.set_constraint(scan_condition)
                atoms.set_positions(atoms.get_positions(), apply_constraint=True)
                opt = LBFGS(atoms, trajectory=trajfile, logfile=logfile)
                opt.attach(_check_stop, interval=1)
                opt.run(fmax=params.fmax, steps=params.max_cycles)
                final_energy = atoms.get_potential_energy() / Hartree
                if os.path.exists(trajfile):
                    try:
                        traj = read(f'{trajfile}@0:')
                        write(trafxyzfile, traj, format='xyz')
                    except Exception:
                        print('Trajectory conversion failed.')
                        traceback.print_exc()
                coordinates_list.append(np.copy(atoms.positions))
                energy_list.append(final_energy)
            except UMATerminationError:  # When user stop
                try:
                    if os.path.exists(trajfile):
                        traj = read(f'{trajfile}@0:')
                        write(trafxyzfile, traj, format='xyz')
                    save_xyz_from_atoms(str(scan_dir / 'last.xyz'), atoms, title='USER_STOP')
                finally:
                    calc_success = False
                raise UMATerminationError('Stopped by user')
            except Exception as e:
                raise RuntimeError(f'UMA optimization failed: {workdir}, scan {n+1}') from e

        energy_list = np.array(energy_list, dtype=float)
        coordinates_list = np.array(coordinates_list, dtype=float)
        atom_symbols = atoms.get_chemical_symbols()

        relative_energy_list = (energy_list - np.min(energy_list)) * config.HARTREE_TO_KCAL
        assigned_values_list = [scan.get_values() for scan in scans]
        real_values_list = [[scan.calc_real_value(coord) for coord in coordinates_list] for scan in scans]
        saddle_check_list = saddle.check_saddle_1d(relative_energy_list)

        # output csv data
        header_line = ['#']
        for scan in scans:
            header_line.append('[scan] ' + scan.get_print_name(atom_symbols))
            header_line.append('[real] ' + scan.get_print_name(atom_symbols))
        header_line.extend(['E (au)', 'rel. E (kcal/mol)', 'saddle check'])
        csv_data = [['concerted', str(scans[0].num_step)], header_line]

        for i in range(len(coordinates_list)):
            data_line = ['{:}'.format(i + 1)]
            for k in range(len(scans)):
                data_line.append('{:.4f}'.format(assigned_values_list[k][i]))
                data_line.append('{:.4f}'.format(real_values_list[k][i]))
            data_line.extend(['{:.8f}'.format(energy_list[i]), '{:.4f}'.format(relative_energy_list[i])])
            data_line.append('{:}'.format(saddle_check_list[i]))
            csv_data.append(data_line)

        if os.name == 'nt':
            encoding = 'utf-8-sig'
            newline = ''
        else:
            encoding = 'utf-8'
            newline = '\n'
        with result_csv_file.open(mode='w', encoding=encoding, newline=newline) as f:
            writer = csv.writer(f)
            writer.writerows(csv_data)

        # output xyz data
        title_list = []
        for i in range(len(coordinates_list)):
            title = ', '.join(['# {:}'.format(i + 1), '{:.8f}'.format(energy_list[i])])
            title_list.append(title)
        xyzutils.save_sequential_xyz_file(result_xyz_file, atom_symbols, coordinates_list, title_list)

        calc_success = True

    # common post-process
    finally:
        try:
            atoms.calc = None
        except:
            pass
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def setenv_uma(num_threads: int = 1, memory_per_thread: Optional[str] = None) -> None:
    if not IS_UMA_VALID:
        raise RuntimeError('Import of ASE/torch/fairchem-core failed.')
    assert int(num_threads) > 0
    num_threads = int(num_threads)

    if memory_per_thread is None:
        memory_per_thread = '500M'
    memory_per_thread = str(memory_per_thread)
    if memory_per_thread.lower().endswith('b'):
        memory_per_thread = memory_per_thread[:-1]

    torch.set_num_threads(num_threads)
    os.environ['OMP_NUM_THREADS'] = str(num_threads)
    os.environ['MKL_NUM_THREADS'] = str(num_threads)
    os.environ['OMP_STACKSIZE'] = memory_per_thread
    global CHECK_SETENV_UMA
    CHECK_SETENV_UMA = True


def save_xyz_from_atoms(file, atoms: Atoms, title: str):
    with open(file, 'w') as f:
        f.write(f'{len(atoms)}\n')
        f.write(f'{title}\n')
        for atom in atoms:
            symbol = atom.symbol
            x, y, z = atom.position
            f.write(f'{symbol:>2} {x:>20.12f} {y:>20.12f} {z:>20.12f}\n')