import os
import shutil
import time
import csv
import subprocess as sp
from pathlib import Path

import numpy as np

import const
import config
import utils
import xyzutils

# Global variable to check whether setenv run or not
CHECK_SETENV = False


class XTBTerminationError(Exception):
    pass


class XTBParams:
    def __init__(self, method='gfn2', charge=0, uhf=0, solvation=None, solvent=''):
        # check parameters
        try:
            method = method.lower()
            assert method in ['gfn1', 'gfn2', 'gfnff']
            assert uhf >= 0
            if solvation is not None:
                solvation = solvation.lower()
                assert solvation in ['alpb', 'gbsa']
                assert len(solvent) != ''
        except AssertionError as e:
            raise ValueError('Given XTB parameters are not valid. ' + str(e.args))

        self.method = method
        self.charge = charge
        self.uhf = uhf
        self.solation = solvation
        self.solvent = solvent

    @property
    def args(self):
        _args = ['--' + self.method, '--chrg', str(self.charge), '--uhf', str(self.uhf)]
        if self.solation is not None:
            _args.extend(['--' + self.solation, self.solvent])
        return _args


class XTBConstrain:

    def __init__(self, constrain_type: str, atom_indices, value=None):

        assert constrain_type in ['atoms', 'distance', 'angle', 'dihedral']

        if constrain_type == 'atoms':
            assert value is None
            self.value = None
        if constrain_type != 'atoms':
            if value is None:
                self.value = 'auto'
            elif str(value).lower().strip() in ['', 'a', 'auto']:
                self.value = 'auto'
            else:
                self.value = str(value).strip()

        if constrain_type == 'distance':
            assert len(atom_indices) == 2
        if constrain_type == 'angle':
            assert len(atom_indices) == 3
        if constrain_type == 'dihedral':
            assert len(atom_indices) == 4

        self.constrain_type = constrain_type
        self.atom_indices = atom_indices

    def get_constrain_string(self):
        if self.constrain_type == 'atoms':
            return '  atoms: {}\n'.format(utils.atom_indices_to_string(self.atom_indices))
        if self.constrain_type == 'distance':
            return '  distance: {}, {}, {}\n'.format(self.atom_indices[0] + 1,
                                                     self.atom_indices[1] + 1,
                                                     self.value)
        if self.constrain_type == 'angle':
            return '  angle: {}, {}, {}, {}\n'.format(self.atom_indices[0] + 1,
                                                      self.atom_indices[1] + 1,
                                                      self.atom_indices[2] + 1,
                                                      self.value)
        if self.constrain_type == 'dihedral':
            return '  dihedral: {}, {}, {}, {}, {}\n'.format(self.atom_indices[0] + 1,
                                                             self.atom_indices[1] + 1,
                                                             self.atom_indices[2] + 1,
                                                             self.atom_indices[3] + 1,
                                                             self.value)

    def __str__(self):
        return self.get_constrain_string()

    def get_print_name(self, atoms=None):
        pass


class XTBScan:
    def __init__(self, scan_type: str, atom_indices, start, end, num_step):

        assert int(num_step) >= 2
        assert scan_type in ['distance', 'angle', 'dihedral']

        if scan_type == 'distance':
            assert len(atom_indices) == 2
        if scan_type == 'angle':
            assert len(atom_indices) == 3
        if scan_type == 'dihedral':
            assert len(atom_indices) == 4

        self.scan_type = scan_type
        self.atom_indices = atom_indices
        self.start = start
        self.end = end
        self.num_step = int(num_step)

    def get_scan_string(self):
        scan_string = '  {}: '.format(self.scan_type)
        scan_string += ','.join([str(v + 1) for v in self.atom_indices]) + ','
        scan_string += '{start}; {start},{end},{num_step}\n'.format(start=self.start, end=self.end,
                                                                    num_step=self.num_step)
        return scan_string

    def get_values(self):
        return np.linspace(float(self.start), float(self.end), int(self.num_step), endpoint=True, dtype=const.FLOAT)

    def calc_real_value(self, coordinates: np.ndarray) -> float:
        if self.scan_type == 'distance':
            return xyzutils.calc_distance(coordinates, self.atom_indices)
        elif self.scan_type == 'anlge':
            return xyzutils.calc_angle(coordinates, self.atom_indices)
        elif self.scan_type == 'dihedral':
            return xyzutils.calc_dihedral(coordinates, self.atom_indices)
        else:
            raise RuntimeError('unknown type of scan.')

    def get_print_name(self, atoms=None):
        if atoms is None:
            atoms = np.full('Atom', max(self.atom_indices) + 2)
        name = self.scan_type + ' ' + '-'.join([atoms[i].capitalize() + str(i + 1) for i in self.atom_indices])
        unit = ' (ang.)' if self.scan_type == 'distance' else ' (deg.)'
        return name + unit

    def __str__(self):
        return self.get_scan_string()


def xtbscan(input_xyz_file, xtb_params: XTBParams, scans, constrains,
            force_constant='1.0', concerted=False, keep_log=0):
    # initial check
    if not CHECK_SETENV:
        setenv(num_threads=1, memory_per_thread='500M')
    input_xyz_file = Path(input_xyz_file).resolve()
    if not input_xyz_file.exists():
        raise FileNotFoundError(str(input_xyz_file) + ' not found.')

    # run
    if not scans:
        _opt(input_xyz_file, xtb_params, constrains, force_constant, keep_log)
    elif len(scans) == 1:
        _scan1d(input_xyz_file, xtb_params, scans[0], constrains, force_constant, keep_log)
    elif concerted:
        _scan_concerted(input_xyz_file, xtb_params, scans, constrains, force_constant, keep_log)
    elif len(scans) == 2:
        _scan2d(input_xyz_file, xtb_params, scans[0], scans[1], constrains, force_constant, keep_log)
    else:
        raise ValueError('Only Opt, 1D scan, 2D scan, or concerted scan is available.')


def _opt(input_xyz_file, xtb_params: XTBParams, constrains, force_constant='1.0', keep_log=0):
    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.OPT_XYZ_SUFFIX)
    stop_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.STOP_FILE_SUFFIX)
    workdir = input_xyz_file.parent / (const.WORKDIR_PREFIX + input_xyz_file.stem + const.WORKDIR_SUFFIX)
    workdir = utils.get_unused_directory_path(workdir)
    workdir.mkdir()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / const.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    # main process
    try:
        if len(constrains) > 0:
            _save_input_file(const.INPUT_FILE, [], constrains, force_constant, concerted=False)
            command = [config.XTB_BIN, const.INIT_XYZ_FILE, '--opt', '--input', const.INPUT_FILE]
        else:
            command = [config.XTB_BIN, const.INIT_XYZ_FILE, '--opt']
        command.extend(xtb_params.args)

        with open(const.XTB_LOG_FILE, 'w', encoding='utf-8') as f:
            proc = utils.popen_bg(command, universal_newlines=True, encoding='utf-8', stdout=f, stderr=sp.STDOUT)
            while True:
                time.sleep(const.STOP_CHECK_INTERVAL)
                if stop_file.exists():
                    proc.terminate()
                    raise XTBTerminationError('Stopped by user')
                if proc.poll() is not None:
                    break
        with open(const.XTB_LOG_FILE, 'r', encoding='utf-8') as f:
            if 'normal termination of xtb' not in f.read():
                raise RuntimeError('xtb optimization failed in {:}'.format(workdir))

        shutil.copy(const.XTB_OPT_FILE, result_xyz_file)
        calc_success = True

    # common post-process
    finally:
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def _scan1d(input_xyz_file, xtb_params: XTBParams, scan: XTBScan, constrains, force_constant='1.0', keep_log=0):
    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.RESULT_XYZ_SUFFIX)
    result_csv_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.RESULT_CSV_SUFFIX)
    stop_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.STOP_FILE_SUFFIX)
    workdir = input_xyz_file.parent / (const.WORKDIR_PREFIX + input_xyz_file.stem + const.WORKDIR_SUFFIX)
    workdir = utils.get_unused_directory_path(workdir)
    workdir.mkdir()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / const.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    # main process
    try:
        _save_input_file(const.INPUT_FILE, [scan], constrains, force_constant, concerted=False)
        command = [config.XTB_BIN, const.INIT_XYZ_FILE, '--opt', '--input', const.INPUT_FILE]
        command.extend(xtb_params.args)

        with open(const.XTB_LOG_FILE, 'w', encoding='utf-8') as f:
            proc = utils.popen_bg(command, universal_newlines=True, encoding='utf-8', stdout=f, stderr=sp.STDOUT)
            while True:
                time.sleep(const.STOP_CHECK_INTERVAL)
                if stop_file.exists():
                    proc.terminate()
                    raise XTBTerminationError('Stopped by user')
                if proc.poll() is not None:
                    break
        with open(const.XTB_LOG_FILE, 'r', encoding='utf-8') as f:
            if 'normal termination of xtb' not in f.read():
                raise RuntimeError('xtb scan failed in {:}'.format(workdir))

        atoms, coordinates_list, energy_list = xyzutils.read_xtbscan_file(const.XTB_SCAN_FILE)
        relative_energy_list = (energy_list - np.min(energy_list)) * const.HARTREE_TO_KCAL
        assigned_value_list = scan.get_values()
        real_value_list = [scan.calc_real_value(coord) for coord in coordinates_list]
        saddle_check_list = _check_saddle_1d(energy_list)

        # output csv data
        csv_data = [['1d', str(scan.num_step)],
                    ['#', '[scan] ' + scan.get_print_name(atoms), '[real] ' + scan.get_print_name(atoms),
                     'E (au)', 'rel. E (kcal/mol)', 'saddle check']]

        for i in range(len(coordinates_list)):
            csv_data.append([
                '{:}'.format(i + 1),
                '{:.4f}'.format(assigned_value_list[i]),
                '{:.4f}'.format(real_value_list[i]),
                '{:.8f}'.format(energy_list[i]),
                '{:.4f}'.format(relative_energy_list[i]),
                '{:}'.format(saddle_check_list[i])])
        with result_csv_file.open(mode='w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(csv_data)

        # output xyz data
        title_list = []
        for i in range(len(coordinates_list)):
            title = ', '.join(['# {:}'.format(i + 1), '{:.8f}'.format(energy_list[i])])
            title_list.append(title)
        xyzutils.save_sequential_xyz_file(result_xyz_file, atoms, coordinates_list, title_list)

        calc_success = True

    # common post-process
    finally:
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def _scan2d(input_xyz_file, xtb_params: XTBParams, scan1: XTBScan, scan2: XTBScan, constrains,
            force_constant='1.0', keep_log=0):
    if constrains is None:
        constrains = []
    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.RESULT_XYZ_SUFFIX)
    result_csv_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.RESULT_CSV_SUFFIX)
    stop_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.STOP_FILE_SUFFIX)
    workdir = input_xyz_file.parent / (const.WORKDIR_PREFIX + input_xyz_file.stem + const.WORKDIR_SUFFIX)
    workdir = utils.get_unused_directory_path(workdir)
    workdir.mkdir()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / const.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    # main process
    try:
        # xtb command (common for all scans)
        command = ['xtb', const.INIT_XYZ_FILE, '--opt', '--input', const.INPUT_FILE]
        command.extend(xtb_params.args)

        # first dimension scan
        internal_workdir = workdir / 'scan1'
        internal_workdir.mkdir()
        os.chdir(internal_workdir)
        shutil.copy(workdir / const.INIT_XYZ_FILE, internal_workdir)
        constrain2 = XTBConstrain(constrain_type=scan2.scan_type, atom_indices=scan2.atom_indices, value='auto')
        _save_input_file(const.INPUT_FILE, [scan1], constrains + [constrain2], force_constant)

        with open(const.XTB_LOG_FILE, 'w', encoding='utf-8') as f:
            proc = utils.popen_bg(command, universal_newlines=True, encoding='utf-8', stdout=f, stderr=sp.STDOUT)
            while True:
                time.sleep(const.STOP_CHECK_INTERVAL)
                if stop_file.exists():
                    proc.terminate()
                    raise XTBTerminationError('Stopped by user')
                if proc.poll() is not None:
                    break
        with open(const.XTB_LOG_FILE, 'r', encoding='utf-8') as f:
            if 'normal termination of xtb' not in f.read():
                raise RuntimeError('xtb scan failed in {:}'.format(workdir))

        atoms, first_scanned_coordinates_list, _ = xyzutils.read_xtbscan_file(const.XTB_SCAN_FILE)
        os.chdir(workdir)

        # second dimension scan
        coordinates_list_2d = []
        energy_list_2d = []
        constrain1 = XTBConstrain(constrain_type=scan1.scan_type, atom_indices=scan1.atom_indices, value='auto')
        for i1 in range(scan1.num_step):
            internal_workdir = workdir / 'scan2_{:}'.format(i1 + 1)
            internal_workdir.mkdir()
            os.chdir(internal_workdir)
            xyzutils.save_xyz_file(const.INIT_XYZ_FILE, atoms, first_scanned_coordinates_list[i1], 'init')
            _save_input_file(const.INPUT_FILE, [scan2], constrains + [constrain1], force_constant)

            with open(const.XTB_LOG_FILE, 'w', encoding='utf-8') as f:
                proc = utils.popen_bg(command, universal_newlines=True, encoding='utf-8', stdout=f, stderr=sp.STDOUT)
                while True:
                    time.sleep(const.STOP_CHECK_INTERVAL)
                    if stop_file.exists():
                        proc.terminate()
                        raise XTBTerminationError('Stopped by user')
                    if proc.poll() is not None:
                        break
            with open(const.XTB_LOG_FILE, 'r', encoding='utf-8') as f:
                if 'normal termination of xtb' not in f.read():
                    raise RuntimeError('xtb scan failed in {:}'.format(workdir))

            _, coordinates_list, energy_list = xyzutils.read_xtbscan_file(const.XTB_SCAN_FILE)
            coordinates_list_2d.append(coordinates_list)
            energy_list_2d.append(energy_list)
            os.chdir(workdir)

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
        energy_list = np.array(energy_list, dtype=const.FLOAT)
        relative_energy_list = (energy_list - np.min(energy_list)) * const.HARTREE_TO_KCAL
        saddle_check_list = _check_saddle_2d(relative_energy_list, scan1.num_step, scan2.num_step)

        # output csv data
        csv_data = [['2d', str(scan1.num_step), str(scan2.num_step)],
                    ['#',
                     '[scan] ' + scan1.get_print_name(atoms), '[real] ' + scan1.get_print_name(atoms),
                     '[scan] ' + scan2.get_print_name(atoms), '[real] ' + scan2.get_print_name(atoms),
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
        with result_csv_file.open(mode='w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(csv_data)

        # output xyz data
        title_list = []
        for i in range(len(coordinates_list)):
            title = ', '.join(['# {:}'.format(i + 1), '{:.8f}'.format(energy_list[i])])
            title_list.append(title)
        xyzutils.save_sequential_xyz_file(result_xyz_file, atoms, coordinates_list, title_list)

        calc_success = True

    # common post-process
    finally:
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def _scan_concerted(input_xyz_file, xtb_params: XTBParams, scans, constrains, force_constant='1.0', keep_log=0):
    # common pre-process
    result_xyz_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.RESULT_XYZ_SUFFIX)
    result_csv_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.RESULT_CSV_SUFFIX)
    stop_file: Path = input_xyz_file.parent / (input_xyz_file.stem + const.STOP_FILE_SUFFIX)
    workdir = input_xyz_file.parent / (const.WORKDIR_PREFIX + input_xyz_file.stem + const.WORKDIR_SUFFIX)
    workdir = utils.get_unused_directory_path(workdir)
    workdir.mkdir()
    prevdir = os.getcwd()
    shutil.copy(input_xyz_file, workdir / const.INIT_XYZ_FILE)
    os.chdir(str(workdir))
    calc_success = False

    # main process
    try:
        _save_input_file(const.INPUT_FILE, scans, constrains, force_constant, concerted=True)
        command = ['xtb', const.INIT_XYZ_FILE, '--opt', '--input', const.INPUT_FILE]
        command.extend(xtb_params.args)

        with open(const.XTB_LOG_FILE, 'w', encoding='utf-8') as f:
            proc = utils.popen_bg(command, universal_newlines=True, encoding='utf-8', stdout=f, stderr=sp.STDOUT)
            while True:
                time.sleep(const.STOP_CHECK_INTERVAL)
                if stop_file.exists():
                    proc.terminate()
                    raise XTBTerminationError('Stopped by user')
                if proc.poll() is not None:
                    break
        with open(const.XTB_LOG_FILE, 'r', encoding='utf-8') as f:
            if 'normal termination of xtb' not in f.read():
                raise RuntimeError('xtb scan failed in {:}'.format(workdir))

        atoms, coordinates_list, energy_list = xyzutils.read_xtbscan_file(const.XTB_SCAN_FILE)
        relative_energy_list = (energy_list - np.min(energy_list)) * const.HARTREE_TO_KCAL
        assigned_values_list = [scan.get_values() for scan in scans]
        real_values_list = [[scan.calc_real_value(coord) for coord in coordinates_list] for scan in scans]
        saddle_check_list = _check_saddle_1d(relative_energy_list)

        # output csv data
        header_line = ['#']
        for scan in scans:
            header_line.append('[scan] ' + scan.get_print_name(atoms))
            header_line.append('[real] ' + scan.get_print_name(atoms))
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

        with result_csv_file.open(mode='w', newline='\n') as f:
            writer = csv.writer(f)
            writer.writerows(csv_data)

        # output xyz data
        title_list = []
        for i in range(len(coordinates_list)):
            title = ', '.join(['# {:}'.format(i + 1), '{:.8f}'.format(energy_list[i])])
            title_list.append(title)
        xyzutils.save_sequential_xyz_file(result_xyz_file, atoms, coordinates_list, title_list)

        calc_success = True

    # common post-process
    finally:
        os.chdir(prevdir)
        if (keep_log == 0) or (keep_log == 1 and calc_success):
            shutil.rmtree(workdir, ignore_errors=True)


def setenv(num_threads=1, memory_per_thread=None):
    assert int(num_threads) > 0
    num_threads = str(int(num_threads))

    if memory_per_thread is None:
        memory_per_thread = '500M'
    memory_per_thread = str(memory_per_thread)
    if memory_per_thread.lower().endswith('b'):
        memory_per_thread = memory_per_thread[:-1]

    os.environ['XTBPATH'] = str(config.XTB_PARAM_DIR)
    os.environ['OMP_NUM_THREADS'] = num_threads + ',1'
    os.environ['OMP_STACKSIZE'] = memory_per_thread
    os.environ['MKL_NUM_THREADS'] = num_threads

    if config.XTB_BIN_DIR not in os.environ['PATH']:
        os.environ['PATH'] = '{:};{:}'.format(config.XTB_BIN_DIR, os.environ['PATH'])
    if config.XTB_OTHER_LIB_DIR is not None:
        if config.XTB_PARAM_DIR not in os.environ['PATH']:
            os.environ['PATH'] = '{:};{:}'.format(config.XTB_OTHER_LIB_DIR, os.environ['PATH'])

    global CHECK_SETENV
    CHECK_SETENV = True


def _save_input_file(file, scans, constrains, force_constant='1.0', concerted=False) -> None:
    data = ['$constrain\n', '  force constant={}\n'.format(force_constant)]

    for constrain in constrains:
        data.append(str(constrain))

    if len(scans) > 0:
        data.append('$scan\n')

        # concerted mode. check all the scan steps are same
        if concerted:
            num_scan = int(scans[0].num_step)
            for scan in scans:
                if scan.num_step != num_scan:
                    raise ValueError('For concerted scan, all scan steps must be same.')
            data.append('  mode=concerted\n')

        for scan in scans:
            data.append(str(scan))

    data.append('$end')

    file = Path(file)
    with file.open(mode='w') as f:
        f.writelines(data)


def _check_saddle_1d(energies: np.ndarray) -> np.ndarray:
    saddle_check_list = np.full_like(energies, 0, dtype='int32')
    for i in range(1, len(energies) - 1):
        if energies[i - 1] < energies[i] and energies[i + 1] < energies[i]:
            saddle_check_list[i] = 1
    return saddle_check_list


def _check_saddle_2d(energies: np.ndarray, num_dim1, num_dim2) -> np.ndarray:
    _energies = energies.reshape((num_dim1, num_dim2))
    _saddle_check_list = np.full_like(_energies, 0, dtype='int32')
    for x in range(1, num_dim1 - 1):
        for y in range(1, num_dim2 - 1):
            num_min = 0
            num_max = 0
            num_slope = 0

            for v1, v2, v3 in [(_energies[x - 1, y - 1], _energies[x, y], _energies[x + 1, y + 1]),
                               (_energies[x - 1, y], _energies[x, y], _energies[x + 1, y]),
                               (_energies[x - 1, y + 1], _energies[x, y], _energies[x + 1, y - 1]),
                               (_energies[x, y - 1], _energies[x, y], _energies[x, y + 1])]:
                if v1 < v2 and v3 < v2:
                    num_max += 1
                elif v1 > v2 and v3 > v2:
                    num_min += 1
                else:
                    num_slope += 1

            if num_min > 0 and num_max > 0:
                _saddle_check_list[x, y] = 1

            # too strict?
            """
            if num_slope > 0:
                continue
            elif num_min > 0 and num_max > 0:
                _saddle_check_list[x,y] = 1
            """

    return _saddle_check_list.flatten()


# Additional functions for geomTric Engines #


def xtb_energy_and_gradient(atoms, coordinates, xtb_params: XTBParams, workdir: Path):
    # initial check
    if not CHECK_SETENV:
        setenv(num_threads=1, memory_per_thread='500M')

    if workdir.exists():
        raise RuntimeError('Indicated XTB working directory already exists.')

    workdir.mkdir()
    prevdir = os.getcwd()

    try:
        os.chdir(str(workdir))
        xyzutils.save_xyz_file(const.INIT_XYZ_FILE, atoms, coordinates, 'for energy and gradient')
        command = [config.XTB_BIN, const.INIT_XYZ_FILE, '--grad']
        command.extend(xtb_params.args)
        with open(const.XTB_LOG_FILE, 'w', encoding='utf-8') as f:
            proc = utils.popen_bg(command, universal_newlines=True, encoding='utf-8', stdout=f, stderr=sp.STDOUT)
            proc.wait()
        with open(const.XTB_LOG_FILE, 'r', encoding='utf-8') as f:
            if 'normal termination of xtb' not in f.read():
                raise RuntimeError('xtb optimization failed in {:}'.format(workdir))
        energy = _read_xtb_energy(const.XTB_ENERGY_FILE)
        gradient = _read_xtb_gradient(const.XTB_GRADIENT_FILE, len(atoms))

    except:
        raise

    else:
        return energy, gradient

    finally:
        os.chdir(prevdir)
        try:
            shutil.rmtree(workdir)
        except:
            pass


def xtb_hessian(atoms, coordinates, xtb_params: XTBParams, workdir: Path):
    # initial check
    if not CHECK_SETENV:
        setenv(num_threads=1, memory_per_thread='500M')

    if workdir.exists():
        raise RuntimeError('Indicated XTB working directory already exists.')

    workdir.mkdir()
    prevdir = os.getcwd()

    try:
        os.chdir(str(workdir))
        xyzutils.save_xyz_file(const.INIT_XYZ_FILE, atoms, coordinates, 'for hessian')
        command = [config.XTB_BIN, const.INIT_XYZ_FILE, '--hess']
        command.extend(xtb_params.args)
        with open(const.XTB_LOG_FILE, 'w', encoding='utf-8') as f:
            proc = utils.popen_bg(command, universal_newlines=True, encoding='utf-8', stdout=f, stderr=sp.STDOUT)
            proc.wait()
        with open(const.XTB_LOG_FILE, 'r', encoding='utf-8') as f:
            if 'normal termination of xtb' not in f.read():
                raise RuntimeError('xtb optimization failed in {:}'.format(workdir))
        hess = _read_xtb_hessian(const.XTB_HESSIAN_FILE, len(atoms))

    except:
        raise

    else:
        return hess

    finally:
        os.chdir(prevdir)
        try:
            shutil.rmtree(workdir)
        except:
            pass


def _read_xtb_energy(energy_file) -> float:
    energy_file = Path(energy_file)
    with energy_file.open(mode='r') as f:
        energy_data = f.readlines()
    assert energy_data[0].strip().startswith('$energy')
    return float(energy_data[1].strip().split()[1])


def _read_xtb_gradient(gradient_file, num_atom: int) -> np.ndarray:
    gradient_file = Path(gradient_file)
    with gradient_file.open(mode='r') as f:
        gradient_data = f.readlines()
    assert gradient_data[0].strip().startswith('$grad')
    gradient = []
    for line in gradient_data[2+num_atom:2+2*num_atom]:
        gradient.append(line.strip().split())

    return np.array(gradient, dtype=const.FLOAT)


def _read_xtb_hessian(hessian_file, num_atom: int) -> np.ndarray:
    hessian_file = Path(hessian_file)
    with hessian_file.open(mode='r') as f:
        hessian_data = f.readlines()
    assert hessian_data[0].strip().startswith('$hessian')

    hess_1d = []
    for line in hessian_data[1:]:
        if line.strip() == '':
            continue
        if line.strip().startswith('$end'):
            break
        hess_1d.extend(line.strip().split())

    assert len(hess_1d) == (num_atom * 3) ** 2

    return np.array(hess_1d, dtype=const.FLOAT).reshape((3*num_atom, 3*num_atom))
