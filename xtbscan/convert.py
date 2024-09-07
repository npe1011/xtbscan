from pathlib import Path
import numpy as np

try:
    import cclib
except ImportError:
    CCLIB = False
else:
    CCLIB = True

from config import ATOM_LIST, FLOAT
from xtbscan import xyzutils


def gaussian_input_to_xyz(file: Path) -> Path:
    """
    convert gaussian input to xyz file.
    return Path to the generated xyz file.
    """
    input_data = []
    structure_data = []

    file = Path(file).absolute()

    with file.open(mode='r') as f:

        # preprocess: remove comments and continuous blank lines
        chk_void = False
        line = f.readline()
        while line:
            line = line.lstrip()
            if not line:
                if not chk_void:
                    input_data.append("\n")
                    chk_void = True
                else:
                    pass
            else:
                if line[0] != "!":
                    chk_void = False
                    input_data.append(line)
            line = f.readline()

        # read data
        state = 0
        for data_line in input_data:
            if state == 0:  # before route section
                if data_line[0] == '#':
                    state = 1
            elif state == 1:  # during route section
                if data_line == '\n':  # blank line > title section in the next line
                    state = 2
            elif state == 2:  # title section
                if data_line == '\n':  # blank line > charge, multi, structure from the next line
                    state = 3
            else:
                if data_line != '\n':
                    if data_line[0:2].upper() != 'LP':  # ignore lone pairs
                        structure_data.append(data_line)
                else:  # blank line > end structure.
                    break

    output_file = file.parent / (file.stem + '.xyz')
    with output_file.open(mode='w') as f:
        f.write('{:}\n'.format(len(structure_data)-1))
        f.write('structure read from gaussian input file\n')
        f.writelines(structure_data[1:])

    return output_file


def gaussian_log_to_xyz(file: Path) -> Path:
    """
    convert gaussian log to xyz file (read final structure)
    return Path to the generated xyz file.
    """

    with Path(file).open(mode='r') as f:
        log_data = f.readlines()

        # num_coord、num_charge_multi
        for (num, line) in enumerate(log_data):
            if 'Input orientation:' in line:
                num_coord = num
            if "Standard orientation:" in line:
                num_coord = num
            # elif "Multiplicity =" in line:
                # line_charge_multi = line

        # charge = int(line_charge_multi.strip().split()[2])
        # multi = int(line_charge_multi.strip().split()[5])

        i = num_coord + 5  # start line for structure information

        atoms = []
        coordinates = []
        while log_data[i] and ('------' not in log_data[i]):  # blank line or ---- line means end
            atom_coord = log_data[i].strip().split()
            atom_label = ATOM_LIST[int(atom_coord[1])]
            coord_x = atom_coord[3]
            coord_y = atom_coord[4]
            coord_z = atom_coord[5]
            atoms.append(atom_label)
            coordinates.append([float(coord_x), float(coord_y), float(coord_z)])
            i = i + 1

    atoms = np.array(atoms)
    coordinates = np.array(coordinates, dtype=FLOAT)
    output_file = file.parent / (file.stem + '.xyz')
    xyzutils.save_xyz_file(output_file, atoms, coordinates, 'structure read from gaussian log file')

    return output_file


def other_to_xyz(file: Path) -> Path:
    if not CCLIB:
        raise RuntimeError('cclib is required to parse other than Gaussian files.')
    data = cclib.io.ccread(str(file))

    coordinates = data.atomcoords[-1,::]
    atoms = np.array([ATOM_LIST[n] for n in data.atomnos])
    output_file = file.parent / (file.stem + '.xyz')
    xyzutils.save_xyz_file(output_file, atoms, coordinates, 'structure read from log file via cclib')
    return output_file
