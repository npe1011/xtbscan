from typing import List, Tuple, Union
from pathlib import Path

import numpy as np

from config import XYZ_FORMAT


def save_xyz_file(xyz_file: Union[str, Path], atoms: np.ndarray, coordinates: np.ndarray, title: str = '') -> None:
    xyz_file = Path(xyz_file)
    data = [
        str(len(atoms)) + '\n',
        title.rstrip() + '\n'
    ]
    for i in range(len(atoms)):
        data.append(XYZ_FORMAT.format(atoms[i], coordinates[i][0], coordinates[i][1], coordinates[i][2]))
    with xyz_file.open(mode='w', encoding='utf-8', newline='\n') as f:
        f.writelines(data)


def read_single_xyz_file(xyz_file: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray]:
    """
    :return: tuple(atoms: np.ndarray, coordinates: np.ndarray)
    """
    xyz_file = Path(xyz_file)
    with xyz_file.open(mode='r') as f:
        data = f.readlines()
    atoms = []
    coordinates = []
    num_atoms = int(data[0].strip())
    for line in data[2:2+num_atoms]:
        temp = line.strip().split()
        atom = temp[0].capitalize()
        x = float(temp[1])
        y = float(temp[2])
        z = float(temp[3])
        atoms.append(atom)
        coordinates.append([x, y, z])
    return np.array(atoms), np.array(coordinates, dtype=float)


def read_sequential_xyz_file(file: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    return: tuple(atoms, coordinates_list, title_list)
    atoms: np.ndarray
    coordinates_list np.ndarray 3d (num_conf x num_atoms x 3)
    title_list list of str
    """

    with Path(file).open() as f:
        data = f.readlines()

    atoms = []
    coordinates_list = []
    title_list = []

    num_atoms = int(data[0].strip())

    # check if validity
    start_lines = []
    for i, line in enumerate(data):
        if i % (num_atoms + 2) == 0:
            if line.strip() == str(num_atoms):
                start_lines.append(i)
            else:
                if line.strip() == '':
                    break
                else:
                    raise ValueError(str(file) + 'is not a valid xyz file.')

    # Read each coordinates and title
    for start_line in start_lines:
        title_list.append(data[start_line+1].strip())
        coordinates = []
        for line in data[start_line+2:start_line+num_atoms+2]:
            temp = line.strip().split()
            x = float(temp[1])
            y = float(temp[2])
            z = float(temp[3])
            coordinates.append([x, y, z])
        coordinates_list.append(np.array(coordinates, dtype=float))

    # Read atoms
    for line in data[2:2+num_atoms]:
        atoms.append(line.strip().split()[0].capitalize())

    return np.array(atoms), np.array(coordinates_list, dtype=float), title_list


def read_xtbscan_file(file: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    """
    return: tuple(atoms, coordinates_list, energy_list)
    atoms: np.ndarray
    coordinates_list np.ndarray 3d (num_conf x num_atoms x 3)
    energy_list np.ndarray 1d (num_conf)
    """

    atoms, coordinates_list, title_list = read_sequential_xyz_file(file)

    energy_list = []
    for title in title_list:
        title: str = title.strip().lower()
        if title.startswith('energy'):
            # case
            # energy: -167.316733480097 xtbscan: 6.4.1 (unknown)
            energy_list.append(float(title.split(':', maxsplit=1)[1].strip().split()[0].strip()))
        elif title.startswith('scf done'):
            # case
            # SCF done      -7.33636977
            energy_list.append(title.split('done', maxsplit=1)[1].strip().split()[0].strip())
        else:
            raise ValueError(str(file) + 'is not a valid xtbscan file. Energy value cannot be read.')

    return atoms, coordinates_list, np.array(energy_list, dtype=float)


def save_sequential_xyz_file(xyz_file: Union[str, Path],
                             atoms: np.ndarray,
                             coordinates_list: Union[np.ndarray, List[np.ndarray]],
                             title_list: List[str]) -> None:
    """
    atoms: np.ndarray
    coordinates_list np.ndarray 3d (num_conf x num_atoms x 3)
    title_list: str list
    """
    data = []
    num_atoms = len(atoms)
    for (coordinates, title) in zip(coordinates_list, title_list):
        data.append('{:5d}\n'.format(num_atoms))
        data.append(title.rstrip() + '\n')
        for i in range(len(atoms)):
            data.append(XYZ_FORMAT.format(atoms[i], coordinates[i][0], coordinates[i][1], coordinates[i][2]))

    with Path(xyz_file).open(mode='w', encoding='utf-8', newline='\n') as f:
        f.writelines(data)


def get_xyz_string(atoms: Union[List[str], np.ndarray], coordinates: np.ndarray) -> str:
    data = []
    for i in range(len(atoms)):
        data.append(XYZ_FORMAT.format(atoms[i], coordinates[i][0], coordinates[i][1], coordinates[i][2]))
    return ''.join(data)


def get_num_structures(xyz_file: Union[str, Path]) -> int:
    _, _, title_list = read_sequential_xyz_file(xyz_file)
    return len(title_list)


def calc_distance(coordinates: np.ndarray, atom_indices: Union[np.ndarray, List[int]], string=False) -> Union[float, str]:
    assert len(atom_indices) == 2
    v = float(np.linalg.norm(coordinates[atom_indices[0]] - coordinates[atom_indices[1]]))
    if not string:
        return v
    else:
        return '{:.4f}'.format(v)


def calc_angle(coordinates: np.ndarray, atom_indices: Union[np.ndarray, List[int]], string=False) -> Union[float, str]:
    assert len(atom_indices) == 3
    target_coordinates = coordinates[atom_indices]

    a = np.linalg.norm(target_coordinates[2] - target_coordinates[1])
    b = np.linalg.norm(target_coordinates[0] - target_coordinates[1])
    c = np.linalg.norm(target_coordinates[2] - target_coordinates[0])
    cos = (a*a + b*b - c*c) / (2.0 * a * b)
    v = float(np.degrees(np.arccos(cos)))

    if not string:
        return v
    else:
        return '{:.1f}'.format(v)


def calc_dihedral(coordinates: np.ndarray, atom_indices: Union[np.ndarray, List[int]],
                  string=False) -> Union[float, str]:
    """
    return dihedral angle (degrees) of four atoms specified with atom_indices
    :param coordinates: coordinates
    :param atom_indices:  four atom numbers (0 start)
    :param string:  return formatted string if True, else return in float
    stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    wikipedia formula
    """

    assert(len(atom_indices)) == 4
    target_coordinates = coordinates[atom_indices]

    b0 = -1.0 * (target_coordinates[1] - target_coordinates[0])
    b1 = target_coordinates[2] - target_coordinates[1]
    b2 = target_coordinates[3] - target_coordinates[2]

    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)

    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)
    y = np.dot(b0xb1_x_b1xb2, b1) * (1.0 / np.linalg.norm(b1))

    x = np.dot(b0xb1, b1xb2)

    v = float(np.degrees(np.arctan2(y, x)))

    #  convert 0-360 range
    if v < 0.0:
        v += 360.0

    if not string:
        return v
    else:
        return '{:.1f}'.format(v)
