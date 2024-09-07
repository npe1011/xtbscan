import os
import subprocess as sp
from threading import Thread
from typing import List


def read_atom_list_string(atom_list_string: str) -> List[int]:
    """
    return atom number list (int list) from string as "1,3,5-10" "5 15 16" etc.
    as 0-based number list
    :param atom_list_string:
    :return: int list
    """
    blocks = [x.strip() for x in atom_list_string.strip().replace(',', ' ').split()]
    atoms = []
    for block in blocks:
        if '-' in block:
            start, end = [int(x) for x in block.split('-')]
            for i in range(start, end + 1):
                atoms.append(i-1)
        else:
            atoms.append(int(block)-1)
    return atoms


def expand_atom_list_string(atom_list_string: str) -> str:
    return ','.join([str(x+1) for x in read_atom_list_string(atom_list_string)])


def atom_indices_to_string(atom_indices: List[int]) -> str:
    terms = []
    current_start_index = None
    current_end_index = None
    for i in sorted(atom_indices):
        if current_start_index is None:
            current_start_index = i
            current_end_index = i
            continue
        else:
            if current_end_index == i-1:
                current_end_index = i
            else:
                if current_start_index == current_end_index:
                    terms.append(str(current_start_index+1))
                else:
                    terms.append(str(current_start_index+1) + '-' + str(current_end_index+1))
                current_start_index = i
                current_end_index = i
    if current_start_index is not None:
        if current_start_index == current_end_index:
            terms.append(str(current_start_index+1))
        else:
            terms.append(str(current_start_index+1) + '-' + str(current_end_index+1))

    return ','.join(terms)


def popen_bg(*args, **kwargs):
    if os.name == 'nt':
        startupinfo = sp.STARTUPINFO()
        startupinfo.dwFlags |= sp.STARTF_USESHOWWINDOW
        win_kwargs = {'startupinfo': startupinfo}
        return sp.Popen(*args, **kwargs, **win_kwargs)
    else:
        return sp.Popen(*args, **kwargs)


def async_func(func):
    def wrapper(*args, **kwargs):
        _func = Thread(target=func, args=args, kwargs=kwargs)
        _func.start()
        return _func

    return wrapper
