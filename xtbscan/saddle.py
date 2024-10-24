import os
from typing import Optional
from multiprocessing import Pool, cpu_count

import numpy as np

import config

if config.USE_SCIPY:
    try:
        from scipy.interpolate import RectBivariateSpline
    except:
        IMPORT_SCIPY = False
    else:
        IMPORT_SCIPY = True


def _set_threads(num_procs: Optional[int]):
    if num_procs is None:
        num_procs = cpu_count()
    num_procs = str(num_procs)
    os.environ['OMP_NUM_THREADS'] = num_procs
    os.environ['MKL_NUM_THREADS'] = num_procs


def check_saddle_1d(energies: np.ndarray) -> np.ndarray:
    saddle_check_list = np.full_like(energies, 0, dtype=int)
    for i in range(1, len(energies) - 1):
        if energies[i - 1] < energies[i] and energies[i + 1] < energies[i]:
            saddle_check_list[i] = 1
    return saddle_check_list


def check_saddle_2d(energies: np.ndarray, num_dim1: int, num_dim2: int,
                    grad_tol: Optional[float], num_procs: Optional[int]) -> np.ndarray:

    if config.USE_SCIPY and IMPORT_SCIPY:
        result = _check_saddle_2d_with_spline_fit(energies, num_dim1, num_dim2, grad_tol, num_procs)
        return result
    elif config.USE_SCIPY:
        raise ImportError('config.USE_SCIPY is True but import failed. Install scipy or set USE_SCIPY=False.')
    else:
        return _check_saddle_2d_without_spline_fit(energies, num_dim1, num_dim2)


def _check_saddle_2d_with_spline_fit(energies: np.ndarray, num_dim1: int, num_dim2: int,
                                     grad_tol: Optional[float], num_procs: Optional[int]) \
        -> np.ndarray:

    # energies are standardized and reshaped
    _energies = energies - np.mean(energies)
    scale = np.max(np.abs(_energies))
    if not np.isclose(scale, 0.0):
        _energies /= scale
    _energies = _energies.reshape(num_dim1, num_dim2)

    saddle_check_list = np.full_like(_energies, 0, dtype=int)

    # Some constants to adjust
    param_max = 10.0
    num_mesh = 1000

    if grad_tol is None:
        grad_tol = config.CHECK_SADDLE2D_GRAD_TOL

    xs = np.linspace(0.0, param_max, num_dim1)
    ys = np.linspace(0.0, param_max, num_dim2)
    spline = RectBivariateSpline(xs, ys, _energies)

    mesh_xs = np.linspace(0.0, param_max, num_mesh)
    mesh_ys = np.linspace(0.0, param_max, num_mesh)

    # serialization
    data_list = [(spline, x, y, xs, ys, grad_tol) for x in mesh_xs for y in mesh_ys]
    # calculate
    _set_threads(1)  # multiprocessing with 1 threads
    with Pool(processes=num_procs) as p:
        result = p.starmap(_check_saddle_for_pal, data_list)
    # read result and set saddle_check_list
    for res in result:
        if res is not None:
            saddle_check_list[res[0], res[1]] = 1
    return saddle_check_list.flatten()


def _check_saddle_for_pal(spline, x, y, xs, ys, grad_tol) -> Optional[tuple]:
    # first derivatives and stationary point check
    dx = float(spline(x, y, dx=1, dy=0))
    dy = float(spline(x, y, dx=0, dy=1))
    df_norm = np.sqrt(dx ** 2 + dy ** 2)
    if df_norm > grad_tol:
        return None

    # Hessian and saddle check
    dx2 = float(spline(x, y, dx=2, dy=0))
    dy2 = float(spline(x, y, dx=0, dy=2))
    dxdy = float(spline(x, y, dx=1, dy=1))
    hessian = np.array([[dx2, dxdy], [dxdy, dy2]])
    eigenvalues = np.linalg.eigvals(hessian)
    if np.all(eigenvalues > 0) or np.all(eigenvalues < 0):
        return None

    # Nearest point
    ix = np.abs(xs - x).argmin()
    iy = np.abs(ys - y).argmin()

    return ix, iy


def _check_saddle_2d_without_spline_fit(energies: np.ndarray, num_dim1: int, num_dim2: int) -> np.ndarray:

    _energies = energies.reshape((num_dim1, num_dim2))
    saddle_check_list = np.full_like(_energies, 0, dtype=int)

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
                saddle_check_list[x, y] = 1

    return saddle_check_list.flatten()
