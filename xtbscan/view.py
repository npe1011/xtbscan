import csv
from pathlib import Path
import subprocess as sp
import tempfile
from typing import Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D

import config
from xtbscan import xyzutils, utils


def plot_scan(csv_file: Union[str, Path], annotation: bool = True) -> None:
    csv_file = Path(csv_file)
    with csv_file.open(mode='r') as f:
        mode = f.readline().lower().split(',')[0].strip()
        if mode == '1d':
            _plot_scan_1d(csv_file, annotation=annotation)
        elif mode == '2d':
            _plot_scan_2d(csv_file, annotation=annotation)
        else:
            _plot_scan_concerted(csv_file, annotation=annotation)


def _plot_scan_1d(csv_file: Path, annotation: bool) -> None:
    energies = []
    parameters = []
    saddle_check_list = []
    with csv_file.open(mode='r') as f:
        reader = csv.reader(f)
        for i, line in enumerate(reader):
            if i == 0:
                continue
            if i == 1:
                parameter_name = line[2].split(maxsplit=1)[1].strip()  # line[2] ~ "[scan] distance N10-C25 (ang.)"
            else:
                parameters.append(float(line[2]))
                energies.append(float(line[-2]))
                saddle_check_list.append(int(line[-1]))

    names = [str(n + 1) for n in range(len(energies))]

    fig = plt.figure(figsize=(config.SCAN_PLOT_1D_WIDTH, config.SCAN_PLOT_1D_HEIGHT))
    try:
        fig.canvas.manager.set_window_title(csv_file.name)
    except:
        pass
    ax = fig.add_subplot(111)

    if annotation:
        for (x, y, name, saddle) in zip(parameters, energies, names, saddle_check_list):
            if saddle == 1:
                ax.plot(x, y, 'o', color=config.SCAN_PLOT_1D_TS_COLOR)
                ax.annotate(name, xy=(x, y), color=config.SCAN_PLOT_1D_TS_COLOR)
            else:
                ax.plot(x, y, 'o', color=config.SCAN_PLOT_1D_COLOR)
                ax.annotate(name, xy=(x, y))
    else:
        ax.scatter(parameters, energies, marker='o', color=config.SCAN_PLOT_1D_COLOR)

    ax.set_xlabel(parameter_name)
    ax.set_ylabel('rel. energy (kcal/mol)')

    plt.tight_layout()
    plt.show()


def _plot_scan_2d(csv_file: Path, annotation: bool) -> None:
    energies = []
    parameters1 = []
    parameters2 = []
    saddle_check_list = []
    with csv_file.open(mode='r') as f:
        reader = csv.reader(f)
        for i, line in enumerate(reader):
            if i == 0:
                continue
            if i == 1:
                parameter1_name = line[2].split(maxsplit=1)[1].strip()  # line[2] ~ "[scan] distance N10-C25 (ang.)"
                parameter2_name = line[4].split(maxsplit=1)[1].strip()  # line[4] ~ "[scan] distance N10-C25 (ang.)"
            else:
                parameters1.append(float(line[2]))
                parameters2.append(float(line[4]))
                energies.append(float(line[-2]))
                saddle_check_list.append(int(line[-1]))

    names = [str(n + 1) for n in range(len(energies))]

    fig = plt.figure(figsize=(config.SCAN_PLOT_2D_WIDTH, config.SCAN_PLOT_2D_HEIGHT))
    try:
        fig.canvas.manager.set_window_title(csv_file.name)
    except:
        pass
    ax = fig.add_subplot(111, projection='3d')

    for (x, y, z, name, saddle) in zip(parameters1, parameters2, energies, names, saddle_check_list):
        if annotation:
            if saddle == 1:
                ax.scatter3D(x, y, z, marker='o', color=config.SCAN_PLOT_2D_TS_COLOR)
                ax.text(x, y, z, name, color=config.SCAN_PLOT_2D_TS_COLOR)
            else:
                ax.scatter3D(x, y, z, marker='o', color=config.SCAN_PLOT_2D_COLOR)
                ax.text(x, y, z, name)
        else:
            ax.scatter3D(x, y, z, marker='o', color=config.SCAN_PLOT_2D_COLOR)

    max_z = max(energies)
    min_z = min(energies)
    buff = (max_z - min_z)*0.05
    ax.set_zlim(min_z-buff, max_z+buff)

    ax.set_xlabel(parameter1_name)
    ax.set_ylabel(parameter2_name)
    ax.set_zlabel('rel. energy (kcal/mol)')

    plt.tight_layout()
    plt.show()


def _plot_scan_concerted(csv_file: Path, annotation: bool) -> None:
    energies = []
    saddle_check_list = []
    with csv_file.open(mode='r') as f:
        reader = csv.reader(f)
        for i, line in enumerate(reader):
            if i <= 1:
                continue
            else:
                energies.append(float(line[-2]))
                saddle_check_list.append(int(line[-1]))

    parameters = [n + 1 for n in range(len(energies))]
    names = [str(n) for n in parameters]

    fig = plt.figure(figsize=(config.SCAN_PLOT_1D_WIDTH, config.SCAN_PLOT_1D_HEIGHT))
    try:
        fig.canvas.manager.set_window_title(csv_file.name)
    except:
        pass
    ax = fig.add_subplot(111)

    if annotation:
        for (x, y, name, saddle) in zip(parameters, energies, names, saddle_check_list):
            if saddle == 1:
                ax.plot(x, y, 'o', color=config.SCAN_PLOT_1D_TS_COLOR)
                ax.annotate(name, xy=(x, y), color=config.SCAN_PLOT_1D_TS_COLOR)
            else:
                ax.plot(x, y, 'o', color=config.SCAN_PLOT_1D_COLOR)
                ax.annotate(name, xy=(x, y))
    else:
        ax.scatter(parameters, energies, marker='o', color=config.SCAN_PLOT_1D_COLOR)

    ax.set_xlabel('# Steps')
    ax.set_ylabel('rel. energy (kcal/mol)')

    plt.tight_layout()
    plt.show()


def plot_surface(csv_file: Union[str, Path]) -> None:
    csv_file = Path(csv_file)
    energies = []
    parameters1 = []
    parameters2 = []
    with csv_file.open(mode='r') as f:
        reader = csv.reader(f)
        for i, line in enumerate(reader):
            if i == 0:
                if line[0].lower() != '2d':
                    raise ValueError('Only 2D scan file can be visualized by surface plot.')
                num_dim1 = int(line[1])
                num_dim2 = int(line[2])
                continue
            elif i == 1:
                parameter1_name = line[2].split(maxsplit=1)[1].strip()  # line[2] ~ "[scan] distance N10-C25 (ang.)"
                parameter2_name = line[4].split(maxsplit=1)[1].strip()  # line[4] ~ "[scan] distance N10-C25 (ang.)"
            else:
                parameters1.append(float(line[2]))
                parameters2.append(float(line[4]))
                energies.append(float(line[-2]))

    fig = plt.figure(figsize=(config.SCAN_PLOT_2D_WIDTH, config.SCAN_PLOT_2D_HEIGHT))
    try:
        fig.canvas.manager.set_window_title(csv_file.name)
    except:
        pass
    ax = fig.add_subplot(111, projection='3d')

    x1 = np.array(parameters1, dtype=config.FLOAT).reshape((num_dim1, num_dim2))
    x2 = np.array(parameters2, dtype=config.FLOAT).reshape((num_dim1, num_dim2))
    z = np.array(energies, dtype=config.FLOAT).reshape((num_dim1, num_dim2))

    ax.plot_surface(x1, x2, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    max_z = max(energies)
    min_z = min(energies)
    buff = (max_z - min_z) * 0.05
    ax.set_zlim(min_z - buff, max_z + buff)

    ax.set_xlabel(parameter1_name)
    ax.set_ylabel(parameter2_name)
    ax.set_zlabel('rel. energy (kcal/mol)')

    plt.tight_layout()
    plt.show()


def view_xyz_file(xyz_file: Union[str, Path]) -> None:
    sp.Popen([config.VIEWER_PATH, str(xyz_file)])


@utils.async_func
def view_xyz_structure(atoms: np.ndarray, coordinates: np.ndarray, title: str = ''):
    _, file = tempfile.mkstemp(suffix='.xyz', text=True)
    xyzutils.save_xyz_file(file, atoms, coordinates, title)
    sp.run([config.VIEWER_PATH, str(file)])
    try:
        Path(file).unlink()
    except:
        pass
