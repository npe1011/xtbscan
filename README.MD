# GUI Tool for PES Scan with xTB or uma/omol25

- [GUI Tool for PES Scan with xTB or uma/omol25](#gui-tool-for-pes-scan-with-xtb-or-umaomol25)
  - [1. Introduction](#1-introduction)
    - [1.1. Overview](#11-overview)
    - [1.2. Update History](#12-update-history)
      - [2025/8/14](#2025814)
      - [2024/10/24](#20241024)
      - [2024/10/18](#20241018)
      - [2024/9/6](#202496)
      - [2022/11/11](#20221111)
  - [2. Requirements and Launch](#2-requirements-and-launch)
    - [2.1. Tested Environment](#21-tested-environment)
    - [2.2. Setup and configuration of xTB](#22-setup-and-configuration-of-xtb)
    - [2.3. Setup and configuration of uma/omol25](#23-setup-and-configuration-of-umaomol25)
  - [3. Usage](#3-usage)
    - [3.1. Launch](#31-launch)
    - [3.2. Load Initial Structure File: Input XYZ](#32-load-initial-structure-file-input-xyz)
    - [3.3.1 Setting Calculation Conditions (XTB)](#331-setting-calculation-conditions-xtb)
    - [3.3.2 Setting Calculation Conditions (UMA)](#332-setting-calculation-conditions-uma)
    - [3.4. Setting Scan Conditions](#34-setting-scan-conditions)
    - [3.4. Additional Parameter Fixing: Constrain Conditions](#34-additional-parameter-fixing-constrain-conditions)
    - [3.5. Run Calculations](#35-run-calculations)
    - [3.6. Visualize Results](#36-visualize-results)
  - [4. Other Notes](#4-other-notes)


## 1. Introduction

### 1.1. Overview

This is a Windows GUI tool for performing and analyzing so-called Scan
(relaxed scan) calculations using Grimme's xTB or uma/omol25 (ASE +
torch + fairchem-core). It supports 1D and 2D scans, as well as
multi-parameter concerted scans.

### 1.2. Update History

#### 2025/8/14

- Added: Calculation wiht uma/omol25 (via ASE).
- Added: fix atomic coordinates (constraint: atoms).
  This function does not always work well especially when combined with other constrains (not recommended).
- Fixed: Behavior when stopping calculations.
- Changed: the default CSV encoding (BOM-UTF-8 on Windows).

#### 2024/10/24

- Modified: 2D scan saddle point detection is re-performed with
  different thresholds when plotting. This is partially parallelized..

#### 2024/10/18

- Changed: 2D scan saddle point detection now uses SciPy's spline fit.
  More reliable results are obtained than before but it takes a longer time.

#### 2024/9/6

- Changed: GUI design and refactoring
- Added: reading input files other than Gaussian or xyz. 
- Changed: When reading xyz files as input, the last structure is used, instead of the first one.

#### 2022/11/11

- Initial version

## 2. Requirements and Launch

### 2.1. Tested Environment

- Windows 11
- Python 3.12.8
- xtb-6.6.1 for Windows
- wxPython==4.2.2
- numpy==2.2.1
- scipy==1.15.0
- cclib==1.8.1
- matplotlib==3.10.0

For uma/omol25,
- ase==3.26.0
- torch==2.6.0
- fairchem-core==2.3.0

```
python -m pip install wxpython cclib matplotlib scipy
```

may work. Errors may occur with older Python versions.

When using uma/omol25 (CPU mode), also install as following.

```
python -m pip install torch==2.6.0 fairchem-core==2.3.0
```

If you want GPU support, install the GPU version of torch==2.6.0 by yourself (not tested).


### 2.2. Setup and configuration of xTB

Download the Windows binary from the "release" at
https://github.com/grimme-lab/xtb and place it in a suitable directory.
Here we assume `D:\programs\xtb-6.6.1`. Also setup a suitable molecular viewer (jmol recommended).

Open `config.py` in the xtbscan folder with a text editor and update the parameters to match your environment.

```
VIEWER_PATH = 'D:/programs/jmol/jmol.bat'
XTB_BIN = 'D:/programs/xtb-6.6.1/bin/xtb.exe'
XTB_PARAM_DIR = 'D:/programs/xtb-6.6.1/share/xtb'
```

To disable SciPy-based saddle point detection in 2D scans, set the following to False.
In that case, SciPy is not required. The next parameter is the default criterion for saddle detection (smaller is
stricter):

```
USE_SCIPY = True
CHECK_SADDLE2D_GRAD_TOL = 0.001
```

### 2.3. Setup and configuration of uma/omol25

Download the parameter file from HuggingFace from [https://huggingface.co/facebook/UMA\] (registration required) and place it in a suitable location.
The small model is recommended. A high-performance GPU might allow the use of the medium model.

Specify the full path to the parameter file in `config.py`.
Set the following variable to True if you have a GPU environment and want to use it:

```
UMA_PARAM_PATH = 'D:/programs/uma/uma-s-1p1.pt'
UMA_USE_GPU = False
```

## 3. Usage

### 3.1. Launch

Just run `xtbscan.pyw`.

### 3.2. Load Initial Structure File: Input XYZ

- Load a file by dropping it or by using the top-left `...` button.
- An initial structure file can be xyz, Gaussian Job/Log files, or any
  file readable via cclib [https://cclib.github.io/].
- Files with extensions gjf/com/log/out are treated as Gaussian files;
  xyz files are treated as xyz; all others are read with cclib.
- If multiple structures are present in the file, the last structure is used.
- When reading Non-xyz input files, the corresponding xyz file is automatically generated and set as the input file.
- For multi-structure xyz files, a `name_last.xyz` file is generated and used.
- In these cases, a file with the same name will be overwritten.

### 3.3.1 Setting Calculation Conditions (XTB)

- `Method`: xtb parameters. `gfn2` is generally recommended. Some transition metal complexes favors `gfn1`. 
- `Charge`: Total charge.
- `UHF`: The number of unpaired electrons (0 for closed-shell, 1 for simple radical, etc.). This is different from spin multiplicity specification used in Gaussian or ORCA.
- `Solvation`: ALPB or GBSA; ALPB (newer) is recommended if implicit solvation is necessary.
- `force`: The constraint force (default 1.0). With xTB, geometrical parameters are not full fixed but constrained with a strong potential. Not necessary to change in most cases.

### 3.3.2 Setting Calculation Conditions (UMA)

- `Charge`: Total charge.
- `Multiplicity`: Spin multiplicity (same as Gaussian or ORCA).
- `Fmax`: Optimization convergence threshold (smaller is stricter)
- `Max Cycles`: The maximum steps of optimization.

### 3.4. Setting Scan Conditions

- You can scan distances (2 atoms), angles (3 atoms), and dihedrals (4 atoms).
- Enter parameters and click `Add` to add them to the table.
- Atoms are indexed from 1, separated by spaces or commas.
- `Start`/`End`/`Num Step` specify the initial value, final value, and the number of scan steps.
- With a structure loaded, you can click `view` button to open the viewer, which helps to input atom indices.
- Click `current` with a loaded structure sets `Start` with the current value.
- One condition: 1D scan, Two conditions: 2D scan. 3D or higher dimension scans are not possible.
- Checking `concerted mode` enables concerted scans, where multiple geometrical parameters change synchronously
  in a 1D scan. `Num Step` must be identical for all scan paramters.
- If no scan conditions are given, a simple geometry optimization is performed.
- With xTB, scans are run via the built-in function of the xTB binary.
  With uma/omol25, scans are done by repeated ASE constraint optimizations, so behaviors may differ.

### 3.4. Additional Parameter Fixing: Constrain Conditions

- Additional geometrical parameters can be fixed.
- For fixing at the current value, you can use current or auto (auto is not available in scan settings).
- Complex constraint/scan combinations are not well tested. Atoms constraints may not work when combined with
  other scan or constraints (especially when u UMA).

### 3.5. Run Calculations

- Set `CPUs` and `Memory per CPU`, then click `Run`.
- `CPUs` should equal or be less than physical CPU cores (e.g., max 4 for 4-core/8-thread CPUs).
- Total Memory usage = input memory × CPU cores.
- `Keep Log` specifies whether to keep intermediate xTB or ASE optimization logs. `When fail` keeps them only on error termination.
- Click `Run` to start. Click `Stop` to cancel (cannot resume from partial progress).

### 3.6. Visualize Results

- When performing simple optimization, only a optimized structure xyz file is given. 
- When successful termination of scans, you are asked whether to load results.
- When the results loaded, the CSV file name appears in the top-right CSV text field.
- Results are saved as `JOBNAME.csv` (energies) and `JOBNAME.xyz` (structures).
- Previous results can be loaded by opening CSV files (drag-and-drop).
- Buttons under the CSV field allow visualization:
  - `view (all)`: open all structures with the viewer.
  - `show table`: display the CSV contents in a table. You can also see the file with Excel or other softwares.
  - Specify a sacn step number and the click `view` or `copy` or `save` to see the structure, copy the coordinates to clipboard, or save as a file.
  - `plot`: Plot scaned parameter/energy graphs (1D scan: normal graph, 2D: 3D graph). The energies are given as the relative values in kcal/mol.
  - `annotation`: Add structure number labels when plotting.
    Saddle-like points are shown as red points.
    In 1D scans, red points are just local maxima.
    In 2D scans, spline fitting and gradient/Hessian-based saddle detection is performed (if USE_SCIPY=True).
    In 2D scans, you can re-run saddle detection with a custom `grad. tol.` when plotting.
  - `surface`: shows PES-like surface for 2D scans.
  - `grad. tol.` is a scaled value without physical units.

## 4. Other Notes

- A dropped xyz file is loaded as the input file. 
- When a CSV file loaded, the corresponding (sane name) xyz file is loaded as the result file.
- In the result CSV, `scan` is the parameter value from settings, `real` is from the actual optimized geometry. 
  Because parameters are just constrained via strong potentials when using xtb, there may be some deviations. `real` values are used when plotting.
