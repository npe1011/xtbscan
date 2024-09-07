# Required PATHs settings
VIEWER_PATH = 'D:/programs/jmol/jmol.bat'
XTB_BIN = 'D:/programs/xtb-6.6.1/bin/xtb.exe'
XTB_PARAM_DIR = 'D:/programs/xtb-6.6.1/share/xtb'
XTB_OTHER_LIB_DIR = None  # in case other library directory (MinGW etc.) has  to be in path to run xtbscan.

# Customization of plotting
SCAN_PLOT_1D_WIDTH = 6
SCAN_PLOT_1D_HEIGHT = 4
SCAN_PLOT_1D_COLOR = 'skyblue'
SCAN_PLOT_1D_TS_COLOR = 'firebrick'

SCAN_PLOT_2D_WIDTH = 6
SCAN_PLOT_2D_HEIGHT = 4
SCAN_PLOT_2D_COLOR = 'skyblue'
SCAN_PLOT_2D_TS_COLOR = 'firebrick'

# Constants (not necessary to edit)
FLOAT = float
JOULE_TO_KCAL = 2.390E-4
HARTREE_TO_JOULE_PER_MOL = 2.6255E6
HARTREE_TO_KCAL = 627.51
ANG_TO_BOHR = 1.8897259886

XTB_NAME = 'xtbscan'
XTB_OPT_FILE = 'xtbopt.xyz'
XTB_SCAN_FILE = 'xtbscan.log'
XTB_LOG_FILE = 'xtblog.txt'
XTB_INPUT_FILE = 'xtbinp.txt'

XYZ_FORMAT = '{:<2s}  {:>12.8f}  {:>12.8f}  {:>12.8f}\n'

INIT_XYZ_FILE = 'init.xyz'
INPUT_FILE = 'input.txt'
STOP_FILE_SUFFIX = '_stopmessage.dat'
STOP_CHECK_INTERVAL = 1  # sec

XTB_SOLVENT_LIST = [
    'Acetone',
    'Acetonitrile',
    'Aniline',
    'Benzaldehyde',
    'Benzene',
    'CH2Cl2',
    'CHCl3',
    'CS2',
    'Dioxane',
    'DMF',
    'DMSO',
    'Ether',
    'Ethylacetate',
    'Furane',
    'Hexadecane',
    'Hexane',
    'Methanol',
    'Nitromethane',
    'Octanol',
    'Phenol',
    'Toluene',
    'THF',
    'Water'
]

ATOM_LIST = ['bq', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
             'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
             'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
             'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
             'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
             'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs',
             'Mt', 'Ds', 'Rg', 'Cn']
