FLOAT = 'float64'
JOULE_TO_KCAL = 2.390E-4
HARTREE_TO_JOULE_PER_MOL = 2.6255E6
HARTREE_TO_KCAL = 627.51

XTB_NAME = 'xtb'
XTB_OPT_FILE = 'xtbopt.xyz'
XTB_SCAN_FILE = 'xtbscan.log'
XTB_LOG_FILE = 'xtblog.txt'

XYZ_FORMAT = '{:<2s}  {:>12.8f}  {:>12.8f}  {:>12.8f}\n'

INIT_XYZ_FILE = 'init.xyz'
INPUT_FILE = 'input.txt'

OPT_XYZ_SUFFIX = '_xtbopt.xyz'
RESULT_XYZ_SUFFIX = '_xtbscan.xyz'
RESULT_CSV_SUFFIX = '_xtbscan.csv'
STOP_FILE_SUFFIX = '_stopmessage.dat'
WORKDIR_PREFIX = '__'
WORKDIR_SUFFIX = '_XTB__'
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
