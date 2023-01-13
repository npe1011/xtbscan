FLOAT = 'float64'
JOULE_TO_KCAL = 2.390E-4
HARTREE_TO_JOULE_PER_MOL = 2.6255E6
HARTREE_TO_KCAL = 627.51
ANG_TO_BOHR = 1.8897259886

XTB_NAME = 'xtb'
XTB_OPT_FILE = 'xtbopt.xyz'
XTB_SCAN_FILE = 'xtbscan.log'
XTB_LOG_FILE = 'xtblog.txt'
XTB_ENERGY_FILE = 'energy'
XTB_GRADIENT_FILE = 'gradient'
XTB_HESSIAN_FILE = 'hessian'

XYZ_FORMAT = '{:<2s}  {:>12.8f}  {:>12.8f}  {:>12.8f}\n'

INIT_XYZ_FILE = 'init.xyz'
INPUT_FILE = 'input.txt'

OPT_XYZ_SUFFIX = '_xtbopt.xyz'

OPTTS_XYZ_SUFFIX = '_xtboptts.xyz'
OPTTS_TRAJXYZ_SUFFIX = '_xtboptts_traj.xyz'
OPTTS_FINAL_HESSIAN_SUFFIX = '_xtboptts_hess.txt'
OPTTS_GEOMETRIC_LOG_SUFFIX = '_xtboptts.log'
OPTTS_RESULT_CSV_SUFFIX = '_xtboptts.csv'

RESULT_XYZ_SUFFIX = '_xtbscan.xyz'
RESULT_CSV_SUFFIX = '_xtbscan.csv'
STOP_FILE_SUFFIX = '_stopmessage.dat'
WORKDIR_PREFIX = '__'
WORKDIR_SUFFIX = '_XTB__'
GEOMETRIC_DIR_SUFFIX = '_geometric'
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
