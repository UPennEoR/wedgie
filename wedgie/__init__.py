"""init file for wedgie"""
import os
import sys
from .cosmo_utils import *
from .gen_utils import *
from .wedge_utils import *

ROOT = os.path.abspath(os.path.dirname(__file__))
CAL_DIR = os.path.join(ROOT, 'calibrations')
sys.path.append(CAL_DIR)

__version__ = 0.1
