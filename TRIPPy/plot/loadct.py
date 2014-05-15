import scipy
from scipy.io import readsav
import matplotlib.colors
import os
import sys
import inspect
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))

def loadct(num, **kwargs):
    output = readsav(cmd_folder+'/idl_colors.sav').output.T
    return matplotlib.colors.LinearSegmentedColormap.from_list('idl', output[num].T, **kwargs)
