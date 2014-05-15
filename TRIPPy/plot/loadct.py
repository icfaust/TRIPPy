import scipy
import matplotlib.colors
import os
import sys
import inspect
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))

def loadct(num, **kwargs):
    output = scipy.genfromtxt(cmd_folder+'/idl_colors.txt',
                              skip_header=256*num,
                              skip_footer=(38-num)*256)/255.
    return matplotlib.colors.LinearSegmentedColormap.from_list('idl', output, **kwargs)
