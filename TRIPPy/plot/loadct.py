import scipy
import matplotlib.colors
import os
import sys
import inspect
import matplotlib.pyplot as plt

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))

def loadct(num, r=False, mayavi=False, **kwargs):
    if not mayavi:

        output = scipy.genfromtxt(cmd_folder+'/idl_colors.txt',
                                  skip_header=256*num,
                                  skip_footer=(39-num)*256)/255.
        if r:
            output = output[::-1]
            
        return matplotlib.colors.LinearSegmentedColormap.from_list('idl', output, **kwargs)
    else:
        output = scipy.ones((256,4),dtype=int)
        output[:,0:3] = scipy.genfromtxt(cmd_folder+'/idl_colors.txt',
                                         skip_header=256*num,
                                         skip_footer=(39-num)*256,dtype=int)
        if r:
            output = output[::-1]

        return output

def loadgs(num, r=False, mayavi=False, **kwargs):
    if not mayavi:

        output = scipy.genfromtxt(cmd_folder+'/gs_colors.txt',
                                  skip_header=256*num,
                                  skip_footer=(15-num)*256)
        if r:
            output = output[::-1]
            
        return matplotlib.colors.LinearSegmentedColormap.from_list('gs', output,**kwargs)
    else:
        output = scipy.ones((256,4),dtype=int)
        output[:,0:3] = scipy.genfromtxt(cmd_folder+'/gs_colors.txt',
                                         skip_header=256*num,
                                         skip_footer=(15-num)*256,dtype=int)
        if r:
            output = output[::-1]

        return output

def showct():
    print('not yet')
