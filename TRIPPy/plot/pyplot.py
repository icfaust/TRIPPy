import scipy
import TRIPPy
import matplotlib.pyplot as plt


def plotTokamak(tokamak, axis=True, **kwargs):

    plt.plot(tokamak.sagi.s,tokamak.norm.s,**kwargs)
    if axis:
        plt.axis('equal')

def plotLine(line, invessel=True, ds=2.5e-3, **kwargs):
    
    try:
        if invessel:
            temp = line(scipy.mgrid[line.norm.s[-2]:line.norm.s[-1]:ds]).r()
        else:
            temp = line(scipy.mgrid[line.norm.s[0]:line.norm.s[-1]:ds]).r()
        plt.plot(temp[0],temp[2],**kwargs)
    except AttributeError:
        for i in line:
            plotLine(i, invessel=invessel, **kwargs)

def sinogram(beam, r, z, invessel=True, **kwargs):
    
    try:
        if invessel:
            temp = line(scipy.mgrid[line.norm.s[-2]:line.norm.s[-1]:ds]).t(r, z)
        else:
            temp = line(scipy.mgrid[line.norm.s[0]:line.norm.s[-1]:ds]).t(r, z)
        plt.plot(temp[2],temp[0],**kwargs)
    except AttributeError:
        for i in line:
            sinogram(i, r, z, invessel=invessel, **kwargs)
