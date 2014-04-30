import scipy
import TRIPPy
import matplotlib.pyplot as plt


def plotTokamak(tokamak, axis=True, *args, **kwargs):

    plt.plot(tokamak.sagi.s, tokamak.norm.s, *args, **kwargs)
    if axis:
        plt.axis('equal')

def plotLine(line, invessel=True, ds=2.5e-3, *args, **kwargs):
    
    try:
        if invessel:
            temp = line(scipy.mgrid[line.norm.s[-2]:line.norm.s[-1]:ds])
        else:
            temp = line(scipy.mgrid[line.norm.s[0]:line.norm.s[-1]:ds])
        plt.plot(temp.r0(), temp.x2(), *args, **kwargs)
    except AttributeError:
        for i in line:
            plotLine(i, invessel=invessel, *args, **kwargs)

def sinogram(beam, r, z, invessel=True, ds=2.5e-3, *args, **kwargs):
    
    try:
        if invessel:
            temp = beam(scipy.mgrid[beam.norm.s[-2]:beam.norm.s[-1]:ds])
        else:
            temp = beam(scipy.mgrid[beam.norm.s[0]:beam.norm.s[-1]:ds])
        plt.plot(temp.t2(r, z), temp.t0(r, z), *args, **kwargs)
    except AttributeError:
        for i in beam:
            sinogram(i, r, z, invessel=invessel, *args, **kwargs)
