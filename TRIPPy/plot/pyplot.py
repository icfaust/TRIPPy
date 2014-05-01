import scipy
import TRIPPy
import matplotlib.pyplot as plt

def plotTokamak(tokamak, axis=True, *args, **kwargs):

    plt.plot(tokamak.sagi.s, tokamak.norm.s, *args, **kwargs)
    if axis:
        plt.axis('equal')

def plotLine(line, invessel=True, ds=2.5e-3, pargs=None, **kwargs):
    
    try:
        if invessel:
            temp = line(scipy.mgrid[line.norm.s[-2]:line.norm.s[-1]:ds])
        else:
            temp = line(scipy.mgrid[line.norm.s[0]:line.norm.s[-1]:ds])
        
        if not pargs is None:
            plt.plot(temp.r0(), temp.x2(), pargs, **kwargs)
        else:
            plt.plot(temp.r0(), temp.x2(), **kwargs)

    except AttributeError:
        for i in line:
            plotLine(i, invessel=invessel, pargs=pargs, **kwargs)

def sinogram(beam, r, z, invessel=True, ds=2.5e-3, pargs=None, **kwargs):
    
    try:
        if invessel:
            temp = beam(scipy.mgrid[beam.norm.s[-2]:beam.norm.s[-1]:ds])
        else:
            temp = beam(scipy.mgrid[beam.norm.s[0]:beam.norm.s[-1]:ds])

        # deal with branch cut
        temp0 = temp.t0(r, z)
        temp2 = temp.t2(r, z)
        temp = scipy.arange(temp0.size)[abs(temp2[1:] - temp2[:-1]) > scipy.pi]
        if len(temp) > 0:
            temp0 = scipy.insert(temp0, temp, None)
            temp2 = scipy.insert(temp2, temp, None)

        if not pargs is None:
            plt.plot(temp2,temp0, pargs, **kwargs)
        else:
            plt.plot(temp2,temp0, **kwargs)

    except AttributeError:
        for i in beam:
            sinogram(i, r, z, invessel=invessel, pargs=pargs, **kwargs)
