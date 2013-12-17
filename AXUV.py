# setup file necessary to generate the AXUV20 geometry for BP-LY and will
# act as a rosetta stone for generating other detector geometries.
import geometry
import surface
import scipy
# aperature position in RZ coordinates, will be redefined later with CAD work
#


def AXUV20(temp):
    diodes = 20*[0]
    area = [4e-3,7.5e-4]
    Vec = [geometry.Vecx((0.,1.,0.)),geometry.Vecx((0.,0.,1.))]
    spacing = 9.5e-4
    pos = scipy.mgrid[-9.5:9.5:1]*spacing
    #pos = scipy.linspace(-9.025e-3,9.025e-3,20)
    
    for i in range(len(diodes)):
        diodes[i] = surface.Rect((pos[i],0.,0.),temp,area,Vec=Vec)
    return diodes

def AXUV22(temp):
    """ temp """
    diodes = 22*[0]
    area = [4.4e-3,1e-3]
    spacing = 2e-3
    Vec = [geometry.Vecx((0.,1.,0.)),geometry.Vecx((0.,0.,1.))]
    pos = scipy.mgrid[-10.5:10.5:1]*spacing 
    #pos = scipy.linspace(-2.1e-2,2.1e-2,22)
    for i in range(len(diodes)):
        diodes[i] = surface.Rect((pos[i],0.,0.),temp,area,Vec=Vec)
    return diodes
