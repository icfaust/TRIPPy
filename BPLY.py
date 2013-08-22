import AXUV20
import surface
import geometry
import scipy

def BPLY(temp):

    area = [5e-3,2e-3]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),temp,area,Vec=Vec)
    offset = geometry.Origin((-2.5e-3,0.,-7.389e-2),aperature,Vec=Vec)
    diodes = AXUV20.AXUV22(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    return diodes
