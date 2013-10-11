import AXUV
import surface
import geometry
import scipy

def BPLY(temp, place=(1.87,0,.157277), angle=(0,-.17453-scipy.pi/2,1.62385-scipy.pi/2)):


    pos = geometry.Origin(place,temp,angle=angle)
    area = [4e-3,3e-3]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    offset = geometry.Origin((-2.5e-3,0.,-7.389e-2),aperature,Vec=Vec)
    diodes = AXUV.AXUV22(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes
