import geometry,surface,AXUV20
import scipy

def diode1(temp,place=(2.0,0,2e-2), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815,0,.5))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-2.25e-2),aperature,Vec=Vec)
    diodes = AXUV20.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes

def diode2(temp,place=(2.0,0,-2e-2), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815,0,-.5))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-2.25e-2),aperature,Vec=Vec)
    diodes = AXUV20.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes

def diode3(temp,place=(.8,0,1.8), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815-.075,0,0))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-1.25e-2),aperature,Vec=Vec)
    diodes = AXUV20.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes

def diode4(temp,place=(1.5,0,-1.5), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815-.4,0,-1.175))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-2.5e-2),aperature,Vec=Vec)
    diodes = AXUV20.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes
