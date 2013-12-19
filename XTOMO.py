import geometry
import surface


def XTOMOchip(temp):
    diodes = 38*[0]
    area = [4e-3,.94e-3]
    Vec = [geometry.Vecx((0.,1.,0.)),geometry.Vecx((0.,0.,1.))]
    spacing = 1e-3 
    pos = scipy.mgrid[-16.5:16.5:38j]*spacing
    #pos = scipy.linspace(e-3,e-3,38)
    for i in range(len(diodes)):
        vecin = geometry.Vecx((0.,pos[i],0.))
        diodes[i] = surface.Rect(vecin,temp,area,Vec=Vec)
    return diodes

def XTOMO1(temp, place=(,0,),angle =(,,)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place,temp,angle=angle)
    area = [  ]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    offset = goemetry.Origin((,,),aperature,Vec=Vec)
    diodes = XTOMOchip(offset)
    for i in diodes:
        i.redefine(temp)
        diodes.append(aperature)
        diodes[-1].redefine(temp)
    return diodes

def XTOMO2(temp, place=(,0,),angle =(,,)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place,temp,angle=angle)
    area = [  ]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    offset = goemetry.Origin((,,),aperature,Vec=Vec)
    diodes = XTOMOchip(offset)
    for i in diodes:
        i.redefine(temp)
        diodes.append(aperature)
        diodes[-1].redefine(temp)
    return diodes

def XTOMO3(temp, place=(,0,),angle =(,,)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place,temp,angle=angle)
    area = [  ]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    offset = goemetry.Origin((,,),aperature,Vec=Vec)
    diodes = XTOMOchip(offset)
    for i in diodes:
        i.redefine(temp)
        diodes.append(aperature)
        diodes[-1].redefine(temp)
    return diodes

def XTOMO4(temp, place=(,0,),angle =(,,)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place,temp,angle=angle)
    area = [  ]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    offset = goemetry.Origin((,,),aperature,Vec=Vec)
    diodes = XTOMOchip(offset)
    for i in diodes:
        i.redefine(temp)
        diodes.append(aperature)
        diodes[-1].redefine(temp)
    return diodes

def XTOMO5(temp, place=(,0,),angle =(,,)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place,temp,angle=angle)
    area = [  ]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    offset = goemetry.Origin((,,),aperature,Vec=Vec)
    diodes = XTOMOchip(offset)
    for i in diodes:
        i.redefine(temp)
        diodes.append(aperature)
        diodes[-1].redefine(temp)
    return diodes
