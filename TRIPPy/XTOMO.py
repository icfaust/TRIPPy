import geometry
import surface
import scipy
import beam

def XTOMOchip(temp):
    diodes = 38*[0]
    area = [4e-3,.94e-3]
    Vec = [geometry.Vecx((0.,1.,0.)),geometry.Vecx((0.,0.,1.))]
    spacing = 1e-3 
    pos = scipy.mgrid[-18.5:18.5:38j]*spacing
    #pos = scipy.linspace(e-3,e-3,38)
    for i in range(len(diodes)):
        vecin = geometry.Vecx((0.,pos[i],0.))
        diodes[i] = surface.Rect(vecin,
                                 temp,
                                 area,
                                 vec=Vec)
    return diodes

def XTOMO1(temp, place1=(0.73674,0.,0.46479), place2=(.747,0.,.4947), angle=(0.,-(180.5*scipy.pi/180.),scipy.pi/2.)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place1,temp,angle=angle)
    area = [1e-3,3e-3]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperture = surface.Rect((0.,0.,0.),pos,area,vec=Vec)
    offset = geometry.Origin(place2,temp,angle=angle)
    diodes = XTOMOchip(offset)
    diodes[0].redefine(temp)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperture)
    diodes[-1].redefine(temp)
    return diodes

def XTOMO1beam(plasma):
    temp = XTOMO1(plasma)
    output = beam.multiBeam(temp[:-1],temp[-1])
    plasma.trace(output)
    return output


#def XTOMO1(temp, place=(.73674,0,.46479),angle =(,,)):
#    # this is the hard part, tree storage values are not done in such a manner
#    pos = geometry.Origin(place,temp,angle=angle)
#    area = [3e-3,1e-3]
#    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
#    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
#    offset = goemetry.Origin((0.,0.,0.),aperature,Vec=Vec)
#    diodes = XTOMOchip(offset)
#    for i in diodes:
#        i.redefine(temp)
#        diodes.append(aperature)
#        diodes[-1].redefine(temp)
#    return diodes
#
#def XTOMO2(temp, place=(,0,),angle =(,,)):
#    # this is the hard part, tree storage values are not done in such a manner
#    pos = geometry.Origin(place,temp,angle=angle)
#    area = [3.6e-3,5e-4]
#    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
#    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
#    offset = goemetry.Origin((,,),aperature,Vec=Vec)
#    diodes = XTOMOchip(offset)
#    for i in diodes:
#        i.redefine(temp)
#        diodes.append(aperature)
#        diodes[-1].redefine(temp)
#    return diodes
#

def XTOMO3(temp, place=(1.0152,0.,0.),angle =(0.,scipy.pi/2.,scipy.pi/2.)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place,temp,angle=angle)
    area = [1e-3,3e-3]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperture = surface.Rect((0.,0.,0.),pos,area,vec=Vec)
    offset = geometry.Origin((0.,0.,-1.8e-2),aperture,vec=Vec)
    diodes = XTOMOchip(offset)
    diodes[0].redefine(temp)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperture)
    diodes[-1].redefine(temp)
    return diodes

def XTOMO3beam(plasma):
    temp = XTOMO3(plasma)
    output = beam.multiBeam(temp[:-1],temp[-1])
    plasma.trace(output)
    return output


#def XTOMO4(temp, place=(,0,),angle =(,,)):
#    # this is the hard part, tree storage values are not done in such a manner
#    pos = geometry.Origin(place,temp,angle=angle)
#    area = [3.6e-3,5e-4]
#    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
#    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
#    offset = goemetry.Origin((,,),aperature,Vec=Vec)
#    diodes = XTOMOchip(offset)
#    for i in diodes:
#        i.redefine(temp)
#        diodes.append(aperature)
#        diodes[-1].redefine(temp)
#    return diodes

def XTOMO5(temp, place=(.8888,0,-.2794),angle =(0.,.79587,scipy.pi/2.)):
    # this is the hard part, tree storage values are not done in such a manner
    pos = geometry.Origin(place,temp,angle=angle)
    area = [1e-3,3e-3]
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperture = surface.Rect((0.,0.,0.),pos,area,vec=Vec)
    offset = geometry.Origin((0.,0.,-2.5e-2),aperture,vec=Vec)
    diodes = XTOMOchip(offset)
    diodes[0].redefine(temp)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperture)
    diodes[-1].redefine(temp)
    return diodes

def XTOMO5beam(plasma):
    temp = XTOMO5(plasma)
    output = beam.multiBeam(temp[:-1],temp[-1])
    traceout = 2*scipy.ones((len(output),),dtype=int)
    traceout[-3:] = 0
    for i in xrange(len(output)):
        plasma.trace(output[i],
                     limiter=traceout[i])
    return output


