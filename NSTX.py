import geometry,surface,AXUV
import scipy

def diode1(temp,place=(2.0,0,2e-2), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle,flag=False)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815,0,.5))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-2.25e-2),aperature,Vec=Vec)
    diodes = AXUV.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes

def diode2(temp,place=(2.0,0,-2e-2), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle,flag=False)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815,0,-.5))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-2.25e-2),aperature,Vec=Vec)
    diodes = AXUV.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes

def diode3(temp,place=(.8,0,1.8), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle,flag=False)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815-.075,0,0))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-1.25e-2),aperature,Vec=Vec)
    diodes = AXUV.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes

def diode4(temp,place=(1.5,0,-1.5), angle=(0,-scipy.pi/2,0)):
    pos = geometry.Origin(place,temp,angle=angle,flag=False)
    area = [4e-3,3e-3]
    centervec = geometry.Vecx(scipy.array((.815-.4,0,-1.175))-place)
    # point toward the nominal center
    Vec = [centervec,geometry.cross(geometry.Vecx((0,1,0)),centervec)]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec)
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0,0.,-2.5e-2),aperature,Vec=Vec)
    diodes = AXUV.AXUV20(offset)
    for i in diodes:
        i.redefine(temp)
    diodes.append(aperature)
    diodes[-1].redefine(temp)
    return diodes

def diode5(center,ap=(0,0,5e-2)):
    pt1 = geometry.Point(scipy.array((84.98,-67.48,0))*2.54e-2,center)
    pt2 = geometry.Point(scipy.array((82.5,-69.44,0))*2.54e-2,center)
    port = geometry.Point(scipy.array((73.0681,-19.3315,0))*2.54e-2,center)
    
    
    d = geometry.pts2Vec(pt2,pt1)
    d.s = d.s/2
    dnorm = geometry.cross(d,geometry.Vecx((0,0,1)))
    
    temp = geometry.Point((pt2.vec + d).x(),center)
    #12 ports
    angle = scipy.pi*4/12 - geometry.angle(port.vec,temp.vec)
    diodeangle = geometry.angle(dnorm,temp.vec)

    cent = temp.vec.c()
    cent.unit[1] = angle
    pos = geometry.Point(cent.c().x(),center) #redefine the position relative to G Side detectors in cartesian coords
    cent.unit[1] = cent.unit[1]-diodeangle+scipy.pi # look in opposite direction

    Vec = [cent.c(),geometry.Vecx((0,0,1))]
    # need to clean this all up and put in unary minus dear lord

    diodeorigin = geometry.Origin(pos.x(),center,Vec=Vec)
    rot = geometry.Origin((0,0,0),diodeorigin,angle=(0,0,0))#angle=(scipy.pi/2,0,scipy.pi/2))
    diodes = AXUV.AXUV20(rot)

    area = [4e-3,3e-3]
    aperature = surface.Rect(ap,rot,area,angle=(0,0,0))
    for i in diodes:
        i.redefine(center)
    diodes.append(aperature)
    diodes[-1].redefine(center)
    return diodes

def diode6(center,ap=(0,0,5e-2)):
    pt1 = geometry.Point(scipy.array((80.01,-71.29,0))*2.54e-2,center)
    pt2 = geometry.Point(scipy.array((77.53,-73.18,0))*2.54e-2,center)
    port = geometry.Point(scipy.array((73.0681,-19.3315,0))*2.54e-2,center)
    
    
    d = geometry.pts2Vec(pt2,pt1)
    d.s = d.s/2
    dnorm = geometry.cross(d,geometry.Vecx((0,0,1)))
    
    temp = geometry.Point((pt2.vec + d).x(),center)
    #12 ports
    angle = scipy.pi*4/12 - geometry.angle(port.vec,temp.vec)
    diodeangle = geometry.angle(dnorm,temp.vec)
    print(diodeangle*180/scipy.pi)
    cent = temp.vec.c()
    print(cent.c().x())
    cent.unit[1] = angle
    print(cent.c().x()) 
    pos = geometry.Point(cent.c().x(),center) #redefine the position relative to G Side detectors in cartesian coords
    cent.unit[1] += scipy.pi +diodeangle # look in opposite direction

    Vec = [cent.c(),geometry.Vecx((0,0,1))]
    # need to clean this all up and put in unary minus dear lord

    diodeorigin = geometry.Origin(pos.x(),center,Vec=Vec)
    rot = geometry.Origin((0,0,0),diodeorigin,angle=(0,0,0))#(scipy.pi/2,0,scipy.pi/2))
    diodes = AXUV.AXUV20(rot)

    area = [4e-3,3e-3]
    aperature = surface.Rect(ap,rot,area,angle=(0,0,0))
    for i in diodes:
        i.redefine(center)
    diodes.append(aperature)
    diodes[-1].redefine(center)
    print(diodes[0].x())
    print(diodes[-1].x())
    return diodes
