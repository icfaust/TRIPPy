import AXUV
import surface, geometry, scipy
import beam as beamin
import scipy.interpolate
import matplotlib.pyplot as plt
import time

def BPLY(temp, place=(1.87,0,.157277), angle=(0,.17453+scipy.pi/2,-1.62385+scipy.pi/2)):


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


def BPLYbeam(alcator):

    temp = BPLY(alcator)
    output = 22*[0]
    for i in xrange(len(output)):
        output[i] = beam.Beam(temp[i],temp[-1])
        output[i].trace(alcator)
    return output

def getBeamFluxSpline(beam,plasma,t,lim1,lim2,points = 1000):
    """ generates a spline off of the beampath.  Assumes
    that the change in flux is MONOTONIC"""

    lim = beam.norm.s

    beam.norm.s = scipy.linspace(0,lim[-1],points)
    
    psi = plasma.eq.rz2rmid(beam.x()[0],beam.x()[2],t) #evaluates all psi's at once
    outspline = len(t)*[0]
    inspline = len(t)*[0]
    for i in xrange(t.size):
        mask = scipy.isfinite(psi[i])
        minpos = scipy.argmin(psi[i][mask])
        sizer = psi[i][mask].size
    #plt.plot(beam.x()[0][mask][0:minpos],psi[mask][0:minpos],beam.x()[0][mask][minpos:],psi[mask][minpos:])
    #plt.show()
        #limout = scipy.insert(lim,(2,2),(beam.norm.s[mask][minpos],beam.norm.s[mask][minpos]))  # add minimum flux s for bound testing
        outspline[i] = beam.norm.s[mask][minpos::-1][scipy.clip(scipy.digitize((lim1,lim2),psi[i][mask][minpos::-1]),0,minpos)]
        inspline[i] = beam.norm.s[mask][minpos:][scipy.clip(scipy.digitize((lim1,lim2),psi[i][mask][minpos:]),0,sizer-minpos-1)]
    #outspline = scipy.interpolate.interp1d(psi[i][mask][minpos::-1],
    #                                       beam.norm.s[mask][minpos::-1],
    #                                       bounds_error = False)((lim1,lim2))
    #inspline = scipy.interpolate.interp1d(psi[i][mask][minpos:],
    #                                      beam.norm.s[mask][minpos:],
    #                                      bounds_error = False)((lim1,lim2))



    return (outspline,inspline)

def calcArea(points):
    val = 0 
    for i in scipy.arange(len(points))-1:
        val += points[i].x0()*points[i+1].x2() - points[i].x2()*points[i+1].x0()
    return val/2

def viewPoints(surf1,surf2,plasma,t,lim1 = .88,lim2 = .92,fillorder = True):
    h=time.time()
    beam = beamin.Beam(surf1,surf2)
    ray1 = beamin.Ray(surf1.edge().split(plasma)[0][0],surf2.edge().split(plasma)[1][1])
    ray2 = beamin.Ray(surf1.edge().split(plasma)[1][1],surf2.edge().split(plasma)[0][0])
    print(time.time()-h,'zero')
    h= time.time()
    beam.trace(plasma,step=1e-3)
    ray1.trace(plasma,step=1e-3) #there has to be a way to improve this /only calculate this once
    ray2.trace(plasma,step=1e-3)
    blim = beam.norm.s
    r1lim = ray1.norm.s
    r2lim = ray2.norm.s

    
    output = t.size * [0]
    outermid,innermid = getBeamFluxSpline(beam,plasma,t,lim1,lim2)
    outertop,innertop = getBeamFluxSpline(ray1,plasma,t,lim1,lim2)
    outerbot,innerbot = getBeamFluxSpline(ray2,plasma,t,lim1,lim2)

    #condition inputs for area calculations
    for i in xrange(t.size):
        print(time.time()-h,'two')
        h= time.time()
        segment = 3*[0]
        #beam and ray masking values are already written to their norm.s values

        segment[0] = scipy.array([outertop[i][1],outertop[i][0],innertop[i][0],innertop[i][1]])
        segment[1] = scipy.array([outermid[i][1],outermid[i][0],innermid[i][0],innermid[i][1]])
        segment[2] = scipy.array([outerbot[i][1],outerbot[i][0],innerbot[i][0],innerbot[i][1]])
        print(time.time()-h,'three')
        h= time.time()
        # compare and mask/replace
        #scipy.copyto(segment[0],ray1.norm.s,where=~scipy.isfinite(segment[0]))
        #scipy.copyto(segment[1],beam.norm.s,where=~scipy.isfinite(segment[1]))
       # scipy.copyto(segment[2],ray2.norm.s,where=~scipy.isfinite(segment[2]))
 
        #turn into points
        ray1.norm.s = segment[0]
        beam.norm.s = segment[1]
        ray2.norm.s = segment[2]
        print(time.time()-h,'four')
        h= time.time()

        temp1 = ray1.split(plasma,obj=geometry.Point)
        temp2 = beam.split(plasma,obj=geometry.Point)
        temp3 = ray2.split(plasma,obj=geometry.Point)

        print(time.time()-h,'five')
        h= time.time()
        
        if fillorder:
            output[i] = []
            for j in ((0,1,1,1,0,0),(2,3,3,3,2,2)):
                output[i] += [scipy.array([temp1[j[0]].x(),
                                           temp1[j[1]].x(),
                                           temp2[j[2]].x(),
                                           temp3[j[3]].x(),
                                           temp3[j[4]].x(),
                                           temp2[j[5]].x()]).T]
        else:
            output[i] = [temp1,temp2,temp3]
    return output
    





def effectiveHeight(surf1, surf2, plasma, t, lim1=.88, lim2=.92):
    output = (calcArea(segments[0]) + calcArea(segments[1]))/(inlen.s+outlen.s)               
    return output
