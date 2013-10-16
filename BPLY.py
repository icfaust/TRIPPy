import AXUV
import surface, geometry, scipy
import beam as beamin
import scipy.interpolate

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

def getBeamFluxSpline(beam,plasma,t,points = 250):
    """ generates a spline off of the beampath.  Assumes
    that the change in flux is MONOTONIC"""

    lim = beam.norm.s
    beam.norm.s = scipy.linspace(lim[-2],lim[-1],points)
    psi = plasma.eq.rz2rmid(beam.x()[0],beam.x()[2],t)
    minpos = psi.argmin()
    lim = scipy.insert(lim,(-2,-2),(beam.norm.s[minpos],beam.norm.s[minpos]))  # add minimum flux s for bound testing
    outspline = scipy.interpolate.interp1d(psi[0:minpos],
                                           beam.norm.s[0:minpos],
                                           kind = 'cubic',
                                           bounds_error = False)
    inspline = scipy.interpolate.interp1d(psi[minpos:],
                                          beam.norm.s[minpos:],
                                          kind = 'cubic',
                                          bounds_error = False)
    beam.norm.s = lim.remove(0)
    return (outspline,inspline)

def calcArea(lower,upper):
    """ it is assumed that the inner (rmid < lcfs) point can not be calculated.
    thus, initially the value is set to the split (inner vs outer SOL position)
    this will increase the size of the polygon while not increasing the area"""
    

def Area(points):
    val = 0 
    for i in scipy.arange(len(points))-1:
        val += points[i].x0()*points[i+1].x2() - points[i].x2()*points[i+1].x0()
    return val/2

def viewPoints(surf1,surf2,plasma,t,lim1 = .88,lim2 = .92):

    beam = beamin.Beam(surf1,surf2)
    ray1 = beamin.Ray(surf1.edge().split(plasma)[0][0],surf1.edge().split(plasma)[1][1])
    ray2 = beamin.Ray(surf1.edge().split(plasma)[1][0],surf1.edge().split(plasma)[0][1])
    beam.trace(plasma)
    ray1.trace(plasma)
    ray2.trace(plasma)
    blim = beam.norm.s
    r1lim = ray1.norm.s
    r2lim = ray2.norm.s

    output = t.size * [0]

    for i in xrange(t.size):
        
        beam.norm.s = blim
        ray1.norm.s = r1lim
        ray2.norm.s = r2lim

        outermid,innermid = getBeamFluxSpline(beam,plasma,t[i])
        outertop,innertop = getBeamFluxSpline(ray1,plasma,t[i])
        outerbot,innerbot = getBeamFluxSpline(ray2,plasma,t[i])
        segments = 3*[0]
        #beam and ray masking values are already written to their norm.s values
        segment[0] = scipy.array([outertop((lim2,lim1)),innertop((lim1,lim2))]).flatten()
        segment[1] = scipy.array([outermid((lim2,lim1)),innermid((lim1,lim2))]).flatten()
        segment[2] = scipy.array([outerbot((lim2,lim1)),innerbot((lim1,lim2))]).flatten()
        
        # compare and mask/replace
        segment[0] = scipy.replace(segment[0],segment[0] == scipy.nan, ray1.norm.s[1:])
        segment[1] = scipy.replace(segment[1],segment[1] == scipy.nan, beam.norm.s[1:])
        segment[2] = scipy.replace(segment[2],segment[2] == scipy.nan, ray2.norm.s[1:])
        
        #turn into points
        ray1.norm.s = segment[0]
        beam.norm.s = segment[1]
        ray2.norm.s = segment[2]

        

        output[i] = [ray1.split(plasma,obj=geometry.Point),beam.split(plasma,obj=geometry.Point),ray2.split(plasma,obj=geometry.Point)]


    return output


def effectiveHeight(surf1, surf2, plasma, t, lim1=.88, lim2=.92):
    output = (calcArea(segments[0]) + calcArea(segments[1]))/(inlen.s+outlen.s)               
    return output
