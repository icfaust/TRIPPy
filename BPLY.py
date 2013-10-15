import AXUV
import surface, geometry, scipy, beam
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
    lim = scipy.insert(lim,-2,beam.norm.s[minpos])  # add minimum flux s for bound testing
    outspline = scipy.interpolate.interp1d(psi[0:minpos],
                                           beam.norm.s[0:minpos],
                                           kind = 'cubic')
    inspline = scipy.interpolate.interp1d(psi[minpos:],
                                          beam.norm.s[minpos:],
                                          kind = 'cubic')
    beam.norm.s = lim
    return (outspline,inspline)

def calcArea(lower,upper):
""" it is assumed that the inner (rmid < lcfs) point can not be calculated.
thus, initially the value is set to the split (inner vs outer SOL position)
this will increase the size of the polygon while not increasing the area"""
    

def effectiveHeight(beam,ray1,ray2,plasma,t):

    outermid,innermid = getBeamFluxSpline(beam,plasma,t)
    outertop,innertop = getBeamFluxSpline(ray1,plasma,t)
    outerbot,innerbot = getBeamFluxSpline(ray2,plasma,t)
    return (calcArea(outerbot,outertop) + calcArea(innerbot,innertop))/(inlen.s+outlen.s)
    
