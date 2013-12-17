import plasma
import beam
import scipy
import scipy.interpolate
import matplotlib.pyplot as plt


def sens(beams,plasmameth,time,points,step=1e-3):
    """ optimal to give multiple times """
    time = scipy.atleast_1d(time)
    interp = scipy.interpolate.interp1d(points,scipy.arange(len(points)),kind='cubic')
    # initialize output array of sensitivities
    output = scipy.zeros((len(time),len(beams),len(points)))

    for i in xrange(len(beams)):
        temp = beams[i].norm.s      
        beams[i].norm.s = scipy.mgrid[beams[i].norm.s[-2]:beams[i].norm.s[-1]:step]
        mapped = plasmameth(beams[i].r()[0],beams[i].r()[2],time)
 
        # knowing that the last point (point[-1]) is assumed to be a ZERO emissivity point,
        #an additional brightness is added which only sees the last emissivity to yield the zero
        scipy.place(mapped,mapped > points[-1], points[-1])
       
        # for a given point along a chord, use a spline to solve what reconstruction points
        #it is most close to. Then in the weighting matrix, (which is (brightness,points) in shape add the various fractional weighting.
        out = interp(mapped)
        
        # find out point using a floor like command (returns ints) 
        idx = out.astype(int)

        # reduce out to the fraction
        out = (out % 1.)*step
        for j in range(len(idx[0])):
            output[:,i,idx[:,j]] += out[:,j]
            scipy.place(out[:,j], out[:,j] == len(points), len(points) - 1)
            output[:,i,idx[:,j]+1] += step - out[:,j]
        beams[i].norm.s = temp      

    return output

def rhosens(beams,plasma,time,points,step=1e-3,meth='psinorm'):
    """ optimal to give multiple times """
    time = scipy.atleast_1d(time)
    interp = scipy.interpolate.interp1d(points,scipy.arange(len(points)),kind='cubic')
    # initialize output array of sensitivities
    output = scipy.zeros((len(time),len(beams),len(points)))

    for i in xrange(len(beams)):
        temp = beams[i].norm.s
        beams[i].norm.s = scipy.mgrid[beams[i].norm.s[-2]:beams[i].norm.s[-1]:step]
        mapped = plasma.eq.rz2rho(meth,beams[i].r()[0],beams[i].r()[2],scipy.array(time))

        # knowing that the last point (point[-1]) is assumed to be a ZERO emissivity point,
        #an additional brightness is added which only sees the last emissivity to yield the zero
        scipy.place(mapped,mapped > points[-1], points[-1])
       
        # for a given point along a chord, use a spline to solve what reconstruction points
        #it is most close to. Then in the weighting matrix, (which is (brightness,points) in shape add the various fractional weighting.
        out = interp(mapped)
        
        # find out point using a floor like command (returns ints) 
        idx = out.astype(int)

        # reduce out to the fraction
        out = (out % 1.)*step
        for j in range(len(idx[0])):
            output[:,i,idx[:,j]] += out[:,j]
            scipy.place(idx[:,j], idx[:,j] < len(points)-1, idx[:,j]+1)
            output[:,i,idx[:,j]] += step - out[:,j]
        beams[i].norm.s = temp

    return output
