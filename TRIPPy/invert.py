import plasma
import beam
import _beam
import scipy
import scipy.interpolate
import scipy.integrate
import scipy.special
import matplotlib.pyplot as plt
import warnings


# pixel methods using flux surfaces (NEEDS TO BE CLEANED UP MASSIVELY)

def fluxFourierSens(beam, plasmameth, centermeth, time, points, mcos=[0], msin=[], ds=1e-3):
    """ optimal to give multiple times """
    time = scipy.atleast_1d(time)
    interp = scipy.interpolate.interp1d(points,scipy.arange(len(points)),kind='cubic')
    # initialize output array of sensitivities
    m = scipy.unique(mcos+msin)
        
    try:
        
        output = scipy.zeros((len(time),len(points))) 
        temp = beams(scipy.mgrid[beam.norm.s[-2]:beams.norm.s[-1]:ds])
        mapped = plasmameth(temp.r0(),
                            temp.x2(),
                            time)
 
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
            output[:,idx[:,j]] += out[:,j]
            scipy.place(out[:,j], out[:,j] == len(points), len(points) - 1)
            output[:,idx[:,j]+1] += ds - out[:,j]

    except AttributeError:
        output = scipy.zeros((len(time),len(beams),len(points)))
        for i in xrange(len(beams)):
            output[:,i,:] = fluxFourierSens(beam[i],
                                            plasmameth,
                                            time,
                                            points,
                                            mcos=mcos,
                                            msin=msin,
                                            ds=ds)

    return output

def fluxFourierSensRho(beams,plasma,time,points,ds=1e-3,meth='psinorm'):
    """ optimal to give multiple times """
    time = scipy.atleast_1d(time)
    interp = scipy.interpolate.interp1d(points,scipy.arange(len(points)),kind='cubic')
    
    # initialize output array of sensitivities
    output = scipy.zeros((len(time),len(beams),len(points)))

    for i in xrange(len(beams)):
        temp = beams[i](scipy.mgrid[beams[i].norm.s[-2]:beams[i].norm.s[-1]:step])
        mapped = plasma.eq.rz2rho(meth,
                                  temp.r0(),
                                  temp.x2(),
                                  scipy.array(time))

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

    return output

def besselFourierKernel(m, zero, rho):
    # I SHOULD TRY AND VECTORIZE THIS AS MUCH AS POSSIBLE
    jprime = (scipy.special.jn(m+1, zero) - scipy.special.jn(m-1, zero))
    return jprime*scipy.integrate.quad(_beam.bessel_fourier_kernel,
                                       0,
                                       scipy.arccos(rho),
                                       args = [m, zero, rho])

def besselFourierSens(beam, rcent, zcent, l=range(15), mcos=[0], msin=[], rmax):

    # find and store bessel zeros 
    m = scipy.unique(mcos+msin)
    length = len(l)
    zeros = scipy.zeros((len(m),length))
    for i in xrange(len(m)):
        zeros[i] = scipy.special.jn_zeros(i,zeros.shape[1])

    kernel = scipy.zeros(len(m),length)

    # need to modify zeros

    try:
        try:
            output = scipy.zeros((len(rcent), length*(len(mcos)+len(msin))))
        except TypeError:
            output = scipy.zeros((1, length*(len(mcos)+len(msin))))
            rmax = scipy.atleast_1d(rmax)

        for i in xrange(len(rcent)):
            idx = 0
                # returns closest approach vector to plasma center
            temp = beam(beam.tmin(rcent[i], zcent[i])).t(rcent[i], zcent[i])

            if temp[0] > rmax[i]:
                rho = 1.0
                warnings.warn('chord outside of specified designated edge zero emissivity', RuntimeWarning)
            else:
                rho = temp[0]/rmax[i]

            for j in xrange(len(m)):
                for k in xrange(length):
                    kernel[j, k] = rmax[i]*besselFourierKernel(m[j],
                                                                 zeros[j,k],
                                                                 rho)
            # fill sens matrix
            for j in xrange(len(mcos)):
                output[i, idx*length:(idx+1)*length] = scipy.cos(mcos[j]*temp[2])*kernel[scipy.where(m == mcos[j])]
                idx++

            for j in xrange(len(msin)):     
                output[i, idx*length:(idx+1)*length] = scipy.sin(msin[j]*temp[2])*kernel[scipy.where(m == msin[j])]
                idx++
                    
        return output

    except AttributeError:
        #implement parallelization here
        output = scipy.zeros((len(rcent),
                              len(beam),
                              length*len(mcos+msin)))
        
        for i in xrange(len(beam)):
            output[:,i,:] += [besselFourierSens(beam[i],
                                                rcent,
                                                zcent,
                                                l=l,
                                                mcos=mcos,
                                                msin=msin)]
        return output
