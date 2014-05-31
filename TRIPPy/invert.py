import plasma
import beam
import _beam
import geometry
import scipy
import scipy.interpolate
import scipy.integrate
import scipy.special
import scipy.linalg
import matplotlib.pyplot as plt
import warnings

def fluxFourierSens(beam, plasmameth, centermeth, time, points, mcos=[0], msin=[], ds=1e-3):
    """ optimal to give multiple times """
    time = scipy.atleast_1d(time)
    interp = scipy.interpolate.interp1d(points,
                                        scipy.arange(len(points)),
                                        kind='cubic')
    length = len(points)
        
    try:
        
        output = scipy.zeros((len(time),length*len(mcos+msin))) 
        temp = beams(scipy.mgrid[beam.norm.s[-2]:beams.norm.s[-1]:ds])
        mapped = plasmameth(temp.r0(),
                            temp.x2(),
                            time)

        # recover angles of each position in temp vector utlizing the t2 method
        # to geometry.Vec improper vectorization strategy in t2 causes the use
        # of a for loop
        angle = scipy.zeros(mapped.shape)
        for i in xrange(len(time)):
            pt0 = centermeth(time[i])
            angle[i] = temp.t2(pt0[0],pt0[1])
        
        # knowing that the last point (point[-1]) is assumed to be a ZERO 
        # emissivity point an additional brightness is added which only sees 
        # the last emissivity to yield the zero
        scipy.place(mapped, mapped > points[-1], points[-1])

        if mapped.min() < 0:
            warning.warn('chord measures at a parameter below point grid', RuntimeWarning)
        # for a given point along a chord, use a spline to solve what 
        # reconstruction points it is most close to. Then in the weighting
        # matrix, (which is (brightness,points) in shape add the various 
        # fractional weighting.
        out = interp(mapped)
        
        # find out point using a floor like command (returns ints) 
        idx1 = out.astype(int)
        scipy.clip(idx1, 0, length-1, out=idx1)
        
        idx2 = idx1 + 1
        scipy.clip(idx2, 0, length-1, out=idx2)
        # reduce out to the fraction in nearby bins
        out = (out % 1.)*ds
        lim = 0

        for i in mcos:
            for j in range(idx1.shape[1]):
                output[:,lim+idx1[:,j]] += out[:,j]*scipy.cos(i*angle)
                output[:,lim+idx2[:,j]] += (ds - out[:,j])*scipy.cos(i*angle)
                lim += length

        for i in msin:
            for j in range(idx2.shape[1]):
                output[:,lim+idx1[:,j]] += out[:,j]*scipy.sin(i*angle)
                output[:,lim+idx2[:,j]] += (ds - out[:,j])*scipy.sin(i*angle)
                lim += length

    except AttributeError:
        output = scipy.zeros((len(time),len(beams),length*len(mcos+msin)))
        for i in xrange(len(beams)):
            output[:,i,:] = fluxFourierSens(beam[i],
                                            plasmameth,
                                            centermeth,
                                            time,
                                            points,
                                            mcos=mcos,
                                            msin=msin,
                                            ds=ds)

    return output

def fluxFourierSensRho(beams,plasma,time,points,ds=1e-3,meth='psinorm'):
    """ optimal to give multiple times, NEEDS TO BE REWORKED AFTER SIN/COS update """
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
                                       args = (m, zero, rho))[0]

def bessel_fourier_kernel(theta,m,zero,rho):
    return scipy.cos(m*theta)*scipy.sin(zero*(scipy.cos(theta)-rho))

def besselFourierSens(beam, rcent, zcent, rmax, l=range(15), mcos=[0], msin=[], rcond=2e-2):

    # find and store bessel zeros 
    m = scipy.unique(mcos+msin)
    length = len(l)
    zeros = scipy.zeros((len(m),length))
    for i in xrange(len(m)):
        zeros[i] = scipy.special.jn_zeros(m[i],zeros.shape[1])

    kernel = scipy.zeros((len(m),length))

    # need to modify zeros

    try:
        try:
            output = scipy.zeros((len(rcent), length*(len(mcos)+len(msin))))
        except TypeError:
            output = scipy.zeros((1, length*(len(mcos)+len(msin))))
            rmax = scipy.atleast_1d(rmax)

        for i in xrange(len(rcent)):
            idx = 0
            
            # returns closest approach vector to plasma center mod for NSTX-U
            
            #temp = j(j.tmin(rcent,zcent)).t(rcent,zcent)
            cent = geometry.Point(geometry.Vecr([rcent[i],0,zcent[i]]),beam._origin)
            temp2 = beam(beam.smin(cent)) - cent
            temp = [temp2.s,0,scipy.arctan2(temp2.x2(),temp2.x0())]
            
            #temp = beam(beam.tmin(rcent[i], zcent[i])).t(rcent[i], zcent[i])

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
                idx += 1

            for j in xrange(len(msin)):     
                output[i, idx*length:(idx+1)*length] = scipy.sin(msin[j]*temp[2])*kernel[scipy.where(m == msin[j])]
                idx += 1
                    
        return output

    except AttributeError:
        #implement parallelization here
        output = scipy.zeros((len(rcent),
                              len(beam),
                              length*len(mcos+msin)))
        
        for i in xrange(len(beam)):
            output[:,i,:] += besselFourierSens(beam[i],
                                               rcent,
                                               zcent,
                                               rmax,
                                               l=l,
                                               mcos=mcos,
                                               msin=msin)

        return output

    
def bFInvert(beams, bright, rcent, zcent, rmax, l=range(15), mcos=[0], msin=[], zeros=None, plasma=None, rcond=2e-2, out=False):
    """ zeros adds the capability of providing guaranteed edge zeros"""
    bright = bright*4*scipy.pi
    for i in xrange(len(bright)):
        bright[i] = bright[i]/beams[i].etendue
    
    if (not plasma is None) and (not zeros is None):
        beams += beam._genBFEdgeZero(plasma, zeros, rcent, zcent)
        bright = scipy.concatenate((bright,scipy.zeros((zeros,)))) #added zeros to bright
    
    sens = besselFourierSens(beams, rcent, zcent, rmax, l=l, mcos=mcos, msin=msin, rcond=rcond)
    output = scipy.zeros((len(sens),len(l)*len(mcos+msin)))
    for i in xrange(len(sens)):
        output[i] = scipy.dot(scipy.linalg.pinv(sens[i],rcond=rcond),bright)
        
    if out:
        return output,sens
    else:
        return output

def cov(sens):
    return scipy.linalg.pinv(scipy.dot(sens.T,sens),rcond=2e-2)

def err(emiss, bright, sens, beams, num=None):
    """ returns the error of emiss"""
    temp2 = bright[:]
    for i in xrange(len(temp2)):
        temp2[i] *= 4*scipy.pi/beams[i].etendue
    temp = scipy.sum((scipy.dot(sens,emiss)[0:len(bright)]-temp2)**2)
    var = cov(sens)
    if num is None:
        output = scipy.zeros(emiss.shape)
    else:
        output = scipy.zeros((num,))
    
    for i in xrange(len(output)):
        output[i] = var[i,i]

    return scipy.sqrt(abs(temp*output/(bright.size-emiss.size)))
