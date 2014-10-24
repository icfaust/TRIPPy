import scipy
import scipy.special
import matplotlib.pyplot as plt

def plotTokamak(tokamak, pltobj=None, axis=True, pargs=None, **kwargs):
    
    if pltobj is None:
        pltobj = plt

    if pargs is None:
        pltobj.plot(tokamak.sagi.s, tokamak.norm.s, **kwargs)
    else:
        pltobj.plot(tokamak.sagi.s, tokamak.norm.s, pargs, **kwargs)
        

    if axis:
        #need to correct this!!
        plt.gca().set_aspect('equal')
        pltobj.autoscale(tight=True)

def plotLine(line, invessel=True, ds=2.5e-3, pargs=None, pltobj=None, **kwargs):
    
    try:
        if invessel:
            temp = line(scipy.mgrid[line.norm.s[-2]:line.norm.s[-1]:ds])
        else:
            temp = line(scipy.mgrid[line.norm.s[0]:line.norm.s[-1]:ds])
        
        if pltobj is None:
            pltobj = plt

        if not pargs is None:
            pltobj.plot(temp.r0(), temp.x2(), pargs, **kwargs)
        else:
            pltobj.plot(temp.r0(), temp.x2(), **kwargs)

    except AttributeError:
        for i in line:
            plotLine(i, invessel=invessel, pargs=pargs, pltobj=pltobj, **kwargs)

def sinogramLine(beam, r, z, invessel=True, ds=2.5e-3, pargs=None, pltobj=None, **kwargs):
    
    try:
        if invessel:
            temp = beam(scipy.mgrid[beam.norm.s[-2]:beam.norm.s[-1]:ds])
        else:
            temp = beam(scipy.mgrid[beam.norm.s[0]:beam.norm.s[-1]:ds])

        # deal with branch cut
        temp0 = temp.t0(r, z)
        temp2 = temp.t2(r, z)
        temp = scipy.arange(temp0.size)[abs(temp2[1:] - temp2[:-1]) > scipy.pi]
        print(temp)
        if len(temp) > 0:
            temp0 = scipy.insert(temp0, temp+1, None)
            temp2 = scipy.insert(temp2, temp+1, None)

        if pltobj is None:
            pltobj = plt

        if not pargs is None:
            pltobj.plot(temp2,temp0, pargs, **kwargs)
        else:
            pltobj.plot(temp2,temp0, **kwargs)

    except AttributeError:
        for i in beam:
            sinogramLine(i, r, z, invessel=invessel, pargs=pargs, pltobj=pltobj, **kwargs)

def image(r, z, out, pltobj=None, **kwargs):
    if pltobj is None:
        pltobj = plt

    pltobj.imshow(out.T,origin='lower',extent = (r.min(),r.max(),z.min(),z.max()), **kwargs)


def plotBF(data, r, z, rcent, zcent, rmax, l=range(15), mcos=[0], msin=[], **kwargs):
    rgrid,zgrid = scipy.meshgrid(r,z)
    theta = scipy.arctan2(zgrid-zcent,rgrid-rcent)
    rgrid = scipy.sqrt((rgrid-rcent)**2 + (zgrid-zcent)**2)/rmax
    output = scipy.zeros(rgrid.shape)
    idx = 0

    u = scipy.unique(mcos+msin)
    zeros = scipy.zeros((len(u),len(l)))
    for i in xrange(len(u)):
        zeros[i] = scipy.special.jn_zeros(u[i],zeros.shape[1])

    for m in mcos:
        for i in l:
            output += data[idx]*scipy.special.jn(m,zeros[m,i]*rgrid)*scipy.cos(m*theta)
            idx += 1
    
    for m in msin:
        for i in l:
            output += data[idx]*scipy.special.jn(m,zeros[m,i]*rgrid)*scipy.sin(m*theta)
            idx += 1

    scipy.place(output,rgrid > rmax, 0)
    image(r,z,output.T,**kwargs)

def plotBFradial(data, l=range(15), mcos=[0], msin=[], err=None, **kwargs):
    rgrid = scipy.linspace(0,1,1e2)
    output = scipy.zeros(rgrid.shape)
    idx = 0
    idxe = 0
    u = scipy.unique(mcos+msin)
    zeros = scipy.zeros((len(u),len(l)))
    for i in xrange(len(u)):
        zeros[i] = scipy.special.jn_zeros(u[i],zeros.shape[1])

    for m in mcos:
        errorout = scipy.zeros(rgrid.shape)
        output = scipy.zeros(rgrid.shape)
        for i in l:
            output += data[idx]*scipy.special.jn(m,zeros[m,i]*rgrid)
            idx += 1
 
        if m > 1:
            labelin = r'$\cos$'+str(m)+r'$\theta$'
        elif m == 1:
            labelin = r'$\cos \theta$'
        else:
            labelin = r'radial'
            
        plt.plot(rgrid,
                 output,
                 label=labelin,
                 **kwargs)
        
        if not err is None:
            outpute= scipy.zeros(rgrid.shape)
            for i in l:
                outpute += err[idxe]*scipy.special.jn(m,zeros[m,i]*rgrid)
                idxe += 1
                
            plt.fill_between(rgrid,output-outpute,output+outpute,color='k',alpha=.3)



    
    for m in msin:
        output = scipy.zeros(rgrid.shape)
        for i in l:
            output += data[idx]*scipy.special.jn(m,zeros[m,i]*rgrid)
            idx += 1
        if m > 1:
            labelin =r'$\sin$'+str(m)+r'$\theta$'
        else:
            labelin =r'$\sin \theta$'
       
        plt.plot(rgrid,
                 output,
                 label=labelin,
                 **kwargs)

        if not err is None:
            outpute= scipy.zeros(rgrid.shape)
            for i in l:
                outpute += err[idxe]*scipy.special.jn(m,zeros[m,i]*rgrid)
                idxe += 1
                
            plt.fill_between(rgrid,output-outpute,output+outpute,color='k',alpha=.3)

def plotBFbright(sens,bright,prof,beams):
    temp = scipy.dot(sens,prof).T[0:len(bright)]

    for i in xrange(len(temp)):
        temp[i] *= beams[i].etendue/(4*scipy.pi)

    plt.plot(scipy.arange(len(temp)),temp,marker='s',color='k',linestyle=' ',label='reconstruction')
    plt.plot(scipy.arange(len(bright)),bright,marker='$\circ$',linestyle=' ',label='model')

    plt.ylabel(r'Current [$\mu$A]')
    plt.xlabel(r'Chord Number')
    plt.ylim((0,8))
    plt.title(r'Chord Signals')
    plt.text(8.5,1,r'(1)',size=20)
    plt.text(28.5,1,r'(2)',size=20)    
    plt.text(48.5,1,r'(3)',size=20)
    plt.gca().axes.get_xaxis().set_ticks_position('bottom')
    plt.gca().axes.get_yaxis().set_ticks_position('left')
    plt.gca().legend(loc='upper right',numpoints=1)

def test():
    print('nothing')
