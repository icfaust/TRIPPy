import AXUV
import surface, geometry, scipy
import beam as beamin
import scipy.interpolate
import matplotlib.pyplot as plt
import time
import multiprocessing as multiproc
import MDSplus
import gc
import plasma, eqtools

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
    h = time.time()
    psi = plasma.eq.rz2rmid(beam.x()[0],beam.x()[2],t) #evaluates all psi's at once
    print(time.time()-h)
    outspline = len(t)*[0]
    inspline = len(t)*[0]
    for i in xrange(t.size):
        mask = scipy.logical_and(scipy.isfinite(psi[i]),psi[i] < lim2+.02)
        try:
            minpos = scipy.argmin(psi[i][mask])
        except ValueError:
            return ((0,0),(0,0))
        sizer = psi[i][mask].size
    #plt.plot(beam.x()[0][mask][0:minpos],psi[mask][0:minpos],beam.x()[0][mask][minpos:],psi[mask][minpos:])
    #plt.show()
        #limout = scipy.insert(lim,(2,2),(beam.norm.s[mask][minpos],beam.norm.s[mask][minpos]))  # add minimum flux s for bound testing

        try:
            temp1 = scipy.clip(scipy.digitize((lim1,lim2),psi[i][mask][minpos::-1]),0,minpos)
            outspline[i] = beam.norm.s[mask][minpos::-1][temp1]
        except ValueError:
           
            #plt.plot(beam.x()[0][mask][minpos:maxpos:-1],psi[i][mask][minpos:maxpos:-1])
            #plt.plot(beam.x()[0][mask][minpos:],psi[i][mask][minpos:])
            #plt.show()
            outspline[i] = scipy.array(2*[beam.norm.s[mask][minpos]])
        except IndexError:
            plt.plot(beam.x()[0][mask],psi[i][mask])
            plt.show()

        try:
            temp2 = scipy.clip(scipy.digitize((lim1,lim2),psi[i][mask][minpos:]),0,sizer-minpos-1)
            inspline[i] = beam.norm.s[mask][minpos:][temp2]
        except ValueError:
            inspline[i] = scipy.array([beam.norm.s[mask][minpos],beam.norm.s[mask][-1]])
    


    #outspline = scipy.interpolate.interp1d(psi[i][mask][minpos::-1],
    #                                       beam.norm.s[mask][minpos::-1],
    #                                       bounds_error = False)((lim1,lim2))
    #inspline = scipy.interpolate.interp1d(psi[i][mask][minpos:],
    #                                      beam.norm.s[mask][minpos:],
    #                                      bounds_error = False)((lim1,lim2))



    return (outspline,inspline)

def calcArea(points):
    val = 0 
    for i in scipy.arange(len(points[0]))-1:
        val += points[0][i]*points[2][i+1] - points[2][i]*points[0][i+1]
    return val/2

def viewPoints(surf1,surf2,plasma,t,lim1 = .88,lim2 = .92,fillorder = True):
    h=time.time()
    beam = beamin.Beam(surf1,surf2)
    ray1 = beamin.Ray(surf1.edge().split(plasma)[0][0],surf2.edge().split(plasma)[1][1])
    ray2 = beamin.Ray(surf1.edge().split(plasma)[1][1],surf2.edge().split(plasma)[0][0])
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

        segment = 3*[0]
        #beam and ray masking values are already written to their norm.s values
        segment[0] = scipy.array([outertop[i][1],outertop[i][0],innertop[i][0],innertop[i][1]])
        segment[1] = scipy.array([outermid[i][1],outermid[i][0],innermid[i][0],innermid[i][1]])
        segment[2] = scipy.array([outerbot[i][1],outerbot[i][0],innerbot[i][0],innerbot[i][1]])

        # compare and mask/replace
        #scipy.copyto(segment[0],ray1.norm.s,where=~scipy.isfinite(segment[0]))
        #scipy.copyto(segment[1],beam.norm.s,where=~scipy.isfinite(segment[1]))
       # scipy.copyto(segment[2],ray2.norm.s,where=~scipy.isfinite(segment[2]))
 
        #turn into points
        ray1.norm.s = segment[0]
        beam.norm.s = segment[1]
        ray2.norm.s = segment[2]

        temp1 = ray1.split(plasma,obj=geometry.Point)
        temp2 = beam.split(plasma,obj=geometry.Point)
        temp3 = ray2.split(plasma,obj=geometry.Point)

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
    """ calculate the effective height of a view through the scrape-off-layer"""

    segments = viewPoints(surf1,surf2,plasma,t,lim1=lim1,lim2=lim2,fillorder=False)
    output = scipy.zeros((len(segments),))
    for i in xrange(len(segments)):
        inlen = geometry.pts2Vec(segments[i][1][0],segments[i][1][1])
        outlen = geometry.pts2Vec(segments[i][1][2],segments[i][1][3])
        temp = []
        for j in ((0,1,1,1,0,0),(2,3,3,3,2,2)):
            temp += [scipy.array([segments[i][0][j[0]].x(),
                                  segments[i][0][j[1]].x(),
                                  segments[i][1][j[2]].x(),
                                  segments[i][2][j[3]].x(),
                                  segments[i][2][j[4]].x(),
                                  segments[i][1][j[5]].x()]).T]

        output[i] = (calcArea(temp[0]) + calcArea(temp[1]))/(inlen.s + outlen.s)               
    return output


def writeToTree(inp):
    shot,chan,surf1,surf2 = inp
    timeout = time.time()
    outstr = str(chan)
    if chan < 10:
        outstr = '0' + outstr
        
    Tree = MDSplus.Tree('spectroscopy',shot)
    b = plasma.Tokamak(eqtools.CModEFITTree(shot,tspline=True)) # I HATE THIS
    output = effectiveHeight(surf1,surf2,b,scipy.mgrid[.2:1.7:1e-4])
    dummy =  MDSplus.Data.compile('build_signal($1,*,$2)',output,scipy.mgrid[.2:1.7:1e-4])
    Tree.getNode('.BOLOMETER.RESULTS.DIODE.BPLY.AREA:CHORD_'+outstr).putData(dummy)
    print('channel '+str(chan)+':\t '+str(time.time()-timeout))


def calctoTree(shot):


    cores = multiproc.cpu_count()
    cores = scipy.array([cores,22]).min() # use up to 32 cores
    print(shot)
    print(cores)
    k = geometry.Center()
    n = BPLY(k)
    pool = multiproc.Pool(cores)
    
    try: 
        print('beginning map.')
        out = pool.map(writeToTree,[[shot,22-i,n[i],n[-1]] for i in scipy.arange(21)])
        pool.close()
        print(shot)

    except MDSplus._treeshr.TreeException:
        print(shot)
        gc.collect()

def globalpowerCalc(shot):

    Tree = MDSplus.Tree('spectroscopy',shot)
    output = None

    for i in scipy.arange(21)+1:
        string = str(i)

        if i < 10:
            string = '0'+string
        try:
            temp = ('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.AREA:CHORD_'+string).data()
            tempt = Tree.getNode('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.AREA:CHORD_'+string).dim_of().data()
            
            temp2 = 2*scipy.pi*(.68)*Tree.getNode('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.SIGNAL:CHORD_'+string).data()
            temp2t = Tree.getNode('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.SIGNAL:CHORD_'+string).dim_of().data()

            a = scipy.digitize(tempt,temp2t)
            if output is None:
                output = temp*temp2[a]
            else:
                output = output + temp*temp2

        return (tempt,output)
            
