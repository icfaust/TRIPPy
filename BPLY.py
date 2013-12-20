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
    # Vec in this case determines the meridonial and normal rays in the 
    Vec = [geometry.Vecx((1.,0.,0.)),geometry.Vecx((0.,0.,1.))]
    aperature = surface.Rect((0.,0.,0.),pos,area,Vec=Vec,flag=False)
    
    Vec = [geometry.Vecx((0.,1.,0.)),geometry.Vecx((0.,0.,1.))]
    offset = geometry.Origin((0.,-2.5e-3,-7.389e-2),aperature,Vec=Vec)
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
        output[i] = beamin.Beam(temp[i],temp[-1])
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
        temp = lim1
        mask = scipy.logical_and(scipy.isfinite(psi[i]),psi[i] < lim2+.02)

        try:
            minpos = scipy.argmin(psi[i][mask])
            test = psi[i][mask][minpos]
        except ValueError:
            test = lim2+.03
            
        #plt.plot(beam.x()[0][mask],psi[i][mask])
        #plt.show()
        sizer = psi[i][mask].size
        if not test > lim2:

        #plt.plot(beam.x()[0][mask][0:minpos],psi[i][mask][0:minpos],beam.x()[0][mask][minpos:],psi[i][mask][minpos:])
        #plt.show()
        #limout = scipy.insert(lim,(2,2),(beam.norm.s[mask][minpos],beam.norm.s[mask][minpos]))  # add minimum flux s for bound testing
            if lim1 < test:
                temp = test

            try:
                temp1 = scipy.clip(scipy.digitize((lim1,lim2),psi[i][mask][minpos::-1]),0,minpos)
                outspline[i] = beam.norm.s[mask][minpos::-1][temp1]
            
            except ValueError:
                tempmask = (psi[i][mask] < lim2)[0]
                outspline[i] = scipy.array([beam.norm.s[mask][minpos],beam.norm.s[mask][tempmask]])

            try:
                temp2 = scipy.clip(scipy.digitize((lim1,lim2),psi[i][mask][minpos:]),0,sizer-minpos-1)
                inspline[i] = beam.norm.s[mask][minpos:][temp2]
                
            except ValueError:
                inspline[i] = scipy.array([beam.norm.s[mask][minpos],beam.norm.s[mask][-1]])

        else:
            outspline[i] = scipy.array([[],[]])
            inspline[i] = scipy.array([[],[]])

    return (outspline,inspline)

def calcArea(points):
    val = 0 
    for i in scipy.arange(len(points))-1:
        val += points[i][0]*points[i+1][1] - points[i][1]*points[i+1][0]
    return abs(val/2)

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

    
    output = t.size * [[None]]
    #outermid,innermid = getBeamFluxSpline(beam,plasma,t,lim1,lim2)
    #outertop,innertop = getBeamFluxSpline(ray1,plasma,t,lim1,lim2)
    #outerbot,innerbot = getBeamFluxSpline(ray2,plasma,t,lim1,lim2)
    mid = getBeamFluxSpline(beam,plasma,t,lim1,lim2)
    top = getBeamFluxSpline(ray1,plasma,t,lim1,lim2)
    bot = getBeamFluxSpline(ray2,plasma,t,lim1,lim2)

    #plt.plot(mid[0])
    #plt.plot(mid[1])
    #plt.show()

    #print(mid)
    #print(top)

    #condition inputs for area calculations
    for i in xrange(t.size):
        output[i] = []
        #segment = 3*[0]
        #beam and ray masking values are already written to their norm.s values
#        segment[0] = scipy.array([outertop[i][1],outertop[i][0],innertop[i][0],innertop[i][1]])
#        segment[1] = scipy.array([outermid[i][1],outermid[i][0],innermid[i][0],innermid[i][1]])
#        segment[2] = scipy.array([outerbot[i][1],outerbot[i][0],innerbot[i][0],innerbot[i][1]])
        for j in xrange(2):
            temp = []
            if scipy.any(mid[j][i]):
                for k,ray in (mid[j][i],beam),(top[j][i],ray1),(bot[j][i],ray2):
                    if scipy.any(k):
                        ray.norm.s = k
 
                        if fillorder:
                            #print(ray.x())
                            temp += [ray.x()]
                        else:
                            #print( ray.norm.s)
                            temp += [ray.split(plasma,obj=geometry.Point)]

                output[i] += [temp]



        # compare and mask/replace
        #scipy.copyto(segment[0],ray1.norm.s,where=~scipy.isfinite(segment[0]))
        #scipy.copyto(segment[1],beam.norm.s,where=~scipy.isfinite(segment[1]))
       # scipy.copyto(segment[2],ray2.norm.s,where=~scipy.isfinite(segment[2]))
 
        #turn into points
        #ray1.norm.s = segment[0]
        #beam.norm.s = segment[1]
        #ray2.norm.s = segment[2]

        #temp1 = ray1.split(plasma,obj=geometry.Point)
        #temp2 = beam.split(plasma,obj=geometry.Point)
        #temp3 = ray2.split(plasma,obj=geometry.Point)

        #if fillorder:
        #    output[i] = []

            
         #   for j in ((0,1,1,1,0,0),(2,3,3,3,2,2)):
         #       output[i] += [scipy.array([temp1[j[0]].x(),
         #                                  temp1[j[1]].x(),
         #                                  temp2[j[2]].x(),
         #                                  temp3[j[3]].x(),
         #                                  temp3[j[4]].x(),
         #                                  temp2[j[5]].x()]).T]
        #else:
        #    output[i] = [temp1,temp2,temp3]
    return output
    





def effectiveHeight(surf1, surf2, plasma, t, lim1=.88, lim2=.92):
    """ calculate the effective height of a view through the scrape-off-layer"""

    segments = viewPoints(surf1,surf2,plasma,t,lim1=lim1,lim2=lim2,fillorder=False)
    output = scipy.zeros((len(segments),))
    for i in xrange(len(segments)):
        area = 0
        if not scipy.all(segments[i] == []):
            
            inlen = geometry.pts2Vec(segments[i][0][0][0],segments[i][0][0][1])
            outlen = geometry.pts2Vec(segments[i][1][0][0],segments[i][1][0][1])

 
            for j in xrange(len(segments[i])): # loop over in vs out

                temp = []
                for k in segments[i][j]: # loop over number of 'rays'
                    for l in k:
                        temp += [l.x()[[0,2]]]
                temp = scipy.array(temp)
        #delaunay
                for k in xrange((len(temp)-2)/2):
                    y = temp[[0,1,2*(k+1),2*(k+1)+1]]
                    if not(scipy.all([y[0] == y[1],y[2] == y[3]])):
                        tri = scipy.spatial.Delaunay(y) #divide into a lower and upper section
               
                #plt.triplot(y[:,0],y[:,1],tri.vertices.copy())
                        for l in tri.points[tri.vertices]:
                            area += calcArea(l)


                output[i] = area/(inlen.s + outlen.s)               
    #plt.show()
    return output


def writeToTree(inp):
    shot,chan,surf1,surf2 = inp
    timeout = time.time()
    outstr = str(chan)
    if chan < 10:
        outstr = '0' + outstr
    try:    
        Tree = MDSplus.Tree('spectroscopy',shot)
        b = plasma.Tokamak(eqtools.CModEFITTree(shot,tspline=True)) # I HATE THIS
        output = effectiveHeight(surf1,surf2,b,scipy.mgrid[b.eq.getTimeBase()[0]:b.eq.getTimeBase()[-1]:1e-3])
        dummy =  MDSplus.Data.compile('build_signal($1,*,$2)',output,scipy.mgrid[b.eq.getTimeBase()[0]:b.eq.getTimeBase()[-1]:1e-3])
        Tree.getNode('.BOLOMETER.RESULTS.DIODE.BPLY.AREA:CHORD_'+outstr).putData(dummy)
        print('channel '+str(chan)+':\t '+str(time.time()-timeout))
    except:
        print('ERROR IN THIS CHANNEL: '+str(chan))

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
        out = pool.map(writeToTree,[[shot,22-i,n[i],n[-1]] for i in scipy.arange(22)])
        pool.close()
        print(shot)

    except MDSplus._treeshr.TreeException:
        print(shot)
        gc.collect()

    except MDSplus._tdishr.TdiException:
        print(shot)
        gc.collect()

    except ValueError:
        print(shot)
        gc.collect()


def globalpowerCalc(shot):

    Tree = MDSplus.Tree('spectroscopy',shot)
    output = None
    temp2 = 2*scipy.pi*(.68)*Tree.getNode('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.BRIGHT').data() #does factor have the 4pi?
    temp2t = Tree.getNode('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.BRIGHT').dim_of().data()
    for i in scipy.arange(20)+2:
        string = str(i)

        if i < 10:
            string = '0'+string
       # try:
        temp = Tree.getNode('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.AREA:CHORD_'+string).data()
        
        tempt = Tree.getNode('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.BPLY.AREA:CHORD_'+string).dim_of().data()
        

        a = scipy.digitize(tempt,temp2t)
        if output is None:
            output = temp*temp2[i][a]
        else:
            output = output + temp*temp2[i][a]
        #except ValueError:
        #    print('no')

    return (tempt,output)

def writeGlobal(shot):
    try:
        out = globalpowerCalc(shot)
        Tree = MDSplus.Tree('spectroscopy',shot)
        dummy =  MDSplus.Data.compile('build_signal($1,*,$2)',out[1],out[0])
        Tree.getNode('.BOLOMETER.RESULTS.DIODE.LYMAN_AVE').putData(dummy)
    except:
        print('error in shot: '+str(shot))

            
def trend(tree,signal,shot):
    temp = MDS.Tree(tree,signal)
    xt = temp.getNode(signal).dim_of().data() 
    x = temp.getNode(signal).dim_of().data() 
    yt,y = globalpowerCalc(shot)
    a = scipy.digitize(xt,yt)
    return x,y[a]
