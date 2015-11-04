import TRIPPy.surface
import TRIPPy
import AXUV
import scipy
import eqtools
import TRIPPy.plot.mayaplot

#data I found from /home/hutch/work/bolo/ this is horrid and there IS NO DOCUMENTATION
# this is all I've been able to glean JESUS FUCKING CHRIST

def ap(plasma, loc=(.9344,.0,-.0311),angle=.863204): #.042962
    
    vloc1 = TRIPPy.Vecr(loc)
    vloc2 = TRIPPy.Vecr((loc[0],scipy.pi+loc[1],0))
    vloc2.s = 1.0
    vloc2.spin(angle)
    meri = TRIPPy.cross(vloc2,TRIPPy.Vecx((0.,0.,1.)))
    area = [.297e-2,.074e-2]
    return TRIPPy.surface.Rect(vloc1, plasma, area, vec=[meri,vloc2])

def det(origin, x0=(1.64e-2,-.717e-2,-3.55e-2), x1=(2.64e-2,-.717e-2,-3.55e-2), x2=(1.64e-2,.284e-2,-3.64e-2)):

    x0vec = TRIPPy.Vecx(x0)
    x1vec = TRIPPy.Vecx(x1)
    x2vec = TRIPPy.Vecx(x2)
    meri = x2vec-x0vec
    sagi = x1vec-x0vec
    meri.s = 1.0
    sagi.s = 1.0
    norm = TRIPPy.cross(meri,sagi)
    org = TRIPPy.Origin(x0vec,origin,vec=[meri,norm])

    return AXUV.AXUV20(org)#actually a circular pinhole, but whatever

def axdiv(plasma):

    apin = ap(plasma)
    temp = det(apin)
    for i in temp:
        i.redefine(plasma)
    temp.append(apin)
    return temp


def det2(origin, x0=(1.64e-2,-.717e-2,-3.55e-2), x1=(2.64e-2,-.717e-2,-3.55e-2), x2=(1.64e-2,.284e-2,-3.64e-2)):

    x0vec = TRIPPy.Vecx(x0)
    x1vec = TRIPPy.Vecx(x1)
    x2vec = TRIPPy.Vecx(x2)
    meri = x2vec-x0vec
    sagi = x1vec-x0vec
    meri.s = 1.0
    sagi.s = 1.0
    norm = TRIPPy.cross(meri,sagi)
    org = TRIPPy.Origin(x0vec,origin,vec=[meri,norm])

    return AXUV.AXUV20(org)#actually a circular pinhole, but whatever

def axdiv2(plasma):

    apin = ap(plasma)
    temp = det2(apin)
    for i in temp:
        i.redefine(plasma)
    temp.append(apin)
    return temp

def volweight(num,numsplit=(3,3), factor=1, fact2=None, eq='/home/ian/python/g1120824019.01400'):

    b =  TRIPPy.Tokamak(eqtools.EqdskReader(gfile=eq))
    
    rgrid = b.eq.getRGrid()
    zgrid = b.eq.getZGrid()
    rgrid = scipy.linspace(rgrid[0],rgrid[-1],len(rgrid)*factor)
    zgrid = scipy.linspace(zgrid[0],zgrid[-1],len(zgrid)*factor)
    
    dmbolo2 = dmbolo(num,b)
    
    surfs = dmbolo2[1].split(numsplit[0],numsplit[1])

    out = scipy.zeros((len(rgrid)-1,len(zgrid)-1))

    print(dmbolo2[0],dmbolo2[1])

    for i in surfs:
        for j in i:
            surf2 = j
            if fact2 is None:
                surf2 = j
            else:
                surf2 = j.split(fact2,fact2)

            beam = TRIPPy.beam.multiBeam(dmbolo2[0],surf2)
            b.trace(beam)
            TRIPPy.plot.mayaplot.plotLine(beam)
            #out += TRIPPy.beam.volWeightBeam(beam,rgrid,zgrid)
    
    return out


temp =[[0.969362,   -0.0147000,     0.387299,     0.404814],
       [0.968710,   -0.0147000,     0.407214,     0.407785],
       [0.968057,   -0.0147000,     0.427165,     0.410602],
       [0.967406,   -0.0147000,     0.447241,     0.413224],
       [0.966755,   -0.0147000,     0.467272,     0.415668],
       [0.966104,   -0.0147000,     0.487278,     0.417908],
       [0.965454,   -0.0147000,     0.507263,     0.419921],
       [0.964804,   -0.0147000,     0.527075,     0.421725],
       [0.964155,   -0.0147000,     0.546720,     0.423297],
       [0.963507,   -0.0147000,     0.566112,     0.424646],
       [0.962858,   -0.0147000,     0.585326,     0.425719],
       [0.962211,   -0.0147000,     0.604192,     0.426555],
       [0.961564,   -0.0147000,     0.622744,     0.427121],
       [0.960917,   -0.0147000,     0.640929,     0.427413],
       [0.960271,   -0.0147000,     0.658644,     0.427462],
       [0.959625,   -0.0147000,     0.675887,     0.427253],
       [0.958980,   -0.0147000,     0.692623,     0.426788],
       [0.958336,   -0.0147000,     0.708846,     0.426052],
       [0.957692,   -0.0147000,     0.724501,     0.425068],
       [0.957048,   -0.0147000,     0.739587,     0.423826]]
