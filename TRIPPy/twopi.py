import TRIPPy.surface
import scipy
import eqtools
#import TRIPPy.plot.mayaplot

#data I found from /home/hutch/work/bolo/ this is horrid and there IS NO DOCUMENTATION
# this is all I've been able to glean JESUS FUCKING CHRIST

def det(origin,ang1=120.,ang2=150.,offset=1e-3,angle=(0.,0.,0.)):
    
    norm = TRIPPy.Vecx((0.,0.,1.))
    meri = TRIPPy.Vecx((0.,1.,0.))
    area = [2*offset*scipy.tan(ang1*scipy.pi/360.),2*offset*scipy.tan(ang2*scipy.pi/360.)]
    return TRIPPy.surface.Rect((0.,0.,-offset), origin, area, angle=angle)

def ap(plasma, angle=(0., scipy.pi/3, 0.), loc=(.922,0,-.261), area=scipy.pi*pow(5e-5,2)):
    return TRIPPy.surface.Rect(loc,plasma,area=[scipy.sqrt(area),scipy.sqrt(area)],angle=angle) #actually a circular pinhole, but whatever

def twopi(plasma,angle=(0.,scipy.pi/3,0.),loc=(.922,0,-.261),area=scipy.pi*pow(5e-5,2)):
    apin = ap(plasma,angle=angle,loc=loc,area=area)
    temp = det(apin)
    temp.redefine(plasma)
    return [temp,apin]


def volweight(numsplit=(3,3), factor=1, fact2=None, eq='/home/ian/python/g1120824019.01400'):

    b =  TRIPPy.Tokamak(eqtools.EqdskReader(gfile=eq))
    
    rgrid = b.eq.getRGrid()
    zgrid = b.eq.getZGrid()
    rgrid = scipy.linspace(rgrid[0],rgrid[-1],len(rgrid)*factor)
    zgrid = scipy.linspace(zgrid[0],zgrid[-1],len(zgrid)*factor)
    
    twopi2 = twopi(b)
    
    surfs = twopi2[0].split(numsplit[0],numsplit[1])

    out = scipy.zeros((len(rgrid)-1,len(zgrid)-1))

    for i in surfs:
        for j in i:
            surf2 = j
            if fact2 is None:
                surf2 = j
            else:
                surf2 = j.split(fact2,fact2)

            beam = TRIPPy.beam.multiBeam(surf2,twopi2[1])
            b.trace(beam)
            #TRIPPy.plot.mayaplot.plotLine(beam)
            out += TRIPPy.beam.volWeightBeam(beam,rgrid,zgrid)

    
    return out
