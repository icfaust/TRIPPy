import TRIPPy.surface
import scipy
import eqtools
#import TRIPPy.plot.mayaplot

#data I found from /home/hutch/work/bolo/ this is horrid and there IS NO DOCUMENTATION
# this is all I've been able to glean JESUS FUCKING CHRIST

def ap(plasma, loc=(1.03220,.198424,.00045),angle=(0.,0.,0.)):
    
    vloc1 = TRIPPy.Vecr(loc)
    vloc2 = TRIPPy.Vecr((loc[0],loc[1],0.))
    vloc2.s = 1.0
    meri = TRIPPy.Vecx((0.,0.,1.))
    area = [.05*.0245,.52*.0245]
    return TRIPPy.surface.Rect(vloc1, plasma, area, vec=[meri,vloc2])

def det(origin, angle=(0.,scipy.pi/2, 0.), loc=(0,0,.0245*.235), radius=50e-6):
    return TRIPPy.surface.Circle(loc, origin, radius=radius, angle=angle) #actually a circular pinhole, but whatever

def dmbolo(num, plasma):
    
    loc = [[1.03220, .198424, .00045],
           [1.0416, 2.452484, .090805],
           [1.02981, 2.89147, -.05969],
           [1.03090, 3.33748, -.002550],
           [1.0305, 5.38974, .00545],
           [1.02872, 6.11321, .00145]]
    
    place = loc[num]

    apin = ap(plasma, place)
    temp = det(apin)
    temp.redefine(plasma)
    return [temp,apin]

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
