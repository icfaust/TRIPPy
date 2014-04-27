import geometry,scipy,eqtools
import surface
import scipy.linalg

class Tokamak(geometry.Center):

    def __init__(self, equilib, flag=True):
        self.eq = equilib
        super(Tokamak, self).__init__(flag=flag)
        self.meri.s = self.eq.getMachineCrossSection()[0] #store R and Z of limiter structure vacuum vessel
        self.sagi.s = self.meri.s #link together.
        self.norm.s = self.eq.getMachineCrossSection()[1]

    def inVessel(self, ptin):
        return bool(eqtools.inPolygon(self.sagi.s,self.norm.s,ptin[0],ptin[2]))
            
    def getVessel(self, idx, pts=250):
        if idx < self.norm.s.size:
            theta = scipy.linspace(0,2*scipy.pi,pts)
            return surface.Circle((0.,0.,self.norm.s[idx]),
                                  self,
                                  self.meri.s[idx],
                                  vec = [self.meri,self.norm])

    def getMachineCrossSection(self):
        return geometry.Point(self.norm+self.sagi,self)

    def pnt2RhoTheta(self, point,t=0, method = 'psinorm', n=0, poloidal_plane=0):
        """ takes r,theta,z and the plasma, and map it to the toroidal position
        this will be replaced by a toroidal point generating system based of
        a seperate vector class"""
        try:
            if not point.flag:
                point = point.c()
        except AttributeError:
            if not point._origin.flag:
                point = point.c()

        r = self.eq.getMagR()
        z = self.eq.getMagZ()

        rho = self.eq.rz2rho(method,point[0],point[2],t)
        theta = (scipy.arctan2(point[2]-z,point[0]-r) + (n*point[1]) + scipy.pi)% (2*scipy.pi) - scipy.pi
        return rho,theta

    def gridWeight(self, beams, rgrid=None, zgrid=None, spacing=1e-3):
        """ this method finds the weighting of the view on a poloidal plane
        it takes a beam and traces it through the plasma"""
        
        if rgrid is None:
            rgrid = self.eq.getRGrid()
    
        if zgrid is None:
            zgrid = self.eq.getZGrid()

        

        print('not implemented')
