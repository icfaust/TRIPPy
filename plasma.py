import geometry,scipy,eqtools

class Tokamak(geometry.Center):

    def __init__(self,equilib,flag=True):
        self.eq = equilib
        super(Tokamak,self).__init__(flag=flag)

    def genWireFrame(self,pts=250):
        theta = scipy.linspace(0,2*scipy.pi,pts)
        theta = scipy.append(theta,0)
        x = scipy.cos(theta)
        y = scipy.sin(theta)
        z = scipy.ones(theta.shape)
        return x,y,z

    def inVessel(self,ptin):
        Rlim,Zlim = self.eq.getMachineCrossSection()
        return bool(eqtools.inPolygon(Rlim,Zlim,ptin[0],ptin[2]))
            
    def getVessel(self,idx,pts=250):
        Rlim,Zlim = self.eq.getMachineCrossSection()
        if idx < Rlim.size:
            theta = scipy.linspace(0,2*scipy.pi,pts)
            return geometry.Point((Rlim[idx]*scipy.ones(theta.size),
                                   theta,
                                   Zlim[idx]*scipy.ones(theta.size)),self)

    # to be added, inversion functions etc! I'm glad I spent the day on this.

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

    def _intersect2d(self,ray,s1,s2):
        print('temp')
