import geometry,scipy,eqtools
import surface
import scipy.linalg
import _beam

class Tokamak(geometry.Center):
    """Tokamak object for handling geometry and equilibria.
    
    Creates a new Point instance which is set to a default 
    coordinate system of cartesian or cylindrical coordinates
    determined from the reference coordinate system. All vector
    mathematics are accomplished in cartesian to simplify 
    computation effort. Cylindrical vectors are converted at
    last step for return.

    Args:
        equilib: An eqtools Equilibrium object which defines the
             flux surfaces (psi) on a regular grid, as well as 
             the wall geometry amongst a host of other important 
             parameters used within the invert.py submodule.
             This requires the installation of the eqtools package
             which handles tokamak equilibria in a standardized 
             pythonic form.  It is the use of different equilibria
             in a standardized form which allows TRIPPy to be 
             abstract, and capable for use on any tokamak.

    Kwargs:
        flag: Origin object.
            Sets the default coordinate nature 
                
    Examples:   
        Accepts all array like (tuples included) inputs, though
        all data is stored in numpy arrays.

        Generate a point 1m in x from point U (xdir)::
            
                xdir = Vec((1.,0.,0.),U) #implicitly in cyl. coords.

        Generate a cartesian point into direction (2,2,2) from Q::
            
                vec2 = Vec((2.,2.,2.), ref=Q)

        Generate a cylindrical radial vector (vec2) into direction (1,0,1)::
            
                cent = Center() #implicitly in cyl. coords.
                vec2 = Vec(vec1, ref=cent)
    """
    
    def __init__(self, equilib, flag=True):
        """
        """
        self.eq = equilib
        super(Tokamak, self).__init__(flag=flag)
        self.meri.s = self.eq.getMachineCrossSection()[0] #store R and Z of limiter structure vacuum vessel
        self.sagi.s = self.meri.s #link together.
        self.norm.s = self.eq.getMachineCrossSection()[1]

    def center(self,time):
        return [self.eq.getMagRSpline()(time),self.eq.getMagZSpline()(time)]

    def inVessel(self, ptin):
        return bool(eqtools.inPolygon(self.sagi.s,self.norm.s,ptin[0],ptin[2]))
            
    def getVessel(self, idx, pts=250):
        if idx < self.norm.s.size:
            theta = scipy.linspace(0,2*scipy.pi,pts)
            return surface.Circle((0.,0.,self.norm.s[idx]),
                                  self,
                                  self.meri.s[idx],
                                  vec = [self.meri.copy(),self.norm.copy()])

    def getMachineCrossSection(self):
        return geometry.Point(self.norm+self.sagi,self)

    def trace(self, ray, limiter=0):
        try:

        # norm vector is modfied following the convention set in geometry.Origin
            if not ray._origin is self:
                ray.redefine(self)
                
            invesselflag = self.inVessel(ray.r()[...,-1])
            
            intersect = _beam.interceptCyl(scipy.atleast_2d(ray.x()[:,-1]), 
                                           scipy.atleast_2d(ray.norm.unit), 
                                           self.meri.s,
                                           self.norm.s) + ray.norm.s[-1]
            if scipy.isfinite(intersect):
                ray.norm.s = scipy.append(ray.norm.s,intersect)

            intersect = _beam.interceptCyl(scipy.atleast_2d(ray.x()[:,-1]),
                                           scipy.atleast_2d(ray.norm.unit),
                                           self.meri.s,
                                           self.norm.s) + ray.norm.s[-1]

            if not invesselflag and scipy.isfinite(intersect):
                ray.norm.s = scipy.append(ray.norm.s,intersect)
        
            # This is used when the actual plasma vessel structure is not as described in the eqdsk
            # example being the Limiter on Alcator C-Mod, which then keys to neglect an intersection,
            # and look for the next as the true wall intersection.
            for i in xrange(limiter):
                intersect = _beam.interceptCyl(scipy.atleast_2d(ray.x()[:,-1]),
                                               scipy.atleast_2d(ray.norm.unit),
                                               self.meri.s,
                                               self.norm.s) + ray.norm.s[-1]
                ray.norm.s[-1] = intersect
        
        except AttributeError:
            for i in ray:
                self.trace(i, limiter=limiter)

        except ValueError:
            # for subBeams, reference main beam
            self.trace(ray.main)
            #self.norm.s = 

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
        theta = (scipy.arctan2(point[2]-z,point[0]-r)
                 + (n*point[1]) + scipy.pi)% (2*scipy.pi) - scipy.pi
        return rho,theta

    def gridWeight(self, beams, rgrid=None, zgrid=None, spacing=1e-3):
        """ this method finds the weighting of the view on a poloidal plane
        it takes a beam and traces it through the plasma"""
        
        if rgrid is None:
            rgrid = self.eq.getRGrid()
    
        if zgrid is None:
            zgrid = self.eq.getZGrid()

        print('not implemented')
