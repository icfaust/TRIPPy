import geometry
import surface
import scipy
import scipy.linalg
import _beam

class Ray(geometry.Point):
    r"""Generates a ray vector object
        
    Uses the definition:
        
    .. math::
    
        \vec{x}= \vec{x}_0 + \vec{x}_1
    
    Args:
        pt1: Point or Point-derived object
            Defines the origin of the Ray, based on the coordinate system
            of the origin.  pt1 position is accessible through Ray(0).

        pt2: Point or Vector-derived object
            Direction of the ray can be defined by a vector object (assumed
            to be in the space of the pt1 origin) from pt1, or a point, which 
            generates a vector pointing from pt1 to pt2.
            
    Returns:
        Ray: Ray object.
        
    Examples:
        Accepts all point and point-derived object inputs, though all data 
        is stored as a python object.

        Generate an y direction Ray in cartesian coords using a Vec from (0,0,1)::
            
                cen = geometry.Center(flag=True)
                ydir = geometry.Vecx((0,1,0))
                zpt = geometry.Point((0,0,1),cen)

    """

    def __init__(self, pt1, inp2):
        """
        """

        try:
            self.norm = geometry.pts2Vec(pt1, inp2)
        except AttributeError:
            self.norm = inp2.copy()
            
        super(Ray,self).__init__(pt1)

        self.norm.s = scipy.atleast_1d(self.norm.s)
        self.norm.s = scipy.insert(self.norm.s,0,0.)

    def x(self):
        """returns array of cartesian coordinate in meters

        Returns:
           numpy array of cartesian coordinates in meters

        """
        return (self + self.norm).x()

    def r(self):
        """return cylindrical coordinate values

        Returns:
            numpy array of cylindrical coordinates in meters and radians

        """
        return (self + self.norm).r()

    def t(self, r, z):
        """return toroidal coordinate values for given cylindrical
        coordinates (r,z) in coordinates of the ray origin.

        Args:
            r: scipy-array of floats or float in meters. r is
            specified in meters.

            z: scipy-array of floats or float in meters. z is
            specified in meters.

        Returns:
            numpy array of cylindrical coordinates in [meters,radians,radians]
            where it is radius in meters, toroidal angle and then poloidal angle.
        """
        return (self + self.norm).t(r, z)

    def rmin(self):
        """rmin returns the s value along the norm vector which minimizes
        the r0() value (the closest position to the origin norm axis)

        Returns:
            numpy array of s values in meters
        """
        return -1*self.s*(self.unit[0]*self.norm.unit[0] - 
                          self.unit[1]*self.norm.unit[1]
                          )/(self.norm.unit[0]**2+self.norm.unit[1]**2)

    def smin(self, point=None):
        """Calculates and returns the s value along the norm vector
        which minimizes the distance from the ray to a point 
        (default is the origin which the ray is defined).

        Kwargs:
            point: Point or Point-derived object, otherwise defaults to ray 
            origin
        
        Returns:
            numpy array of s values in meters
        """

        if point is None:
            return -self.s*(self.norm.unit[0]*self.unit[0]+
                            self.norm.unit[1]*self.unit[1]+
                            self.norm.unit[2]*self.unit[2])
        # test in same coordinate system
        if not self._origin is point._origin :
            raise ValueError('not in same coordinate system, use redefine')
        else:
            return -1*( self.norm.unit[0]*(self.unit[0]*self.s - point.x0() )+
                        self.norm.unit[1]*(self.unit[1]*self.s - point.x1() )+
                        self.norm.unit[2]*(self.unit[2]*self.s - point.x2() ))

    def tmin(self, r, z, trace=False):
        """Calculates and returns the s value along the norm vector
        which minimizes the distance from the ray to a circle defined by
        input (r,z). 

        Args:
            r: value, iterable or scipy.array, radius in meters 
            
            z: value, iterable or scipy.array, z value in meters
        
        Kwargs:
            trace: bool if set true, the ray is assumed to be traced within a
            tokamak.  A further evaluation reduces the value to one within
            the bounds of the vacuum vessel/limiter.

        Returns:
            numpy array of s values in meters
        """
        r = scipy.atleast_1d(r)
        z = scipy.atleast_1d(z)
        params = _beam.lineCirc(self(0).x(),
                                self.norm.unit,
                                r,
                                z)

        sout = scipy.zeros(r.shape)

        for i in xrange(len(params)):
            temp = scipy.roots(params[i])

            # only positive real solutions are taken
            temp = temp[scipy.imag(temp) == 0]
            temp = scipy.real(temp[temp > 0])

            test = self(temp).r()
            if not trace:
                # must decide between local and global minima
                sout[i] = temp[((test[0]-r[i])**2 
                                + (test[2] - z[i])**2).argmin()]
            else:
                #need to implement this such that it searches only in area of interest
                sout[i] = temp[scipy.logical_and(temp > self.norm.s[-2],
                                                 temp < self.norm.s[-1])].min()

        return sout
   
    def __getitem__(self,idx):
        return (self + self.norm)[idx]

    def __call__(self,inp):
        """ call is used to minimize the changes to the norm vector.
        it returns a vector"""
        #temporarily store the norm.s
        temp = self.norm.s

        self.norm.s = inp
        out = self + self.norm
        self.norm.s = temp

        return out

    def redefine(self, neworigin):
        """redefine Ray object or Ray-derived object
        into new coordinate system

        Args:
            neworigin: Origin or Origin-derived object
        """

        lca = self._lca(neworigin)
        self._rotate(lca, neworigin)
        super(Ray,self)._translate(lca, neworigin)

    def _rotate(self, lca, neworigin):
        """ rotates the fundamental vectors of the space"""
        org = lca[0]
        orgnew = lca[1]

        temp1 = self.norm

        for idx in range(len(org)-1,-1,-1):
            # change the _rot coordinates to accurately reflect all of the necessary variation.
            temp1 = org[idx].rot(temp1)

        for idx in range(len(orgnew)):
            # the arot allows for translating into the current coordinate system
            temp1 = orgnew[idx].arot(temp1)

        self.norm = temp1

# generate necessary beams for proper inversion (including etendue, etc)
class Beam(geometry.Origin):
    r"""Generates a Beam vector object assuming macroscopic surfaces. A
    beam is a ray which has finite, invariant etendue. The etendue is 
    derived assuming that the angular divergence /cross-sectional area
    can be parameterized by two surfaces.
        
    Uses the definition:
        
    .. math::
    
        \vec{x}= \vec{x}_0 + \vec{x}_1
    
    Args:
        surf1: Surface or Surface-derived object
            Defines the origin surface, based on the coordinate system
            of the surface.  Center position is accessible through Beam(0).
            Generated beam contains same origin as from surf1.

        surf2: Surface or Surface-derived object
            Defines the aperature surface, based on the coordinate system
            of the surface. Both surfaces must be in the same coordinate
            system.
            
    Returns:
        Beam: Beam object.
        
    Examples:
        Accepts all surface or surface-derived object inputs, though all data 
        is stored as a python object.

        Generate an y direction Ray in cartesian coords using a Vec from (0,0,1)::
            
                cen = geometry.Center(flag=True)
                ydir = geometry.Vecx((0,1,0))
                zpt = geometry.Point((0,0,1),cen)

    """

    def __init__(self, surf1, surf2):
        """
        """

        normal = geometry.pts2Vec(surf1, surf2)
        #orthogonal coordinates based off of connecting normal

        snew = surf1.sagi - normal*((surf1.sagi * normal)*(surf1.sagi.s/normal.s))
        mnew = surf1.meri - normal*((surf1.meri * normal)*(surf1.meri.s/normal.s))
        super(Beam, self).__init__(surf1, surf1._origin, vec=[mnew,normal])
        #calculate area at diode.
        self.sagi.s = snew.s
        a1 = surf1.area(snew.s,mnew.s)

        #calculate area at aperature
        a2 = surf2.area((((self.sagi*surf2.sagi)/self.sagi.s)**2
                         + ((self.meri*surf2.sagi)/self.meri.s)**2)**.5,
                        (((self.sagi*surf2.meri)/self.sagi.s)**2 
                         + ((self.meri*surf2.meri)/self.meri.s)**2)**.5)

        #generate etendue
        self.etendue = a1*a2/(normal.s ** 2)

        # give inital beam, which is two points      
        self.norm.s = scipy.atleast_1d(self.norm.s)
        self.norm.s = scipy.insert(self.norm.s,0,0.)

    def x(self):
        """returns array of cartesian coordinate in meters

        Returns:
           numpy array of cartesian coordinates in meters

        """
        return (self + self.norm).x()

    def r(self):
        """return cylindrical coordinate values
        
        Returns:
            numpy array of cylindrical coordinates in meters and radians

        """
        return (self + self.norm).r()

    def t(self, r, z):
        """return toroidal coordinate values for given cylindrical
        coordinates (r,z) in coordinates of the beam origin.

        Args:
            r: scipy-array of floats or float in meters. r is
            specified in meters.

            z: scipy-array of floats or float in meters. z is
            specified in meters.

        Returns:
            numpy array of cylindrical coordinates in [meters,radians,radians]
            where it is radius in meters, toroidal angle and then poloidal angle.
        """
        return (self + self.norm).t(r, z)
    
    def c(self):
        """Conversion of vector to opposite coordinate system

        Returns:
            copy of vector object with opposite coordinate system
            (set with .flag parameter)

        """
        return (self + self.norm).c()

    def rmin(self):
        """Calculates and returns the s value along the norm vector
        which minimizes the r0() value (the closest position to the
        origin norm axis)

        Returns:
            numpy array of s values in meters
        """
        return -1*self.s*(self.unit[0]*self.norm.unit[0] - 
                          self.unit[1]*self.norm.unit[1]
                          )/(self.norm.unit[0]**2+self.norm.unit[1]**2) 

    def smin(self, point=None):
        """Calculates and returns the s value along the norm vector
        which minimizes the distance from the ray to a point 
        (default is the origin which the ray is defined).

        Kwargs:
            point: Point or Point-derived object, otherwise defaults to ray 
            origin
        
        Returns:
            numpy array of s values in meters
        """

        if point is None:
            return -self.s*( self.norm.unit[0]*self.unit[0] +
                             self.norm.unit[1]*self.unit[1] +
                             self.norm.unit[2]*self.unit[2] )
        # test in same coordinate system
        if not self._origin is point._origin :
            raise ValueError('not in same coordinate system, use redefine')
        else:
            return -1*(self.norm.unit[0]*(self.unit[0]*self.s - point.x0())+
                       self.norm.unit[1]*(self.unit[1]*self.s - point.x1())+
                       self.norm.unit[2]*(self.unit[2]*self.s - point.x2()))

    def tmin(self, r, z, trace=False):
        """Calculates and returns the s value along the norm vector
        which minimizes the distance from the ray to a circle defined by
        input (r,z). 

        Args:
            r: value, iterable or scipy.array, radius in meters 
            
            z: value, iterable or scipy.array, z value in meters
        
        Kwargs:
            trace: bool if set true, the ray is assumed to be traced within a
            tokamak.  A further evaluation reduces the value to one within
            the bounds of the vacuum vessel/limiter.

        Returns:
            numpy array of s values in meters
        """
        r = scipy.atleast_1d(r)
        z = scipy.atleast_1d(z)
        params = _beam.lineCirc(self(0).x(),
                                self.norm.unit,
                                r,
                                z)

        sout = scipy.zeros(r.shape)

        for i in xrange(len(params)):
            temp = scipy.roots(params[i])

            # only positive real solutions are taken
            temp = temp[scipy.imag(temp) == 0]
            temp = scipy.real(temp[temp > 0])

            test = self(temp).r()
            
            # must decide between local and global minima
            sout[i] = temp[((test[0]-r[i])**2 
                            + (test[2] - z[i])**2).argmin()]
            
            if trace and ((sout[i] > self.norm.s[-1]) or (sout[i] < self.norm.s[-2])):
                #need to implement this such that it searches only in area of interest
                sout[i] = None

        return sout

    def __getitem__(self,idx):
        return (self + self.norm)[idx]

    def __call__(self,inp):
        """ call is used to minimize the changes to the norm vector.
        it returns a vector"""
        #temporarily store the norm.s
        temp = self.norm.s

        self.norm.s = inp
        out = self + self.norm
        self.norm.s = temp

        return out

class subBeam(Beam):
    r"""Generate an array of Beam objects from two surface objects
    
    Args:
        surf1: Surface object
            Beam origin surfaces, based on the coordinate system
            of the surfaces.  Center position is accessible through Beam(0),
            Beam.x()[...,0] or Beam.r()[...,0] (last two options create
            numpy arrays, the first generats a geometry.Vec object).

        surf2: Surface object
            Direction of the ray can be defined by a vector object (assumed
            to be in the space of the pt1 origin) from pt1, or a point, which 
            generates a vector pointing from pt1 to pt2.
            
    Kwargs:
        split1: two-element tuple
            Describes how many segments to split surf1 in [sagi,meri]

        split2: two-element tuple
            Describes how many segments to split surf1 in [sagi,meri]

    Returns:
        output: multiBeam object.
        
    Examples:
        Accepts all surface or surface-derived object inputs, though all data 
        is stored as a python object.

        Generate an y direction Ray in cartesian coords using a Vec from (0,0,1)::
            
                cen = geometry.Center(flag=True)
                ydir = geometry.Vecx((0,1,0))
                zpt = geometry.Point((0,0,1),cen)

    """

    def __init__(self, surf1, surf2, split1 = None, split2 = None):
        """
        """
        
        #generating grid off of split1, split2 kwargs
        if split1 is None:
            split1 = [1,1]
        if split2 is None:
            split2 = [1,1]

        ins1 = float((split1[0]-1))/split1[0]
        inm1 = float((split1[1]-1))/split1[1]
        ins2 = float((split2[0]-1))/split2[0]
        inm2 = float((split2[1]-1))/split2[1] 

        grid = scipy.meshgrid(surf1.sagi.s*scipy.linspace(-ins1,
                                                           ins1,
                                                           split1[0]),
                              surf1.meri.s*scipy.linspace(-inm1,
                                                           inm1,
                                                           split1[1]),
                              surf2.sagi.s*scipy.linspace(-ins2,
                                                           ins2,
                                                           split2[0]),
                              surf2.meri.s*scipy.linspace(-inm2,
                                                           inm2,
                                                           split2[1]))
        
        #check to see if they are in same coordinate system
        print((surf1.x0() + surf1.sagi(grid[0]).x0()).shape)
        surf1 + surf1.sagi(grid[0])
        surf1cents = surf1 + surf1.sagi(grid[0]) + surf1.meri(grid[1]) #vectors
        surf2cents = surf2 + surf2.sagi(grid[2]) + surf2.meri(grid[3])

        del grid
        self.norm = surf2cents - surf1cents
        self.s = surf1cents.s
        self.unit = surf1cents.unit
        #orthogonal coordinates based off of connecting normal

        self.sagi = surf1.sagi - self.norm*((surf1.sagi * self.norm)*(surf1.sagi.s/self.norm.s))
        self.meri = surf1.meri - self.norm*((surf1.meri * self.norm)*(surf1.meri.s/self.norm.s))
        
        #reduce calcuations in calling super to inherited classes
        self._origin = surf1._origin
        self._depth = surf1._depth
        self.flag = surf1.flag

        #calculate area at diode.
        self.sagi.s /= split1[0]
        self.meri.s /= split1[1]
        a1 = surf1.area(self.sagi.s,self.meri.s)

        #calculate area at aperature
        a2 = surf2.area((((self.sagi*surf2.sagi)/self.sagi.s)**2
                         + ((self.meri*surf2.sagi)/self.meri.s)**2)**.5,
                        (((self.sagi*surf2.meri)/self.sagi.s)**2 
                         + ((self.meri*surf2.meri)/self.meri.s)**2)**.5)/(split2[0]*split2[1])

        #generate etendue
        self.etendue = a1*a2/(self.norm.s ** 2)
        self.shape = self.s.shape[:]
        self.main = Beam(surf1,surf2) #minimal memory waste, but infinitely useful for minimizing calculations
        #self.flatten()
        # flatten data, but store original shapes.

    def flatten(self):

        self.etendue = self.etendue.ravel()
        self.s = self.s.ravel()
        self.unit = self.unit.reshape((3,self.unit.size/3))
        self.norm.s = self.norm.s.ravel()
        self.norm.unit = self.unit.reshape((3,self.norm.unit.size/3))

    def reshape(self):
        raise NotImplementedError(' not yet')

    def __getitem__(self,idx):
        return geometry.Vec(self.unit[idx],s=self.s[idx]) + self.norm[idx]

    def __call__(self,inp):
        """ call is used to minimize the changes to the norm vector.
        it returns a vector"""
        #temporarily store the norm.s
        temp = self.norm.s

        self.norm.s = inp
        out = self + self.norm
        self.norm.s = temp

        return out

def multiBeam(surf1, surf2, split=None):
    r"""Generate a tuple of Beam objects from tuples of surface objects
    
    Args:
        surf1: tuple of Surfaces or a Surface object
            Beam origin surfaces, based on the coordinate system
            of the surfaces.  Center position is accessible through Beam(0),
            Beam.x()[...,0] or Beam.r()[...,0] (last two options create
            numpy arrays, the first generats a geometry.Vec object).

        surf2: tuple of Surfaces or a Surface object
            Direction of the ray can be defined by a vector object (assumed
            to be in the space of the pt1 origin) from pt1, or a point, which 
            generates a vector pointing from pt1 to pt2.
            
    Returns:
        output: tuple of beam objects.
        
    Examples:
        Accepts all surface or surface-derived object inputs, though all data 
        is stored as a python object.

        Generate an y direction Ray in cartesian coords using a Vec from (0,0,1)::
            
                cen = geometry.Center(flag=True)
                ydir = geometry.Vecx((0,1,0))
                zpt = geometry.Point((0,0,1),cen)

    """
    if not split is None:
        surf1 = surf1.split(split[0],split[1])
        surf2 = surf2.split(split[0],split[1])
    
    output = []
    try:
        output += [Beam(surf1,surf2)]
                
    except AttributeError:
        try:
            for i in surf1:
                try:
                    output += [Beam(i,surf2)]
                except AttributeError:
                    output += multiBeam(i,surf2)
                    
        except TypeError:
            for i in surf2:
                try:
                    output += [Beam(surf1,i)]
                except AttributeError:
                    output += multiBeam(surf1,i)
                
    return output


def volWeightBeam(beam, rgrid, zgrid, trace=True, ds=2e-3, toroidal=None, **kwargs):
    r"""Generate a tuple of Beam objects from tuples of surface objects
    
    Args:
        beam: tuple of Surfaces or a Surface object
            Beam origin surfaces, based on the coordinate system
            of the surfaces.  Center position is accessible through Beam(0),
            Beam.x()[...,0] or Beam.r()[...,0] (last two options create
            numpy arrays, the first generats a geometry.Vec object).

        rgrid: tuple of Surfaces or a Surface object
            Direction of the ray can be defined by a vector object (assumed
            to be in the space of the pt1 origin) from pt1, or a point, which 
            generates a vector pointing from pt1 to pt2.
            
        zgrid: tuple of Surfaces or a Surface object
            Direction of the ray can be defined by a vector object (assumed
            to be in the space of the pt1 origin) from pt1, or a point, which 
            generates a vector pointing from pt1 to pt2.

    Returns:
        output: tuple of beam objects.
        
    Examples:
        Accepts all surface or surface-derived object inputs, though all data 
        is stored as a python object.

        Generate an y direction Ray in cartesian coords using a Vec from (0,0,1)::
            
                cen = geometry.Center(flag=True)
                ydir = geometry.Vecx((0,1,0))
                zpt = geometry.Point((0,0,1),cen)

    """
    out = scipy.zeros((len(rgrid)-1,len(zgrid)-1))
    try:
        if toroidal is None:
            if trace:
                temp = beam(scipy.mgrid[beam.norm.s[-2]:beam.norm.s[-1]:ds]).r()
            else:
                temp = beam(scipy.mgrid[beam.norm.s[0]:beam.norm.s[-1]:ds]).r()
                
            out += scipy.histogram2d(temp[0],
                                     temp[2], 
                                     bins = [rgrid, zgrid],
                                     weights=scipy.ones(temp[0].shape)*beam.etendue*ds)[0]
        else:
            if trace:
                temp = beam(scipy.mgrid[beam.norm.s[-2]:beam.norm.s[-1]:ds]).t(toroidal[0],toroidal[1])
            else:
                temp = beam(scipy.mgrid[beam.norm.s[0]:beam.norm.s[-1]:ds]).t(toroidal[0],toroidal[1])
            
            out += scipy.histogram2d(temp[2],
                                     temp[0],
                                     bins = [rgrid, zgrid],
                                     weights=scipy.ones(temp[0].shape)*beam.etendue*ds)[0]
    except AttributeError:
        for i in beam:
            try:
                out += volWeightBeam(i, rgrid, zgrid, trace=trace, ds=ds, toroidal=toroidal, **kwargs)
            except TypeError:
                pass

    return out


def _genBFEdgeZero(plasma, zeros, rcent, zcent):
    """ this will absolutely need to be rewritten"""

    theta = scipy.linspace(-scipy.pi,scipy.pi,zeros)
    cent = geometry.Point(geometry.Vecr([rcent,0,zcent]),plasma)
    zerobeam = []
    outline = []
    for i in xrange(len(plasma.norm.s)-1):
        outline += [geometry.Vecx([plasma.sagi.s[i],
                                   0,
                                   plasma.norm.s[i]])-cent]
        
    for i in xrange(zeros):
        temp2 = geometry.Vecr([scipy.cos(theta[i]),
                               0,
                               scipy.sin(theta[i])])
        s = 0
        for j in outline:
            temp4 = j*temp2
            if temp4 > s:
                s = temp4

        temp2.s = s
        zerobeam += [Ray(geometry.Point(cent+temp2,
                                        plasma),
                         geometry.Vecr([scipy.sin(theta[i]),
                                        0,
                                        -scipy.cos(theta[i])]))]

    return zerobeam


def pos2Ray(pos, tokamak, angle=None, eps=1e-6):
    r"""Take in GENIE pos vectors and convert it into TRIPPy rays
    
    Args:
        pos: 4 element tuple or 4x scipy-array
            Each pos is assembled into points of (R1,Z1,RT,phi)

        tokamak: 
            Tokamak object in which the pos vectors are defined.
            
    Returns:
        Ray: Ray object or typle of ray objects.
        
    """

    r1 = scipy.array(pos[0])
    z1 = scipy.array(pos[1])
    rt = scipy.array(pos[2])
    phi = scipy.array(pos[3])

    zt = z1 - scipy.tan(phi)*scipy.sqrt(r1**2 - rt**2)
    angle2  = scipy.arccos(rt/r1)

    if angle is None:
        angle = scipy.zeros(r1.shape)

    pt1 = geometry.Point((r1,angle,z1),tokamak)
    pt2 = geometry.Point((rt,angle+angle2,zt),tokamak)

    output = Ray(pt1,pt2)
    output.norm.s[-1] = eps
    tokamak.trace(output)
    return output
    
