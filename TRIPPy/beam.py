import geometry
import surface
import scipy
import scipy.linalg
import _beam
#import matplotlib.pyplot as plt

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
            self.norm = inp2
            
        super(Ray,self).__init__(pt1)
        self.norm.s = scipy.concatenate(([0.],self.norm.s))

    def x0(self):
        """returns cartesian coordinate along first dimension

        Returns:
           numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[0] + self.norm.s*self.norm.unit[0]

    def x1(self):
        """returns cartesian coordinate along second dimension

        Returns:
            numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[1] + self.norm.s*self.norm.unit[1]

    def x2(self):
        """returns cartesian coordinate along third dimension

        Returns:
            numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[2] + self.norm.s*self.norm.unit[2]

    def r(self):
       """return cylindrical coordinate values

        Returns:
            numpy array of cylindrical coordinates in meters and radians

        """
        return (self + self.norm).r()

    def rmin(self):
        """rmin returns the s value along the norm vector which minimizes
        the r0() value (the closest position to the origin norm axis

        Returns:
            numpy array of s values in meters
        """
        return -1*self.s*(self.unit[0]*self.norm.unit[0] - 
                          self.unit[1]*self.norm.unit[1]
                          )/(self.norm.unit[0]**2+self.norm.unit[1]**2)

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

    def tangency(self, point, sigma=False):
        """ returns a vector which points from the point to the closest
        approach of the ray"""


        # test in same coordinate system
        if not ((self._origin is point._origin) or(self._origin is point)) :
            raise ValueError('not in same coordinate system, use redefine')

        # define vector r1, from point to ray origin
        temp = self(0)
        r1 = temp - point

        # define tangency sigma (or length along rtan to tangency point
        sigma = -(r1*self.norm)/(temp.s**2)
        if sigma:
            return sigma
        else:
            return r1 + geometry.Vecx(sigma*self.norm.unit)

    def point(self,err=[]):
        return Point((self+self.norm), self.ref, err=err)


# generate necessary beams for proper inversion (including etendue, etc)
class Beam(geometry.Origin):
    r"""Generates a Beam vector object assuming macroscopic surfaces
        
    Uses the definition:
        
    .. math::
    
        \vec{x}= \vec{x}_0 + \vec{x}_1
    
    Args:
        surf1: Surface or Surface-derived object
            Defines the origin surface, based on the coordinate system
            of the surface.  Center position is accessible through Beam(0).

        surf2: Surface or Surface-derived object
            Direction of the ray can be defined by a vector object (assumed
            to be in the space of the pt1 origin) from pt1, or a point, which 
            generates a vector pointing from pt1 to pt2.
            
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
        a2 = surf2.area((((self.sagi*surf2.sagi)/self.sagi.s)**2 + ((self.meri*surf2.sagi)/self.meri.s)**2)**.5,
                        (((self.sagi*surf2.meri)/self.sagi.s)**2 + ((self.meri*surf2.meri)/self.meri.s)**2)**.5)

        #generate etendue
        self.etendue = a1*a2/(normal.s ** 2)

        # give inital beam, which is two points      
        self.norm.s = scipy.insert(self.norm.s,0,0.)


    def tangent(self,point=None):
        """ returns the point of closest approach of the beam as
        defined by its position and normal vector """
        if point is None:
            point = self._origin


    def intercept(self,surface):
        if self._origin is surface._origin:
            try:
                params = scipy.dot(scipy.inv(scipy.array([self.norm.unit,
                                                          surface.meri.unit,
                                                          surface.sagi.unit])),
                                   (self-surface.vec).x())

                if surface.edgetest(params[1],params[2]):
                    return params[0]
                else:
                    return []

            except ValueError:
                print('no?')
                return []
        else:
            return []

    def x0(self):
        """returns cartesian coordinate along first dimension

        Returns:
           numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[0] + self.norm.s*self.norm.unit[0]

    def x1(self):
        """returns cartesian coordinate along second dimension

        Returns:
            numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[1] + self.norm.s*self.norm.unit[1]

    def x2(self):
        """returns cartesian coordinate along third dimension

        Returns:
            numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[2] + self.norm.s*self.norm.unit[2]

    def c(self):
        """Conversion of vector to opposite coordinate system

        Returns:
            copy of vector object with opposite coordinate system
            (set with .flag parameter)

        """
        return (self + self.norm).c()

    def rmin(self):
        """rmin returns the s value along the norm vector which minimizes
        the r0() value (the closest position to the origin norm axis

        Returns:
            numpy array of s values in meters
        """
        return -1*self.s*(self.unit[0]*self.norm.unit[0] - 
                          self.unit[1]*self.norm.unit[1]
                          )/(self.norm.unit[0]**2+self.norm.unit[1]**2) 

    def r(self):
       """return cylindrical coordinate values

        Returns:
            numpy array of cylindrical coordinates in meters and radians

        """
        return (self + self.norm).r()
    
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


#def multiBeam(surf1,surf2,x=10,y=10):
#    """ This function provides the capability of generating multiple beams
#    from two surfaces using the split functionality.  The total number of
#    beams goes as (x*y)^2, thus at very small numbers of x,y the code slows
#    abruptly"""#
#
#    temp1 = surf1.split(x,y)
#    temp2 = surf2.split(x,y)
#    for 
