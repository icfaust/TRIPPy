import geometry
import surface
import scipy
import scipy.linalg

class Ray(geometry.Point):

    def __init__(self, pt1, inp2, err=[]):

        try:
            self.norm = geometry.pts2Vec(pt1, inp2)
        except AttributeError:
            self.norm = inp2
            
        self._start = None
        self._end = None
        super(Ray,self).__init__(pt1, err=err)

    def x(self):
        return (self + self.norm).x()

    def r(self):
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

    def redefine(self, neworigin):
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

    def trace(self, plasma, step=1e-2):
        """ extends the norm vector into finding viewing length in vessel"""
        
        # norm vector is modfied following the convention set in geometry.Origin
        if not self._origin is plasma:
            self.redefine(plasma)


        
        end = plasma.eq.getMachineCrossSection()[0].max() + self.s
        #step = end/1e2


        temp = self(scipy.mgrid[0:end:step]).r() #start from diode, and trace through to find when it hits the vessel
        sin = temp.size

        # set if in cylindrical coordinates
        idx = 1
        invesselflag = plasma.inVessel(temp[:,0])

        # if diode is not invessel, find when the view is in vessel
        if not invesselflag:
            flag = True
            while flag:

                pntinves = plasma.inVessel(temp[:,idx])
                idx += 1

                if pntinves or idx + 1 == sin:
                    flag = False
                invesselflag = pntinves
            self._start = idx*step

        # find point at which the view escapes vessel
        if invesselflag:
            flag = True
            while flag:
            
                pntinves = plasma.inVessel(temp[:,idx])
                idx += 1

                if not pntinves or idx + 1 == sin:
                    flag = False
                invesselflag = pntinves
            self._end = idx*step
        
        self.norm.s = scipy.array([0])
        if self._start:
            self.norm.s = scipy.append(self.norm.s, self._start)
        if self._end:
            self.norm.s = scipy.append(self.norm.s, self._end)


    def trace2(self,plasma,eps=1e-4):
        """ finds intercepts with vessel surfaces assuming that the vessel is toroidal"""
        pt1 = self(0).r() # pull r,z of diode (pt1)
        pntinves = plasma.inVessel(pt1) # test if invessel

        s = self.norm.s
        dels = 1e2
        norm = self(s).r() # pull r,z of 'aperature' (pt2)
        
        while abs(dels) > eps:

            # find at least 2 intercepts
            s = s + dels

        print('not implemented yet')

    def _intercept(self,pt1,plasma):
        print('not implemented yet')

    def intercept(self,surface):
        if self._origin is surface._origin:
            try:
                params = scipy.dot(scipy.linalg.inv(scipy.array([self.norm.unit,
                                                                 surface.meri.unit,
                                                                 surface.sagi.unit])),
                                   (self-surface).x())

                if surface.edgetest(params[1],params[2]):
                    return params[0]
                else:
                    return []

            except ValueError:
                print('no?')
                return []
        else:
            return []

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
    """ generates an origin with defined etendue and non scalar.
    Beams propagation is assumed to be non-refractive, such that
    the variation in the nature of the etendue analytical. 
    """

    def __init__(self, surf1, surf2, err=[]):
        """ the beam generates a normal vector which points from
        surface 1 to surface 2.  After which the sagittal and 
        meridonial rays are simply the surface 1 sagittal and
        meridonial rays with respect """

        normal = geometry.pts2Vec(surf1, surf2)
        #orthogonal coordinates based off of connecting normal

        snew = surf1.sagi - normal*((surf1.sagi * normal)*(surf1.sagi.s/normal.s))
        mnew = surf1.meri - normal*((surf1.meri * normal)*(surf1.meri.s/normal.s))
        super(Beam, self).__init__(surf1, surf1._origin, Vec = [mnew,normal], err=err)
        #calculate area at diode.
        self.sagi.s = snew.s
        a1 = surf1.area(snew.s,mnew.s)

        #calculate area at aperature
        a2 = surf2.area((((self.sagi*surf2.sagi)/self.sagi.s)**2 + ((self.meri*surf2.sagi)/self.meri.s)**2)**.5,
                        (((self.sagi*surf2.meri)/self.sagi.s)**2 + ((self.meri*surf2.meri)/self.meri.s)**2)**.5)

        #generate etendue
        self.etendue = a1*a2/(normal.s ** 2)
        self._start = scipy.array(0)
        self._end = []
       
    def trace(self, plasma, step=1e-2):
        """ extends the norm vector into finding viewing length in vessel"""

        end = plasma.eq.getMachineCrossSection()[0].max() + self.s
        

        self.redefine(plasma)

        temp = self(scipy.mgrid[0:end:step]).r() #start from diode, and trace through to find when it hits the vessel
        sin = temp.size

        # set if in cylindrical coordinates
        idx = 1
        invesselflag = plasma.inVessel(temp[:,0])
        # if diode is not invessel, find when the view is in vessel
        if not invesselflag:
            flag = True
            while flag:

                pntinves = plasma.inVessel(temp[:,idx])
                idx += 1
                if pntinves or idx + 1 == sin:
                    flag = False
                invesselflag = pntinves

            self._start = self.norm.s[idx]

        # find point at which the view escapes vessel
        if invesselflag:
            flag = True
            while flag:
            
                pntinves = plasma.inVessel(temp[:,idx])
                idx += 1
                if not pntinves or idx + 1 == sin:
                    flag = False
                invesselflag = pntinves
            self._end = self.norm.s[idx]
        
        self.norm.s = scipy.array([0])
        if self._start:
            self.norm.s = scipy.append(self.norm.s, self._start)
        if self._end:
            self.norm.s = scipy.append(self.norm.s, self._end)


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

    def x(self):
        return (self + self.norm).x()

    def c(self):
        return (self + self.norm).c()

    def r(self):
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
