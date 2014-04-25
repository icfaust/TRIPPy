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

    def trace3(self, plasma, limiter=0):
        """ extends the norm vector into finding viewing length in vessel"""
        temp = [0.]
        # norm vector is modfied following the convention set in geometry.Origin
        if not self._origin is plasma:
            self.redefine(plasma)

        invesselflag = plasma.inVessel(self.r()[...,-1])
        vessel = plasma.eq.getMachineCrossSection()
        
        intersect = _beam.interceptCyl(self.x()[...,-1],self.norm.unit,vessel[0],vessel[1])
        if scipy.isfinite(intersect):
            self.norm.s = scipy.append(self.norm.s,intersect)

        intersect = _beam.interceptCyl(self.x()[...,-1],self.norm.unit,vessel[0],vessel[1])
        if not invesselflag and scipy.isfinite(intersect):
            self.norm.s = scipy.append(self.norm.s,intersect)
        
        # This is used when the actual plasma vessel structure is not as described in the eqdsk
        # example being the Limiter on Alcator C-Mod, which then keys to neglect an intersection,
        # and look for the next as the true wall intersection.
        for i in xrange(limiter):
            intersect = _beam.interceptCyl(self.x()[...,-1],self.norm.unit,vessel[0],vessel[1])
            self.norm.s[-1] = intersect

    def trace2(self, plasma, step=1e-2):
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

    def trace(self,plasma,eps=1e-4):
        """ finds intercepts with vessel surfaces assuming that the vessel is toroidal"""
        pt1 = self(0).r() # pull r,z of diode (pt1)
        pntinves = plasma.inVessel(pt1) # test if invessel

        s = eps #self.norm.s
        dels = s
        vessel = plasma.eq.getMachineCrossSection()
        trig = True

        while trig:

            pt2 = self(s).r() # pull vec of 'aperature' (pt2)
            dels *= _beam.intercept2d(pt1,pt2,vessel[0],vessel[1])-1

            # find at least 2 intercepts
            s += dels
            pt1 = pt2
            
            if abs(dels) < eps:
                trig = False

        if not pntinves:

            s1 = s.copy()
            s += eps
            dels = eps

            pt1 = self(s).r()
            trig = True

            while trig:

                pt2 = self(s+eps).r() # pull vec of 'aperature' (pt2)                
                dels *= _beam.intercept2d(pt1,pt2,vessel[0],vessel[1])-1
                
                # find at least 2 intercepts
                s += dels
                pt1 = pt2
            
                if abs(dels) < eps:
                    trig = False

            self.norm.s = scipy.array([0,s1,s])
        else:
            self.norm.s = scipy.array([0,s])
            
    def _intercept(self,pt1,pt2,outline,flag=False):
        r = outline[0]
        z = outline[1]
        delr = pt1[0] - pt2[0]
        delz = pt1[2] - pt2[2]
        params = scipy.zeros((2,len(r)))

        for i in xrange(len(r)-1):
            #mat = scipy.array([[delr, r[i+1]-r[i]],
            #                   [delz, z[i+1]-z[i]]])


            a = delr
            b = r[i+1]-r[i]
            c = delz
            d = z[i+1]-z[i]

            invmat = scipy.array([[d, -b], [-c, a]])/(a*d - b*c)
            params[:,i] = scipy.dot(invmat,scipy.array([pt1[0]-r[i],pt1[2]-z[i]]))

            #params[:,i] = scipy.dot(scipy.linalg.inv(mat),
            #                      scipy.array([pt1[0]-r[i],pt1[2]-z[i]]))

            #probably need a try catch for bad inverses.
        goods = params[0][scipy.logical_and(params[1] > 0,params[1] < 1)]  #index where there are actual intercepts  
        return goods[goods > 0].min()


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
    r"""Generates a Beam vector object assuming macroscopic 
        
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


    def trace(self, plasma, eps=1e-4):
        """ finds intercepts with vessel surfaces assuming that the vessel is toroidal"""
        pt1 = self(0).r() # pull r,z of diode (pt1)
        pntinves = plasma.inVessel(pt1) # test if invessel

        s = eps #self.norm.s
        dels = s
        vessel = plasma.eq.getMachineCrossSection()
        trig = True

        while trig:

            pt2 = self(s).r() # pull vec of 'aperature' (pt2)
            dels *= _beam.intercept2d(pt1,pt2,vessel[0],vessel[1])-1

            # find at least 2 intercepts
            s += dels
            pt1 = pt2
            
            if abs(dels) < eps:
                trig = False

        if not pntinves:

            s1 = s.copy()
            s += eps
            dels = eps

            pt1 = self(s).r()
            trig = True

            while trig:

                pt2 = self(s+eps).r() # pull vec of 'aperature' (pt2)                
                dels *= _beam.intercept2d(pt1,pt2,vessel[0],vessel[1])-1
                
                # find at least 2 intercepts
                s += dels

                pt1 = pt2
            
                if abs(dels) < eps:
                    trig = False

            self.norm.s = scipy.array([0,s1,s])
        else:
            self.norm.s = scipy.array([0,s])
    
    def _intercept(self,pt1,pt2,outline):
        r = outline[0]
        z = outline[1]
        delr = pt1[0] - pt2[0]
        delz = pt1[2] - pt2[2]
        params = scipy.zeros((2,len(r)))

        for i in xrange(len(r)-1):
            #mat = scipy.array([[delr, r[i+1]-r[i]],
            #                   [delz, z[i+1]-z[i]]])

            a = delr
            b = r[i+1]-r[i]
            c = delz
            d = z[i+1]-z[i]

            #e = pt1[0]-r[i]
            #f = pt1[2]-z[i]


            #params[:,i] = scipy.array([d*e - b*f,a*f - c*e])/(a*d - b*c)
            invmat = scipy.array([[d, -b], [-c, a]])/(a*d - b*c)
            params[:,i] = scipy.dot(invmat,scipy.array([pt1[0]-r[i],pt1[2]-z[i]]))


            #params[:,i] = scipy.dot(scipy.linalg.inv(mat),
            #                      scipy.array([pt1[0]-r[i],pt1[2]-z[i]]))

            #probably need a try catch for bad inverses.
        goods = params[0][scipy.logical_and(params[1] > 0,params[1] < 1)]  #index where there are actual intercepts  
        return goods[goods > 0].min()

    def trace3(self, plasma, limiter=0):
        """ extends the norm vector into finding viewing length in vessel"""
        temp = [0.]
        # norm vector is modfied following the convention set in geometry.Origin
        if not self._origin is plasma:
            self.redefine(plasma)

        invesselflag = plasma.inVessel(self.r()[...,-1])
        vessel = plasma.eq.getMachineCrossSection()
        
        intersect = _beam.interceptCyl(self.x()[...,-1],self.norm.unit,vessel[0],vessel[1]) + self.norm.s[-1]
        if scipy.isfinite(intersect):
            self.norm.s = scipy.append(self.norm.s,intersect)

        intersect = _beam.interceptCyl(self.x()[...,-1],self.norm.unit,vessel[0],vessel[1]) + self.norm.s[-1]
        if not invesselflag and scipy.isfinite(intersect):
            self.norm.s = scipy.append(self.norm.s,intersect)
        
        # This is used when the actual plasma vessel structure is not as described in the eqdsk
        # example being the Limiter on Alcator C-Mod, which then keys to neglect an intersection,
        # and look for the next as the true wall intersection.
        for i in xrange(limiter):
            intersect = _beam.interceptCyl(self.x()[...,-1],self.norm.unit,vessel[0],vessel[1]) + self.norm.s[-1]
            self.norm.s[-1] = intersect

    def trace2(self, plasma, step=1e-2):
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


#def multiBeam(surf1,surf2,x=10,y=10):
#    """ This function provides the capability of generating multiple beams
#    from two surfaces using the split functionality.  The total number of
#    beams goes as (x*y)^2, thus at very small numbers of x,y the code slows
#    abruptly"""#
#
#    temp1 = surf1.split(x,y)
#    temp2 = surf2.split(x,y)
#    for 
