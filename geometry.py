import scipy

class Hat(Object):
    """ explicitly just unit vector without error
    which is then defined for the various coordinate
    systems based on classes based on this, a Hat class
    never interacts with another Hat class, only one 
    with a defined coordinate system"""
    def __init__(self, x_hat):

        self.unit = scipy.array(x_hat)

    def _cross(self):
        " matrix necessary for a cross-product calculation"""
        return  scipy.array(((0,-self.unit[2],self.unit[1]),
                              (self.unit[2],0,-self.unit[0]),
                              (-self.unit[1],self.unit[0],0)))


class Vecx(Hat):
    """ explicitly a cartesian unit vector, but can be set as 
    a cylindrical unit vector by setting the flag to true, all
    vector math defaults to first vector"""
    flag = False
    
    def __init__(self, x_hat, s=[]):
        #if r is specified, it is assumed that x_hat has unit length
        if not s:
            s = scipy.sqrt(scipy.sum(x_hat**2))
            x_hat /= s
        super(Vecx,self).__init__(x_hat)
        self.s = SP.array(s)

    def __add__(self, Vec):
        """ Vector addition, convert to cartesian
        add, and possibly convert back"""
        if Vec.flag:
            Vec = Vec.c()
        return Vecx(self.x() + Vec.x())

    def __sub__(self, Vec):
        """ Vector subtraction """
        if Vec.flag:
            Vec = Vec.c()
        return Vecx(self.x() - Vec.x())

    def __mul__(self, Vec):
        """ Dot product """
        if Vec.flag:
            Vec = Vec.c()
        return self.unit.T*Vec.unit

    def c(self):
        """ convert to cylindrical coord """
        return Vecr((scipy.sqrt(self.unit[0]**2+self.unit[1]**2),
                     scipy.arctan2(self.unit[1],self.unit[0]),
                     scipy.unit[2]),
                    s=s)
    
    def x0(self):
        return self.s*self.hat[0]
    
    def x1(self):
        return self.s*self.hat[1]
    
    def x2(self):
        return self.s*self.hat[2]
    
    def x(self):
        return self.s*self.hat

    def point(self,ref,err=[]):
        return Point(self.x(),ref,err=err)
                   
class Vecr(Hat):
    """ explicitly a cylindrical unit vector, but can be set as 
    a cartesian unit vector by using the c call, all
    vector math defaults to first vector"""
    flag = True

    def __init__(self, x_hat, s=[]):
        if not s:
            s = scipy.sqrt(x_hat[0]**2 + x_hat[2]**2)
            x_hat[0] /= s
            x_hat[2] /= s
        super(Vecr,self).__init__(x_hat)
        self.s = SP.array(s)

    def __add__(self, Vec):
        """ Vector addition, convert to cartesian
        add, and possibly convert back"""
        if Vec.flag:
            Vec = Vec.c()
        return Vecr(self.c().x() + Vec.x())

    def __sub__(self, Vec):
        """ Vector subtraction """
        if Vec.flag:
            Vec = Vec.c()
        return Vecr(self.c().x() - Vec.x())
        
    def __mul__(self, Vec):
        """ Dot product """
        if not Vec.flag:
            Vec = Vec.c()
        return self.unit[0]*Vec.unit[0]*scipy.cos(self.unit[1]-Vec.unit[1]) + self.unit[2]*Vec.unit[2]

    def c(self):
        """ convert to cartesian coord """
        return Vecx((self.unit[0]*scipy.cos(self.unit[1]),
                     self.unit[0]*scipy.sin(self.unit[1]),
                     self.unit[2]),
                    s=s)
    
    def x0(self):
        return self.s*self.hat[0]
    
    def x1(self):
        return self.hat[1]
    
    def x2(self):
        return self.s*self.hat[2]
    
    def x(self):
        return scipy.array((self.x0,self.x1,self.x2))
        
    def point(self,ref,err=[]):
        return Point(self.x(),ref,flag=True)

def angle(Vec1, Vec2):
    return scipy.arccos(Vec1 * Vec2) 

def cross(Vec1, Vec2):
    if Vec1.flag == Vec2.flag:
        x_hat = scipy.dot(Vec1._cross(),Vec2.unit)      
    else:
        x_hat = scipy.dot(Vec1._cross(),Vec2.c().unit)
            
    if Vec1.flag:
        return Vecr(x_hat, s=Vec1.s*Vec2.s)
    else:
        return Vecx(x_hat, s=Vec1.s*Vec2.s)

class Point(Object):
    """ a point class can only be defined relative to an origin,
    there will be an additional point class which will be known
    as grid which reduces the redundant reference to origin
    and depth for memory savings, and will order points in 
    such a way for easier calculation."""
    def __init__(self, x_hat, ref, err=[],flag=False):
        
        self.x = scipy.array(x_hat)
        if err:
            self.error = err

        self.flag = flag
        self._origin = ref
        self._depth = ref._depth + 1 # basis origin is depth = 0
         
    def Vec(self, c=False):
        """ c provides the ability to convert coordinate
        systems"""
        if c == self._origin.flag:
            return Vecx(self.x)
        else:
            return Vecr(self.x)

    def redefine(self, neworigin):
        """ changes depth of point by calculating with respect to a new
        origin, for calculations with respect to a flux grid, etc. this
        should reduce caluclation substantially."""
        
        # use _lca to find common ancestor and return tree to common
        org1,org2 = self._lca(neworigin)
        
        # loop over the first 'path' of the current point
        temp = org1[0].vec()
        for idx in range(len(org1)-1):
            temp = temp + org1[idx]
        self.x =         
        self._origin = neworigin
        self._depth = neworigin._depth + 1

    def _genOriginsToParent(self, depth=self._depth):
        """ generate a list of points which leads to the overall basis
        origin of the geometry, the number of elements will be the depth"""
        temp = self.ref
        pnts = depth*[0]
        
        for idx in range(depth):
            temp = temp._origin
            pnts[idx] = temp

        return pnts

    def _lca(self, point2):
        """ recursively solve for the common point of reference
        between two points as given by a depth number starting
        from the base of the tree. It will return a tuple which
        contains the nodes leading to the common point."""
        
        temp = True
        idx = -1
        
        pt1 = self._getOriginsToParent()
        pt2 = point2._getOriginsToParent()
        
        # determine the shorter origins list
        if len(pt1) > len(pt2):
            lim = len(pt2) - 1
        else:
            lim = len(pt1) - 1
            
        # compare origins from the base (which should be the same)
        # and when they don't match store the index.  Return the
        # negative index which corresponds to the last match
        while temp:

            if (pt1[idx] is pt2[idx]) and (lim + idx):
                idx -= 1
            else:
                idx += 1
                temp = False

        return (pt1[idx:],pt2[idx:])

class Grid(Point):
    
    def __init__(self, x_hat, ref, err=scipy.matrix((0,0,0))):
        """ a grid compartmentalizes a large set of points which
        are easily defined on a regular grid. Unlike points, grids
        points cannot change reference frames.  While this might
        be okay for a flat plane of points, for shapes such as 
        spheres, parabolas, ellipsoids and other complicated shapes
        would not be easily defined.  When a reference frame change
        is required, it will revert to a similarly sized array
        of Point Objects, which use order mxn more memory"""


class Origin(Point):

    def __init__(self, x_hat, ref, Vec=[], err=scipy.matrix((0,0,0)), angle=[]):
        """ an Origin is defined by a point and two vectors.
        The two vectors being: 1st the normal to the surface,
        principal axis, z-vector or the (0,0,1) vector of the
        new system defined in the reference system . The 
        second vector along with the first fully defines 
        meridonial ray paths, x-vector or the (1,0,0) vector.
        The sagittal ray path, y-vector or the (0,1,0)
        is defined through a cross product.  Point position
        and rotation matricies are stored at instantiation.

        If the angles alpha, beta, and gamma are specified 
        following the eulerian rotation formalism, it is
        processed in the following manner: alpha is the 
        rotation from the principal axis in the meridonial
        plane, beta is the rotation about the plane normal
        to the meridonial ray, or 2nd specified vector,
        and gamma is the 2nd rotation about the principal
        axis.  This might change based on what is most
        physically intuitive."""
        # test Vec1 and Vec2 for ortho-normality

        if Vec:
            if not Vec[0] * Vec[1]: 
                # generate point based off of previous origin
                super(Origin,self).__init__(x_hat, ref, err=err)

                # generate rotation matrix based off coordinate system matching (this could get very interesting)
                self.rot = scipy.matrix((Vec[0].unit.T,
                                         cross(Vec[0],Vec[1]).unit.T,
                                         Vec[1].unit.T))

        elif angle:
            super(Origin,self).__init__(x_hat, ref, err=err)
            a = scipy.array(angle[0])
            b = scipy.array(angle[1])
            g = scipy.array(angle[2])

            # taken from the wikipedia convention, which will allow for easier modification.
            # to make this extensive rather than intensive might require me flipping alpha
            # and gamma.  This is to be seen when testing.
            # https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
            c1 = scipy.cos(a)
            c2 = scipy.cos(b)
            c3 = scipy.cos(g)
            s1 = scipy.sin(a)
            s2 = scipy.sin(b)
            s3 = scipy.sin(g)

            self.rot = scipy.matrix(((c1*c3 - c2*s1*s3, -c1*s3 - c2*c3*s1, s1*s2),
                                     (c3*s1 + c1*c2*s3, c1*c2*c3 - s1*s3, -c1*s2),
                                     (s2*s3, c3*s2, c2)))
        else:
            raise ValueError
            #throw error here

class Center(Origin):
    """ this is the class which underlies all positional calculation.
    It is located at (0,0,0) and is inherently a cylindrical coordinate
    system.  It is from the translation of inherently cylindrical data
    into toroidal coordinates requires this rosetta stone"""
