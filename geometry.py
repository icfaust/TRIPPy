import scipy

class Hat(object):
    """ explicitly just unit vector without error
    which is then defined for the various coordinate
    systems based on classes based on this, a Hat class
    never interacts with another Hat class, only one 
    with a defined coordinate system"""
    def __init__(self, x_hat):

        self.unit = scipy.array(x_hat)
        if self.unit.shape[0] != 3:
            raise ValueError

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
        xin = scipy.array(x_hat,dtype=float)
        if not len(s):
            s = scipy.sqrt(scipy.sum(xin**2,axis=0))  # this needs to be modified s
            xin /= s

        super(Vecx,self).__init__(xin)
        self.s = scipy.atleast_1d(s)

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

    def __getitem__(self,idx):
        return self.x()[idx]

    def c(self):
        """ convert to cylindrical coordinates """
        return Vecr((scipy.sqrt(self.unit[0]**2+self.unit[1]**2),
                     scipy.arctan2(self.unit[1],self.unit[0]),
                     self.unit[2]),
                    s=self.s)
    
    def x0(self):
        return self.s*self.unit[0]
    
    def x1(self):
        return self.s*self.unit[1]
    
    def x2(self):
        return self.s*self.unit[2]
    
    def x(self):
        return scipy.array([self.x0(),
                            self.x1(),
                            self.x2()])

    def point(self,ref,err=[]):
        return Point(self.x(), ref, err=err)
                   
class Vecr(Hat):
    """ explicitly a cylindrical unit vector, but can be set as 
    a cartesian unit vector by using the c call, all
    vector math defaults to first vector"""
    flag = True

    def __init__(self, x_hat, s=[]):

        xin = scipy.array(x_hat, dtype=float)

        if not len(s):
            s = scipy.sqrt(x_hat[0]**2 + x_hat[2]**2)
            xin[0] /= s
            xin[2] /= s
        super(Vecr,self).__init__(xin)
        self.s = scipy.atleast_1d(s)

    def __add__(self, Vec):
        """ Vector addition, convert to cartesian
        add, and possibly convert back"""
        if Vec.flag:
            Vec = Vec.c()
        
        add = self.c().x() + Vec.x() # computed in cartesian coordinates
        add[1] = scipy.arctan2(add[1], add[0]) #convert to angle
        add[0] = scipy.absolute(add[0]/scipy.cos(add[1])) # generate radius
        return Vecr(add)

    def __sub__(self, Vec):
        """ Vector subtraction """
        if Vec.flag:
            Vec = Vec.c()

        sub = self.c().x() - Vec.x()        
        sub[1] = scipy.arctan2(sub[1], sub[0])
        sub[0] = scipy.absolute(sub[0]/scipy.cos(sub[1]))
        return Vecr(sub) 
        
    def __mul__(self, Vec):
        """ Dot product """
        if not Vec.flag:
            Vec = Vec.c()
        return self.unit[0]*Vec.unit[0]*scipy.cos(self.unit[1]-Vec.unit[1]) + self.unit[2]*Vec.unit[2]

    def __getitem__(self,idx):
        return self.x()[idx]

    def c(self):
        """ convert to cartesian coord """
        return Vecx((self.unit[0]*scipy.cos(self.unit[1]),
                     self.unit[0]*scipy.sin(self.unit[1]),
                     self.unit[2]),
                    s=self.s)
    
    def x0(self):
        return self.s*self.unit[0]
    
    def x1(self):
        return scipy.ones(self.s.shape)*self.unit[1]
    
    def x2(self):
        return self.s*self.unit[2]
    
    def x(self):
        return scipy.array([self.x0(),
                            self.x1(),
                            self.x2()])
        
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

class Point(object):
    """ a point class can only be defined relative to an origin,
    there will be an additional point class which will be known
    as grid which reduces the redundant reference to origin
    and depth for memory savings, and will order points in 
    such a way for easier calculation."""
    def __init__(self, x_hat, ref, err=[]):
        
        self._x = scipy.array(x_hat)
        if len(err):
            self.error = err
            
        self._origin = ref
        self._depth = ref._depth + 1 # basis origin is depth = 0
         
    def Vec(self, c=False):
        """ c provides the ability to convert coordinate
         systems"""
        if c == self._origin.flag:
            return Vecx(self._x)
        else:
            return Vecr(self._x)

    def redefine(self, neworigin):
        """ changes depth of point by calculating with respect to a new
        origin, for calculations with respect to a flux grid, etc. this
        should reduce caluclation substantially."""
        
        # use _lca to find common ancestor and return tree to common
        lca = self._lca(neworigin)
        
        # this will allows for the matrix math of the rotation to behave properly
        self._translate(lca)

    def _translate(self, lca):
        org = lca[0]
        orgnew = lca[1]
        shape = self._x.shape
        self._x = self._x.reshape(3,self._x.size/3)

        # loop over the first 'path' of the current point
        temp = self.Vec()
        for idx in range(len(org)-1,-1,-1):
            # a origin's point is defined by its recursive coordinate system
            # thus the value addition is not occuring in current origin system
            # put in fact the coordinate system that defines the origin.
            # thus this requires using the rotation matrix of the current origin
            # system to define it in the 'parent' coordinate system
            temp = org[idx].Vec() + org[idx].rot(temp)

        for idx in range(len(orgnew)):
            # the arot allows for translating into the current coordinate system
            temp = orgnew[idx].arot(temp - orgnew[idx].Vec())

        temp = neworigin.arot(temp - neworigin.Vec())
        # what is the vector which points from the new origin to the point?
        self._origin = neworigin
        self._depth = neworigin._depth + 1

        # convert vector to proper coordinate system matching new origin and save
        # arot forces the coordinate system to that of the new origin
        self._x = temp.x().reshape(shape)
        
    def x(self):
        """ heavily redundant, but will smooth out variational differences
        from the grid function"""
        return self._x

    def _genOriginsToParent(self):
        """ generate a list of points which leads to the overall basis
        origin of the geometry, the number of elements will be the depth"""
        temp = self._origin
        pnts = self._depth*[0]

        # as index increases, the closer to the point it becomes.
        for idx in range(self._depth-1,-1,-1):
            pnts[idx] = temp
            temp = temp._origin

        return pnts

    def _lca(self, point2):
        """ recursively solve for the common point of reference
        between two points as given by a depth number starting
        from the base of the tree. It will return a tuple which
        contains the nodes leading to the common point."""
        
        temp = True
        idx = 0
        pt1 = self._genOriginsToParent()
        pt2 = point2._genOriginsToParent()

        # determine the shorter origins list
        if self._depth > point2._depth:
            lim = point2._depth
        else:
            lim = self._depth
            

        # compare origins from the base (which should be the same)
        # and when they don't match store the index.  Return the
        # negative index which corresponds to the last match
        while temp:
            if not (lim - idx):
                idx += 1 # this takes care of ambiguity associated with the centerpoint
                temp = False
            elif (pt1[idx] is pt2[idx]):                
                idx += 1
            else:
                temp = False

        print(idx)
        return (pt1[idx:],pt2[idx:])


class Grid(Point):
    
    def __init__(self, x_hat, ref, mask=(False,False,False), err=scipy.array((0,0,0))):
        """ a grid compartmentalizes a large set of points which
        are easily defined on a regular grid. Unlike points, grids
        points cannot change reference frames.  While this might
        be okay for a flat plane of points, for shapes such as 
        spheres, parabolas, ellipsoids and other complicated shapes
        would not be easily defined.  When a reference frame change
        is required, it will revert to a similarly sized array
        of Point Objects, which use order mxn more memory"""
        
        # x_hat is expected to be a tuple containing three entities.
        # the entities are described with the mask variable which provides
        # a basic descriptor of the functional dependence.  At most,
        # only two variables may be dependent on the third.
        self._x = x_hat
        
        # mask provides an understanding of the nature of the grid, whether
        # it is planar or follows a funciton in a specific dimension this
        # should allow for various shapes to easily be implemented
        self._mask = mask

    def redefine(self,neworigin):
        """ redefine will break the simplicity of the grid, have it spit back
        a tuple of Points at memory cost"""
        raise ValueError

    def x(self):

        # access the mask, and determine if the variable is a function.

        # it is a ssumed that the masked variables are functions of the
        # other inputs. 
        return self._x # its not this simple


class Origin(Point):

    def __init__(self, x_hat, ref, Vec=[], err=scipy.array((0,0,0)), angle=[], flag=[]):
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
                self.norm = Vec[1]
                self.saggital = Vec[0]
                temp = cross(self.norm,self.saggital).unit.T
                # generate rotation matrix based off coordinate system matching (this could get very interesting)

        elif len(angle):
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


            self.norm = Vecx((c1*c3 - c2*s1*s3, -c1*s3 - c2*c3*s1, s1*s2),s=1)
            self.saggital =  Vecx((s2*s3, c3*s2, c2),s=1)
            temp = (c3*s1 + c1*c2*s3, c1*c2*c3 - s1*s3, -c1*s2)
        else:
            raise ValueError
            #throw error here
        self._rot = scipy.array((self.norm.unit.T,
                                 temp,
                                 self.saggital.unit.T))
        # set to coordinate system only if specified. otherwise inherit based on reference
        if flag:
            self.flag = flag
        else:
            self.flag = ref.flag


    def redefine(self,neworigin):
        """ rotation matrix also needs to be redefined """

        lca = self._lca(neworigin)
        super(Point,self)._translate(lca)
        self._rotate(lca)

    def _rotate(self,lca):
        temp = self._rot
        for idx in range(len(org)-1,-1,-1):
            # change the _rot coordinates to accurately reflect all of the necessary variation.
            temp = org[idx].rot(temp)

        for idx in range(len(orgnew)):
            # the arot allows for translating into the current coordinate system
            temp = orgnew[idx].arot(temp)

        self._rot = neworigin.arot(temp)


    def rot(self,vec):
        if not self._origin.flag == vec.flag:
            vec = vec.c()

        if self.flag:
            return Vecr(scipy.dot(self._rot,vec.unit),s=vec.s)
        else:
            return Vecx(scipy.dot(self._rot,vec.unit),s=vec.s)
                
    def arot(self,vec):
        if not self.flag == vec.flag:
            vec = vec.c()

        if self.flag:
            return Vecr(scipy.dot(self._rot.T,vec.unit),s=vec.s)
        else:
            return Vecx(scipy.dot(self._rot.T,vec.unit),s=vec.s)

class Center(Origin):
    """ this is the class which underlies all positional calculation.
    It is located at (0,0,0) and is inherently a cylindrical coordinate
    system (unless flag set otherwise). 
    It is from the translation of inherently cylindrical data
    into toroidal coordinates requires this rosetta stone, it can
    be dynamically set to becoming an origin given a specification
    of another origin."""


    def __init__(self, flag=True):
        self._x = scipy.array((0.,0.,0.))
        self._depth = 0
        self._origin = []
        self._rot = scipy.eye(3)
        self.flag = flag
        # could not use super due to the problem in defining the value of 
        # the depth.  This is simple, though slightly redundant.
        # large number of empty values provide knowledge that there are no
        # lower references or rotations to this, the main coordinate system

    def Vec(self,c=False):
        """ c provides the ability to convert coordinate
        systems"""
        if c == self.flag:
            return Vecx(self._x,s=1)
        else:
            return Vecr(self._x,s=1)

    def rot(self,vec):

        if self.flag:
            return Vecr(vec.unit,s=vec.s)
        else:
            return Vecx(vec.unit,s=vec.s)
