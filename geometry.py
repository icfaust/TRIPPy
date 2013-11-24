import scipy,warnings

one = scipy.array([1.0])

def unit(x_hat):
    """ explicitly just unit vector without error
    which is then defined for the various coordinate
    systems based on classes based on this, a Hat class
    never interacts with another Hat class, only one 
    with a defined coordinate system"""
    temp = scipy.squeeze(x_hat)
 
    if len(temp) != 3 or temp.max() > 1 or temp.min() < -1:
        raise ValueError
    return temp


def Vecx(x_hat, s=None):
    """ explicitly a cartesian unit vector, but can be set as 
    a cylindrical unit vector by setting the flag to true, all
    vector math defaults to first vector"""
    flag = False
    xin = scipy.array(x_hat, dtype=float)

    if s is None:
        s = scipy.sqrt(scipy.sum(xin**2, axis=0)) 
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            xin = scipy.where(s == 0, xin, xin/s)

    return Vec(xin, s, flag=flag)


def Vecr(x_hat, s=None):
    """ explicitly a cylindrical unit vector, but can be set as 
    a cartesian unit vector by using the c call, all
    vector math defaults to first vector"""

    flag = True    
    xin = scipy.array(x_hat, dtype=float)
        
    if s is None:
        s = scipy.sqrt(x_hat[0]**2 + x_hat[2]**2)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning)
            xin[0] = scipy.where(s == 0, xin[0], xin[0]/s)
            xin[2] = scipy.where(s == 0, xin[2], xin[2]/s)
            
    return Vec((xin[0]*scipy.cos(xin[1]),
                xin[0]*scipy.sin(xin[1]),
                xin[2])
               ,s, flag=flag)       


class Vec(object):
    
    def __init__(self, x_hat, s, flag=False):
        """ inherently a cartesian set of coordinates given by unit,
        and a length set by s. The flag tells if it is a cartesian
        or cylindrial vector"""

        self.unit = unit(x_hat)
        self.s = scipy.squeeze(s)
        self.flag = flag

    def __add__(self, Vec):
        """ Vector addition, if to vectors go a walkin'
        the first does the talking and sets coordinate system"""
        x0 = self.x0() + Vec.x0()
        x1 = self.x1() + Vec.x1()
        x2 = self.x2() + Vec.x2()
        
        new = Vecx((x0,x1,x2))
        new.flag = self.flag
        return new

    def __sub__(self, Vec):
        """ Vector subtraction """
        x0 = self.x0() - Vec.x0()
        x1 = self.x1() - Vec.x1()
        x2 = self.x2() - Vec.x2()
        
        new = Vecx((x0,x1,x2))    
        new.flag = self.flag
        return new

    def __neg__(self):
        """ uniary minus"""
        return Vec(-self.unit,self.s)

    def __mul__(self, Vec):
        """ Dot product """
        try:
            return scipy.dot(self.unit.T,Vec.unit)
        except AttributeError:
            return Vec(self.unit,Vec*self.s)

    def __div__(self, val):
        return Vec(self.unit,self.s/val)

    def __getitem__(self,idx):
        if flag:
            return self.r()[idx]
        else:
            return self.x()[idx]

    def _cross(self):
        " matrix necessary for a cross-product calculation"""
        return  scipy.array(((0,-self.unit[2],self.unit[1]),
                             (self.unit[2],0,-self.unit[0]),
                             (-self.unit[1],self.unit[0],0)))

    def x0(self):
        """x0,x1,x2 is taking care of an array multiplication problem"""
        return self.s*self.unit[0]

    def x1(self):
        return self.s*self.unit[1]
    
    def x2(self):
        return self.s*self.unit[2]

    def c(self):
        return Vec(self.unit,self.s, flag = (not self.flag))

    def x(self):
        """ cartesian full vector"""
        return scipy.squeeze([self.x0(),
                              self.x1(),
                              self.x2()])
    
    def r(self):
        """ cylindrical full vector"""
        return scipy.squeeze([self.s*scipy.sqrt(self.unit[0]**2+self.unit[1]**2),
                              scipy.arctan2(self.unit[1],self.unit[0]),
                              self.s*self.unit[2]])

    def point(self,ref,err=[]):
        return Point(self, ref, err=err)
                   

def angle(Vec1, Vec2):
    """angle between two vectors"""
    return scipy.arccos(Vec1 * Vec2) 

def cross(Vec1, Vec2):
    """cross product"""
    new = Vec(scipy.dot(Vec1._cross(),Vec2.unit),Vec1.s*Vec2.s)
    new.flag = Vec1.flag
    return new

class Point(Vec):
    """ a point class can only be defined relative to an origin,
    there will be an additional point class which will be known
    as grid which reduces the redundant reference to origin
    and depth for memory savings, and will order points in 
    such a way for easier calculation."""
    def __init__(self, x_hat, ref, err=[]):        
        

        if type(x_hat) is Vec:
            temp = x_hat
        elif ref.flag:
            temp = Vecr(x_hat)
        else:
            temp = Vecx(x_hat)

        if len(err):
            self.err = err
            
        self._origin = ref
        self._depth = ref._depth + 1 # basis origin is depth = 0
        super(Point,self).__init__(temp.unit, temp.s, flag=ref.flag)
         
    def redefine(self, neworigin):
        """ changes depth of point by calculating with respect to a new
        origin, for calculations with respect to a flux grid, etc. this
        should reduce caluclation substantially."""
        
        # use _lca to find common ancestor and return tree to common
        lca = self._lca(neworigin)
        
        # this will allows for the matrix math of the rotation to behave properly
        self._translate(lca, neworigin)

    def _translate(self, lca, neworigin):
        org = lca[0]
        orgnew = lca[1]
        shape = self.unit.shape
        if len(shape) > 2:    
            sshape = self.s.shape
            self.s = self.s.flatten()
            self.unit = self.unit.reshape(3,self.unit.size/3)

        # loop over the first 'path' of the current point
        temp = Vec(self.unit,self.s,flag=self.flag)

        for idx in range(len(org)-1,-1,-1):
            # a origin's point is defined by its recursive coordinate system
            # thus the value addition is not occuring in current origin system
            # put in fact the coordinate system that defines the origin.
            # thus this requires using the rotation matrix of the current origin
            # system to define it in the 'parent' coordinate system
            
            temp = org[idx] + org[idx].rot(temp)
            # vector addition is really tricky

        for idx in range(len(orgnew)):
            # the arot allows for translating into the current coordinate system
            temp = orgnew[idx].arot(temp - orgnew[idx])

        # what is the vector which points from the new origin to the point?
        self._origin = neworigin
        self._depth = neworigin._depth + 1

        # convert vector to proper coordinate system matching new origin and save
        # arot forces the coordinate system to that of the new origin
        
        if len(shape) > 2:
            self.unit = temp.unit.reshape(shape)
            self.s = temp.s.reshape(sshape)
        else:
            self.unit = temp.unit
            self.s = temp.s

    def split(self, *args, **kwargs):
        obj = kwargs.pop('obj', None)
        if obj is None:
            obj = type(self)
        temp = self.x()
        if temp.size > 3:
            # initialize
            return fill(obj,temp[0],temp[1],temp[2], *args, **kwargs)

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
        pt2.append(point2)
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

        
        return (pt1[idx:],pt2[idx:])
    
    def c(self):
        raise NotImplementedError('point locked by reference coordinate system')


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

    def __init__(self, x_hat, ref, Vec=[], err=scipy.array((0,0,0)), angle=[], flag=None):
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
        super(Origin,self).__init__(x_hat, ref, err=err)
        if Vec:
            # generate point based off of previous origin

            self.norm = Vec[1]
            self.meri = Vec[0]

            self.sagi = cross(self.meri,self.norm)
            # generate rotation matrix based off coordinate system matching (this could get very interesting)

        elif len(angle):

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

            self.sagi = Vec((c3*c2, s3*c1 - c3*s2*s1, c3*s2*c1 + s3*s1), one)
            self.meri = Vec((-s3*c2, c3*c1 + s3*s2*s1, c3*s1 - s3*s2*c1), one)
            self.norm = Vec((-s2, -c2*s1, c2*c1), one)

        else:
            raise ValueError("rotation matrix cannot be specified without a normal"
                             " and meridonial ray, please specify either two"
                             " vectors or a set of euclidean rotation angles")
            #throw error here
        self._rot = [self.sagi.unit,
                     self.meri.unit,
                     self.norm.unit]
        # set to coordinate system only if specified. otherwise inherit based on reference
        if not flag is None:
            self.flag = flag
        else:
            self.flag = ref.flag


    def redefine(self, neworigin):
        """ rotation matrix also needs to be redefined """

        lca = self._lca(neworigin)
        super(Origin,self)._translate(lca, neworigin)
        self._rotate(lca, neworigin)


    def _rotate(self, lca, neworigin):
        """ rotates the fundamental vectors of the space"""
        org = lca[0]
        orgnew = lca[1]

        temp1 = self.norm
        temp2 = self.sagi
        mtemp = self.meri.s

        for idx in range(len(org)-1,-1,-1):
            # change the _rot coordinates to accurately reflect all of the necessary variation.
            temp1 = org[idx].rot(temp1)
            temp2 = org[idx].rot(temp2)

        for idx in range(len(orgnew)):
            # the arot allows for translating into the current coordinate system
            temp1 = orgnew[idx].arot(temp1)
            temp2 = orgnew[idx].arot(temp2)

        self.norm = temp1
        self.sagi = temp2
        self.meri = cross(self.norm, self.sagi)
        self.meri.s = mtemp

    def rot(self,vec):

        temp = Vec(scipy.dot(scipy.array(self._rot).T, vec.unit),vec.s)
        temp.flag = vec.flag
        return temp                
    
    def arot(self,vec):

        temp = Vec(scipy.dot(self._rot, vec.unit), vec.s)      
        temp.flag = vec.flag
        return temp

    def split(self, *args, **kwargs):
        
        try:
            obj = kwargs['obj']
            return super(Origin,self).split(*args,**kwargs)
        except TypeError:
            return super(Origin,self).split(self,self._origin, obj=type(self), Vec=[self.meri,self.norm], flag=self.flag)

class Center(Origin):
    """ this is the class which underlies all positional calculation.
    It is located at (0,0,0) and is inherently a cylindrical coordinate
    system (unless flag set otherwise). 
    It is from the translation of inherently cylindrical data
    into toroidal coordinates requires this rosetta stone, it can
    be dynamically set to becoming an origin given a specification
    of another origin."""

    _depth = 0
    _origin = []
    sagi = Vec((1.,0.,0.),one)
    meri = Vec((0.,1.,0.),one)
    norm = Vec((0.,0.,1.),one)
    _rot = [sagi.unit,
            meri.unit,
            norm.unit]

    def __init__(self, flag=True):
        
        temp = Vec((0.,0.,0.), one, flag=flag)
        self.unit = temp.unit
        self.s = temp.s
        self.flag = flag
        self.sagi.flag = self.flag
        self.meri.flag = self.flag
        self.norm.flag = self.flag

        # could not use super due to the problem in defining the value of 
        # the depth.  This is simple, though slightly redundant.
        # large number of empty values provide knowledge that there are no
        # lower references or rotations to this, the main coordinate system

def pts2Vec(pt1,pt2):
    """
    pts2Vec creates a vector from pt1 pointing to pt2
    """
    if pt1._origin is pt2._origin:
        return pt2 - pt1
    else:
        raise ValueError("points must exist in same coordinate system")

def fill(funtype,x0,x1,x2,*args,**kwargs):
    try:
        temp = []
        for i in xrange(x0.shape[0]):
            temp+= [fill(funtype,x0[i],x1[i],x2[i],*args,**kwargs)]
        return temp
    except IndexError:
        return funtype((x0,x1,x2),*args,**kwargs)
