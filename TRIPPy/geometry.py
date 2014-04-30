import scipy
import warnings
import copy
one = scipy.array([1.0])

def unit(x_hat):
    r"""Checks vector input to be of correct array dimensionality.
    
    Numpy array dimensionality can often create unexpected array sizes
    in multiplications. Unit also forces the input to follow unit 
    vector conventions. It checks the expected cartesian unit 
    vector to be 3xN, and that all elements of x_hat are within the
    range [-1,1]
        
    Args:
        x_hat: Array-like, 3 or 3xN.
            3 dimensional cartesian vector input, which is stores the 
            direction and magnitude as seperate values.  
   
    Returns:
        Vec: Vector object.
        
    Examples:
        Accepts all array like (tuples included) inputs, though all data 
        is stored in numpy arrays.

        Check [1.,0.,0.]::
            
                xdir = unit([1.,0.,0.])

        Check [[1.,2.],[0.,0.],[0.,0]]::
           
                xdir = unit([[1.,2.],[0.,0.],[0.,0]])
               
    Raises:
        ValueError: If any of the dimensions of x_hat are not 3xN or not
            of unit length.
    """

    temp = scipy.squeeze(x_hat)
    if len(temp) != 3 or temp.max() > 1 or temp.min() < -1:
        raise ValueError
    return temp


def Vecx(x_hat, s=None):        
    r"""Generates a cartesian coordinate vector object
        
    Uses the definition:
        
    .. math::
    
        \vec{x}= \texttt{xhat}[0]\hat{x} + \texttt{xhat}[1]\hat{y} + \texttt{xhat}[2]\hat{z}
    
    Capable of storing multiple directions and lengths as a single
    vector, but highly discouraged (from POLS).
        
    Args:
        x_hat: Array-like, 3 or 3xN.
            3 dimensional cartesian vector input, which is stores the 
            direction and magnitude as seperate values.  
 
    Kwargs:
        s: Array-like or scalar float.
            Vector magnitude in meters. When specified, it is assumed that 
            x_hat is representative of direction only and is of unit
            length. Saves in computation as length calculation is avoided.
            
    Returns:
        Vec: Vector object.
        
    Examples:
        Accepts all array like (tuples included) inputs, though all data 
        is stored in numpy arrays.

        Generate an x direction unit vector (xdir)::
            
                xdir = Vecx((1.,0.,0.))

        Generate a cartesian vector (vec1) into direction (1,3,-4)::
            
                vec1 = Vecx(scipy.array([1.,3.,-4.]))

        Generate a cartesian vector (vec2) into direction (2,2,2)::
            
                vec2 = Vecx(scipy.array([1.,1.,1.])/3.0,s=scipy.array(2.0))

        Generate a cartesian vector (vec3) into direction (3,3,3)::
            
                vec3 = Vecx(vec2.unit,s=scipy.array(3.0))
    """
    flag = False
    xin = scipy.array(x_hat, dtype=float)

    if s is None:
        s = scipy.sqrt(scipy.sum(xin**2, axis=0)) 
        #in the case that s is 0, avoid /0 problems
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            xin = scipy.where(s == 0, xin, xin/s)

    return Vec(xin, s, flag=flag)


def Vecr(x_hat, s=None):
    r"""Generates a cylindrical coordinate vector object
        
    Uses the definition:
        
    .. math::
    
        \vec{x}= \texttt{xhat}[0]\hat{r} + \texttt{xhat}[1]\hat{\theta} + \texttt{xhat}[2]\hat{z}
    
    Capable of storing multiple directions and lengths as a single
    vector, but highly discouraged (from POLS).
        
    Args:
        x_hat: Array-like, 3 or 3xN.
            3 dimensional cylindrical vector input, which is stores the 
            direction and magnitude as seperate values.  All values of 
            theta will be aliased to :math:`(\pi,\pi]`
 
    Kwargs:
        s: Array-like or scalar float.
            Vector magnitude in meters. When specified, it is assumed that 
            x_hat is representative of direction only and is of unit
            length. Saves in computation as length calculation avoided.
            
    Returns:
        Vec: Vector object.
        
    Examples:
        Accepts all array like (tuples included) inputs, though all data 
        is stored in numpy arrays.

        Generate an y direction unit vector in cylindrical coords (ydir)::
            
                ydir = Vecr((1.,scipy.pi/2,0.))

        Generate a cartesian vector (vec1) into direction (1,pi/3,-4)::
            
                vec1 = Vecr(scipy.array([1.,scipy.pi/3,-4.]))

        Generate a cartesian vector (vec2) into direction (6,0,8)::
            
                vec2 = Vecr(scipy.array([3.,0.,4.])/5.0,s=scipy.array(10.0))

        Generate a cartesian vector (vec3) into direction (.3,0,.4)::
            
                vec3 = Vecr(vec2.r()/vec2.s,s=scipy.array(.1))
    """

    flag = True    
    xin = scipy.array(x_hat, dtype=float)
        
    if s is None:
        s = scipy.sqrt(x_hat[0]**2 + x_hat[2]**2)
        #in the case that s is 0, avoid /0 problems
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning)
            xin[0] = scipy.where(s == 0, xin[0], xin[0]/s)
            xin[2] = scipy.where(s == 0, xin[2], xin[2]/s)
            
    return Vec((xin[0]*scipy.cos(xin[1]),
                xin[0]*scipy.sin(xin[1]),
                xin[2])
               ,s, flag=flag)       


class Vec(object):
    """Vector object with inherent cartesian backend mathematics.
    
    Creates a new Vec instance which can be set to a default 
    coordinate system of cartesian or cylindrical coordinates.
    All vector mathematics are accomplished in cartesian to 
    simplify computation effort. Cylindrical vectors are
    converted at last step for return.
    
    It is highly recommended to utilize the Vecx and Vecr 
    functions which allow for proper data checks in generating
    vectors.

    Args:
        x_hat: Array-like of size 3 or 3xN in cartesian.
            for all i, x_hat[0][i]**2 + x_hat[1][i]**2 + x_hat[2][i]**2
            is equal to 1.

        s: Array-like of size 1 or shape N.
            Values of the positions of the 2nd
            dimension of f. Must be monotonic without duplicates.

    Kwargs:
        flag: Boolean.
            Sets the default coordinate nature of the vector to 
            cartesian if False, or cylindrical if True.
                
    Examples:   
        Accepts all array like (tuples included) inputs, though
        all data is stored in numpy arrays.

        Generate an x direction unit vector (xdir)::
            
                xdir = Vec((1.,0.,0.),1.0)

        Generate a cartesian vector (vec1) into direction (2,2,2)::
            
                vec1 = Vec(scipy.array([1.,1.,1.])/3.0,scipy.array(2.0))

        Generate a cartesian vector (vec2) into direction (4,4,4)::
            
                vec2 = Vec(vec1.unit,vec1.s*2)
    """

    def __init__(self, x_hat, s, flag=False):
        """ 
        """
        s = scipy.array(s)
        self.unit = unit(x_hat)
        self.s = scipy.squeeze(s)
        self.flag = flag

    def __add__(self, vec):
        """vector addition, x.__add__(y) <==> x+y

        Args:
            vec: Vector object

        Returns:
            Vector object with coordinate system of first argument

        """
        # corrects for some matrix math oddities of numpy
        x0 = self.x0() + vec.x0()
        x1 = self.x1() + vec.x1()
        x2 = self.x2() + vec.x2()
        
        new = Vecx((x0,x1,x2))
        new.flag = self.flag
        return new

    def __sub__(self, vec):
        """vector addition, x.__sub__(y) <==> x-y

        Args:
            vec: Vector object

        Returns:
            Vector object with coordinate system of first argument

        """
        # corrects for some matrix math oddities of numpy
        x0 = self.x0() - vec.x0()
        x1 = self.x1() - vec.x1()
        x2 = self.x2() - vec.x2()
        
        new = Vecx((x0,x1,x2))    
        new.flag = self.flag
        return new

    def __neg__(self):
        """vector negation, x.__neg__() <==> -x

        Returns:
            Vector object with coordinate system of object

        """
        return Vec(-self.unit, self.s, flag=self.flag)

    def __mul__(self, vec):
        """vector dot product, x.__mul__(y) <==> x*y or x.s = x.s*y

        Args:
            val: Vector object or float or numpy array

        Returns:
            Vector object with coordinate system of first argument,
            If the second argument is not a vector object, it 
            modifies the vector magnitude by value vec.    

        """

        try:
            return (self.s*vec.s)*scipy.dot(self.unit.T, vec.unit)
        except AttributeError:
            return Vec(self.unit, vec*self.s, flag=self.flag)

    def __div__(self, val):
        """vector magnitude division, x.__div__(y) <==> x.s = x.s/y

        Args:
            val: float or numpy array
                
        Returns:
            Vector object with coordinate system of object.    

        """

        return Vec(self.unit, self.s/val, flag=self.flag)

    def __getitem__(self, idx):
        """coordinates at index, x.__getitem__(y) <==> x[y]

        Args:
            idx: int or int array
                index of interest.

        Returns:
            numpy array of cartesian or cylindrical coordinate
            values at index.

        """

        if self.flag:
            return self.r()[idx]
        else:
            return self.x()[idx]

    def copy(self):
        """copy of object

        Returns:
            copy of object 
        """


        return copy.deepcopy(self)

    def _cross(self):
        """returns matrix multiplication form of vector cross product

        Returns:
            numpy square 3x3 array of cartesian coordinate
            values at index.  It is assumed that while there
            might be multiple magnitudes to the vector, that
            there are is only a singular direction.

        """

        return  scipy.array(((0,-self.unit[2],self.unit[1]),
                             (self.unit[2],0,-self.unit[0]),
                             (-self.unit[1],self.unit[0],0)))

    def x0(self):
        """returns cartesian coordinate along first dimension

        Returns:
           numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[0]

    def x1(self):        
        """returns cartesian coordinate along second dimension

        Returns:
            numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[1]
    
    def x2(self):
        """returns cartesian coordinate along third dimension

        Returns:
            numpy array of cartesian coordinates in meters

        """
        return self.s*self.unit[2]

    def r0(self):
        """returns cylindrical coordinate along first dimension

        Returns:
           numpy array of cylindrical coordinates in meters

        """
        return self.s*scipy.sqrt(self.unit[0]**2+self.unit[1]**2)

    def r1(self):        
        """returns cylindrical coordinate along second dimension

        Returns:
            numpy array of cylindrical coordinates in radians

        """
        return scipy.arctan2(self.unit[1],self.unit[0])
    
    def r2(self):
        """returns cylindrical coordinate along third dimension

        Returns:
            numpy array of cylindrical coordinates in meters

        """
        return self.x2()

    def t0(self, r, z)
        """returns toroidal distance given cylindrical
        coordinates (r,z).

        Args:
            r: scipy-array of floats or float in meters. r is
                specified in meters.

            z: scipy-array of floats or float in meters. z is
                specified in meters.

        Returns:
            numpy array of cylindrical coordinates in meters

        """
        return scipy.sqrt((self.r0() - r)**2 + (self.x2() - z)**2)


    def t2(self, r, z):
        """returns poloidal angle given cylindrical
        coordinates (r,z)

        Args:
            r: scipy-array of floats or float in meters. r is
                specified in meters.

            z: scipy-array of floats or float in meters. z is
                specified in meters.

        Returns:
            numpy array of cylindrical coordinates in meters

        """
        return scipy.arctan2(self.x2() - z,self.r0() - r)

    def c(self):
        """Conversion of vector to opposite coordinate system

        Returns:
            copy of vector object with opposite coordinate system
            (set with .flag parameter)

        """
        return Vec(self.unit,self.s, flag = (not self.flag))

    def x(self):        
        """return cartesian coordinate values

        Returns:
            numpy array of cartesian coordinates in meters
        """
        return scipy.squeeze([self.x0(),
                              self.x1(),
                              self.x2()])
    
    def r(self):
        """return cylindrical coordinate values

        Returns:
            numpy array of cylindrical coordinates in meters and radians

        """
        return scipy.squeeze([self.r0(),
                              self.r1(),
                              self.x2()])

    def t(self, r, z):
        """return toroidal coordinate values for given cylindrical
        coordinates (r,z).

        Args:
            r: scipy-array of floats or float in meters. r is
                specified in meters.

            z: scipy-array of floats or float in meters. z is
                specified in meters.

        Returns:
            numpy array of cylindrical coordinates in [meters,radians,radians]
            where it is radius in meters, toroidal angle and then poloidal angle.
        """
        return scipy.squeeze([self.t0(),
                              self.r1(),
                              self.t2()])

    def spin(self,angle):
        """Spin vector object about the cylindrical (0,0,1)/norm vector
        axis. This function is different from rot.

        Args:
            angle: Singule element float or scipy array.
        """
        temp = self.r0()/self.s
        self.unit[0] = temp*scipy.cos(angle)
        self.unit[1] = temp*scipy.sin(angle)


    def point(self,ref):
        """returns point based off of vector

        Args:
            ref: reference origin

        Returns:
            Point object

        """
        return Point(self, ref)
                   

def angle(Vec1, Vec2):
    """Returns angle between two vectors.

    Args:
        Vec1: Vec Object

        Vec2: Vec Object

    Returns:
        angle in radians :math:`[0,\pi]` seperating the two vectors
        on the plane defined by the two.
     
    """
    return scipy.arccos((Vec1 * Vec2)/(Vec1.s * Vec2.s)) 

def cross(Vec1, Vec2):
    """Returns angle between two vectors.

    Args:
        Vec1: Vec Object

        Vec2: Vec Object

    Returns:
        Vec object of the vector cross product of the two vectors.
        It is in the coordinate system of the first argument.
     
    """
    new = Vec(scipy.dot(Vec1._cross(),Vec2.unit),Vec1.s*Vec2.s)
    new.flag = Vec1.flag
    return new

class Point(Vec):
    """Point object with inherent cartesian backend mathematics.
     
    Creates a new Point instance which is set to a default 
    coordinate system of cartesian or cylindrical coordinates
    determined from the reference coordinate system. All vector
    mathematics are accomplished in cartesian to simplify 
    computation effort. Cylindrical vectors are converted at
    last step for return.

    Args:
        x_hat: geometry object or 3xN coordinate system values.
             Input coordinate system value options assume that
             it matches the coordinate system of the reference
             origin.

    Kwargs:
        ref: Origin object.
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

    def __init__(self, x_hat, ref=None):        
        """
        """

        try:
            super(Point,self).__init__(x_hat.unit, x_hat.s, flag=x_hat.flag)

        except AttributeError:
            try:
                if ref.flag:
                    x_hat = Vecr(x_hat)
                else:
                    x_hat = Vecx(x_hat)
                    
                super(Point,self).__init__(x_hat.unit, x_hat.s, flag=x_hat.flag)
                    
            except AttributeError:
                raise ValueError('reference not specified')
            
        if ref is None: 
            self._origin = x_hat._origin
            self._depth = x_hat._depth
        else:
            self._origin = ref
            self._depth = ref._depth + 1
         
    def redefine(self, neworigin):        
        """redefine Point object or Point-derived object
        into new coordinate system

        Args:
            neworigin: Origin or Origin-derived object
        """
       
        # use _lca to find common ancestor and return tree to common
        lca = self._lca(neworigin)
        
        # this will allows for the matrix math of the rotation to behave properly
        self._translate(lca, neworigin)

    def _translate(self, lca, neworigin):
        """performs necessary rotations and translations of point.

        Args:
            lca: tuple of Origins or Origin-derived Objects
                Contains a tuple of Origins in order to common
                ancestor of current Object and another tuple of 
                Origins to the new base Origin (neworigin).
             
            neworigin: Origin or Origin-derived object
                New origin for the Point or Point-derived object.

        """

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
        """split coordinate values into seperate objects

        Kwargs:
            obj: geometry-derived object which to form from data

        Returns:
            Object tuple size N for N points. Works to arbitrary
            dimension.

        """
        
        obj = kwargs.pop('obj', None)
        if obj is None:
            obj = type(self)
        temp = self.x()

        if temp.size > 3:
            # initialize
            return fill(obj,temp[0],temp[1],temp[2], *args, **kwargs)

    def _genOriginsToParent(self):
        """Tuple of Origins to Center of space.

        Returns:
            Tuple of Origin or Origin derived objects in order
            of increasing depth
        """

        temp = self._origin
        pnts = self._depth*[0]

        # as index increases, the closer to the point it becomes.
        for idx in range(self._depth-1,-1,-1):
            pnts[idx] = temp
            temp = temp._origin

        return pnts
    
    def _lca(self, point2):
        """Lowest Common Ancestor

        Args:
            point2: Point or Point-derived object

        Returns:
            (pt1,pt2) Tuple of tuples which contain all Origin or
            Origin-derived Objects to lowest common ancestor.
        """
        
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


class Origin(Point):
    """Origin object with inherent cartesian backend mathematics.
     
    Creates a new Origin instance which can be set to a default 
    coordinate system of cartesian or cylindrical coordinates.
    All vector mathematics are accomplished in cartesian to 
    simplify computation effort. Cylindrical vectors are
    converted at last step for return.
    
    An Origin is defined by a point and two vectors. The two 
    vectors being: 1st the normal to the surface, principal axis,
    z-vector or the (0,0,1) vector of the new system defined in
    the reference system. The second vector along with the 
    first fully defines meridonial ray paths, y-vector or the
    (0,1,0) vector (.meri). The sagittal ray path, x-vector or
    the (1,0,0) is defined through a cross product (.sagi).
    Point position and rotation matricies are stored at
    instantiation.

    These conventions of norm to z, meri to y axis and sagi to
    x axis are exactly as perscribed in the OSLO optical code,
    allowing for easier translation from its data into Toroidal
    systems.
           
    If the angles alpha, beta, and gamma are specified following
    the eulerian rotation formalism, it is processed in the 
    following manner: alpha is the rotation from the principal
    axis in the meridonial plane, beta is the rotation about the
    plane normal to the meridonial ray, or 2nd specified vector,
    and gamma is the 2nd rotation about the principal axis. 
    This might change based on what is most physically intuitive.
    These are converted to vectors and stored as attributes.
    
    Args:
        x_hat: geometry-derived object or Array-like of size 3 or 3xN.
            Position in the coordinate system defined by origin, 
            which, if it is a scipy array, follows the input convention
            of the origin. Specified vector will be converted as necessary.

    Kwargs:
        ref: Origin or Origin-derived object.
            Sets the default coordinate nature of the vector to 
            cartesian if False, or cylindrical if True.

        vec: Tuple of two Vector objects
            The two vectors describe the normal (or z) axis and
            the meridonial (or y) axis. Inputs should follow
            [meri,normal]. If not specified, it assumed that angle
            is specified.

        angle: tuple or array of 3 floats
            alpha, beta and gamma are eulerian rotation angles which
            describe the rotation and thus the sagittal and 
            meridonial rays.
            
        flag: Boolean.
            Sets the default coordinate nature of the vector to 
            cartesian if False, or cylindrical if True.
                
    Examples:   
        Accepts all array like (tuples included) inputs, though all data 
        is stored in numpy arrays.

        Generate an origin at (0,0,0) with a :math:`\pi/2` rotation:
            
                cent = Center() #implicitly in cyl. coords.
                newy = Vecr((1.,scipy.pi,0.))
                z = Vecr((0.,0.,1.))
                ex = Origin((0.,0.,0.), cent, vec=[newy,z])

        Generate an origin at (0,0,0) with a :math:`\pi/2` rotation:
            
                cent = Center() #implicitly in cyl. coords.
                ex1 = Origin((0.,0.,0.), cent, angle=(scipy.pi/2,0.,0.))

        Generate an origin at (1,10,-7) with a cartesian coord system:

                cent = Center() #implicitly in cyl. coords.
                place = Vecx((1.,10.,-7.))
                ex2 = Origin(place, cent, angle=(0.,0.,0.), flag=False)

        Generate an origin at (1,1,1) with a cartesian coord system:

                cent = Center(flag=False) #cartesian coords.
                ex3 = Origin((1.,1.,1.), cent, angle=(0.,0.,0.))

    """

    def __init__(self, x_hat, ref=None, vec=None, angle=None, flag=None):
        """ 
        """
        # test Vec1 and Vec2 for ortho-normality
        super(Origin,self).__init__(x_hat, ref=ref)
        if not vec is None:
            # generate point based off of previous origin
            self.norm = vec[1]
            self.meri = vec[0]

            self.sagi = cross(self.meri,self.norm)
            # generate rotation matrix based off coordinate system matching (this could get very interesting)

        elif not angle is None:

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

            self.sagi = Vecx([c3*c2, s3*c1 - c3*s2*s1, c3*s2*c1 + s3*s1], one)
            self.meri = Vecx([-s3*c2, c3*c1 + s3*s2*s1, c3*s1 - s3*s2*c1], one)
            self.norm = Vecx([-s2, -c2*s1, c2*c1], one)

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
        """redefine Origin object or Origin-derived object
        into new coordinate system

        Args:
            neworigin: Origin or Origin-derived object
        """

        lca = self._lca(neworigin)
        super(Origin,self)._translate(lca, neworigin)
        self._rotate(lca, neworigin)


    def _rotate(self, lca, neworigin):
        """performs necessary rotations and translations of Origin.

        Origin or Origin-derived Object requires that the coordinate
        system basis vectors be accurately modified for the new
        origin.  This requires a set of rotations and anti-rotations
        through the lowest common ancestor.

        Args:
            lca: tuple of Origins or Origin-derived Objects
                Contains a tuple of Origins in order to common
                ancestor of current Object and another tuple of 
                Origins to the new base Origin (neworigin).
             
            neworigin: Origin or Origin-derived object
                New origin for the Point or Point-derived object.

        """
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

    def spin(self,angle):
        """Spin vector or vector-derived object around Origin
        about the cylindrical (0,0,1)/norm vector axis. This function
        is different from rot.

        Args:
            vec: Vector or Vector-derived object
        """
        super(Point,self).spin(angle)
        self.sagi.spin(angle)
        self.meri.spin(angle)
        
    def rot(self,vec):
        """Rotate input vector objects into coordinates of Origin.

        Args:
            vec: Vector object
        """
        temp = Vec(scipy.dot(scipy.array(self._rot).T, vec.unit),vec.s)
        temp.flag = vec.flag
        return temp                
    
    def arot(self,vec):
        """Rotate vector out of coordinates of Origin.
        (Anti-Rotate)

        Inverse of rot() function

        Args:
            vec: Vector or Vector-derived object
        """
        temp = Vec(scipy.dot(self._rot, vec.unit), vec.s)      
        temp.flag = vec.flag
        return temp

    def split(self, *args, **kwargs):
        """split coordinate values into seperate objects

        Kwargs:
            obj: geometry-derived object which to form from data.
                If not specified, returns a tuple of Origin objects.

        Returns:
            Object tuple size N for N points. Works to arbitrary
            dimension.

        """

        try:
            obj = kwargs['obj']
            return super(Origin,self).split(*args,**kwargs)
        except KeyError:
            return super(Origin,self).split(self._origin, obj=type(self), Vec=[self.meri,self.norm], flag=self.flag)

        
class Center(Origin):
    """Center object with inherent cartesian backend mathematics.
    
    Creates a new Center instance which can be set to a default 
    coordinate system of cartesian or cylindrical coordinates.
    All vector mathematics are accomplished in cartesian to 
    simplify computation effort. Cylindrical vectors are
    converted at last step for return.  It defaults to cylindrical
    coordinates
    
    The Center class which underlies all positional calculation.
    It is located at (0,0,0) and is inherently a cylindrical coordinate
    system (unless flag set otherwise). It is from the translation
    of inherently cylindrical data into toroidal coordinates 
    requires this rosetta stone, it can be dynamically set to 
    becoming an origin given a specification of another origin.
    
    Kwargs:
        flag: Boolean.
            Sets the default coordinate nature of the vector to 
            cartesian if False, or cylindrical if True.
                
    Examples:   
        Accepts all array like (tuples included) inputs, though all data 
        is stored in numpy arrays.

        Generate a Center in cylindrical coordinates:
            
                cent = Center() #implicitly in cyl. coords.

        Generate a Center in cartesian coordinates:
            
                cent = Center(flag=False) 
    """

    _depth = 0
    _origin = []
    sagi = Vec((1.,0.,0.),one)
    meri = Vec((0.,1.,0.),one)
    norm = Vec((0.,0.,1.),one)
    _rot = [sagi.unit,
            meri.unit,
            norm.unit]

    def __init__(self, flag=True):
        """
        """

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
    """Returns angle between two vectors.

    Args:
        pt1: geometry Object with reference origin

        pt2: geometry Object with reference origin

    Returns:
        Vector object: Vector points from pt1 to pt2.
    
    """
    if pt1._origin is pt2._origin:
        return pt2 - pt1
    else:
        raise ValueError("points must exist in same coordinate system")

def fill(funtype,x0,x1,x2,*args,**kwargs):
    """Recursive function to generate TRIPPy Objects

    Args:
        funtype: Object type to replicate
        
        x0: coordinate of 1st direction

        x1: coordinate of 2md direction
        
        x2: coordinate of 3rd direction

    Returns:
        Tuple or object of type funtype.
    """

    try:
        temp = []
        for i in xrange(x0.shape[0]):
            temp+= [fill(funtype,x0[i],x1[i],x2[i],*args,**kwargs)]
        return temp
    except IndexError:
        return funtype(Vecx((x0,x1,x2)),*args,**kwargs)
