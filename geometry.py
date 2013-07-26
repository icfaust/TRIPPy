import scipy

class Hat(Object):
    """ explicitly just unit vector without error
    which is then defined for the various coordinate
    systems based on classes based on this, a Hat class
    never interacts with another Hat class, only one 
    with a defined coordinate system"""
    def __init__(self, x_hat):

        self.unit = scipy.matrix(x_hat)
        self.flag = []

    def __mul__(self,hat):
        """ Dot product """
        if self.flag == hat.flag:
            return self.unit.T*hat.unit
        else:
            return self.unit.T*hat.c().unit

    def _cross(self):
        " matrix necessary for a cross-product calculation"""
        return  scipy.matrix(((0,-self.unit[2],self.unit[1]),
                              (self.unit[2],0,-self.unit[0]),
                              (-self.unit[1],self.unit[0],0)))


class Vecx(Hat):
    """ explicitly a cartesian unit vector, but can be set as 
    a cylindrical unit vector by setting the flag to true, all
    vector math defaults to first vector"""

    def __init__(self, x_hat,r=[]):
        #if r is specified, it is assumed that x_hat has unit length
        if not r:
            r = scipy.sqrt(scipy.sum(x_hat**2))
            x_hat /= r
        super(CartHat,self).__init__(x_hat)
        self.flag = False
        self.r = r
        
    def c(self):
        """ convert to cylindrical coord """
        return Vecr((scipy.sqrt(self.unit[0]**2+self.unit[1]**2),
                     scipy.arctan2(self.unit[1],self.unit[0]),
                     scipy.unit[2]),
                    r=r)
                       

class Vecr(Hat):
    """ explicitly a cylindrical unit vector, but can be set as 
    a cartesian unit vector by using the c call, all
    vector math defaults to first vector"""
    
    def __init__(self, x_hat, r=[]):
        if not r:
            r = scipy.sqrt(x_hat[0]**2 + x_hat[2]**2)
            x_hat /= r
        super(CylHat,self).__init__(x_hat) # not correct, this will modify the angular variable.
        self.flag = True
        self.r = r
        
    def c(self):
        """ convert to cartesian coord """
        return Vecx((self.unit[0]*scipy.cos(self.unit[1]),
                    self.unit[0]*scipy.sin(self.unit[1]),
                    self.unit[2]),
                    r=r)

def angle(Vec1,Vec2):
    return scipy.arccos(Vec1.hat * Vec2.hat) 

def cross(Vec1,Vec2):
    if Vec1.flag == Vec2.flag:
        return (Vec1.r*Vec2.r)*(Vec1._cross() * Vec2)      
    else:
        return (Vec1.r*Vec2.r)*(Vec1._cross() * Vec2.c())
            

class Point(Vector):
    """ a point class can only be defined relative to an origin"""
    def __init__(self, x_hat, ref, err=scipy.array((0,0,0))):
        
        self.loc = scipy.array(x_hat)
        self.error = scipy.array(err)
        self._origin = ref
        self._depth = ref._depth + 1

    def dist(self, origin):
        common = self._lca(origin)
        for i in 

    def err(self,origin):
        

    def r(self, origin):
        return scipy.sqrt(scipy.sum(self.dist(origin)**2))

    def r_err(self, origin):
        return scipy.sqrt(scipy.sum((self.dist(origin)*self.err(origin))**2)/self.r(origin))
        
    def redefine(self, neworigin):
        """ changes depth of point by calculating with respect to a new
        origin, for calculations with respect to a flux grid, etc. this
        should reduce caluclation substantially."""
        
        self.loc += self.dist(neworigin)
        self.error += self.err(neworigin)
        
        self._origin = neworigin
        self._depth = neworigin._depth + 1

    def _genPointsToParent(self, depth=self._depth):
        """ generate a list of points which leads to the overall basis
        origin of the geometry, its length will be depth"""
        temp = self
        pnts = depth*[0]
        
        for idx in range(depth):
            temp = temp._origin
            pnts[idx] = temp

        return pnts

    def _lca(self, point2):
        """ recursively solve for the common point of reference between two points
        as given by a depth number starting from the base of the tree. It will return
        a number which then provides the looping for various computations."""
        
        pt1 = self._getPointsToParent()
        pt2 = point2._getPointsToParent()
        lim = scipy.min((point2._depth,self._depth))
        found = []

        # find first uncommon ancestor
        idx = -1
        while not found:       
            if pt1[idx] is pt2[idx]:
                idx += -1
                if not lim + idx:
                    found = lim
            else:
                found = -1 - idx
        
        return found
