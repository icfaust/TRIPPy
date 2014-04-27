import geometry
import scipy

edges = scipy.array(([-1,-1],
                     [-1, 1],
                     [ 1,-1],
                     [ 1, 1]))

class Surf(geometry.Origin):

    def __init__(self, x_hat, ref, area, vec=None, angle=None, flag=None):
        
        if flag is None:
            flag = ref.flag

        super(Surf,self).__init__(x_hat, ref, vec=vec, angle=angle, flag=flag)
        self.sagi.s = scipy.atleast_1d(area[0])/2
        self.meri.s = scipy.atleast_1d(area[1])/2 
        # this utilizes an unused attribute of the geometry.Origin where the length
        # of the defining coordinate system unit vectors are used to define the
        # cross sectional area of the rectangle



class Rect(Surf):

    def area(self,sagi = None, meri = None):
        if not sagi is None:
            sagi = self.sagi.s
        if not meri is None:
            meri = self.meri.s

        return sagi*meri*4

    def edge(self):
        """ return points at the edge of rectangle """
        temp1 = self.sagi.x()
        temp2 = self.meri.x()
        return geometry.Point((self + geometry.Vecx(scipy.dot(edges,[temp1,temp2]).T)),self._origin)

    def edgetest(self,sagi,meri):

        if abs(meri) <= self.meri and abs(sagi) <= self.sagi:
            return True
        else:
            return False

    def split(self,sagi,meri):
        """ utilizes geometry.grid to change the rectangle into a generalized surface,
        it is specified with a single set of basis vectors to describe the meridonial,
        normal, and sagittal planes."""
        ins = float((sagi-1))/sagi
        inm = float((meri-1))/meri
        stemp = self.sagi.s/sagi
        mtemp = self.meri.s/meri

        self.sagi.s,self.meri.s = scipy.meshgrid(scipy.linspace(-self.sagi.s*ins,self.sagi.s*ins,sagi),
                                                 scipy.linspace(-self.meri.s*inm,self.meri.s*inm,meri))

        x_hat = self + self.sagi + self.meri #creates a vector which includes all the centers of the subsurface
        self.sagi.s = stemp*sagi #returns values to previous numbers
        self.meri.s = mtemp*meri

        temp = Rect(x_hat,
                    self._origin,
                    [2*stemp,2*mtemp],
                    vec=[self.meri.copy(), self.norm.copy()],
                    flag=self.flag)
        #return temp

        return super(Rect, temp).split(temp._origin,
                                       [2*stemp,2*mtemp],
                                       vec=[temp.meri,temp.norm],
                                       flag=temp.flag,
                                       obj=type(temp))
"""
class Parabola(Surf):

class Cylinder(Surf):

class Sphere(Surf):
"""
class Ellipse(Surf):

    def __init__(self, x_hat, ref, area, vec=None, angle=None, flag=None):
 
        super(Surf,self).__init__(x_hat, ref, vec=vec, angle=angle, flag=flag)
        self.sagi.s = scipy.atleast_1d(area[0])
        self.meri.s = scipy.atleast_1d(area[1])


    def edge(self, angle=[0,2*scipy.pi], pts=250):
        theta = scipy.linspace(angle[0], angle[1], pts)
        temp = geometry.Point(geometry.Vecr(scipy.sqrt((self.sagi.s*scipy.cos(theta))**2 +(self.meri.s*scipy.cos(theta))**2)*scipy.ones(theta.shape),
                                            theta,
                                            scipy.zeros(theta.shape))
                              ,self)

        temp.redefine(self._origin)
        return temp

    def area(self, sagi=None, meri=None):
        if not sagi is None:
            sagi = self.sagi.s
        if not meri is None:
            meri = self.meri.s

        return scipy.pi*sagi*meri

    def edgetest(self, meri, sagi):
        if (meri/self.meri)**2+(sagi/self.sagi)**2 <= 1:
            return True
        else:
            return False

class Circle(Ellipse):

    def __init__(self, x_hat, ref, radius, vec=None, angle=None, flag=None):
 
        super(Surf,self).__init__(x_hat, ref, vec=vec, angle=angle, flag=flag)
        self.sagi.s = scipy.atleast_1d(radius)
        self.meri.s = scipy.atleast_1d(radius)

    def area(self, radius=None):
        if not radius is None:
            radius = self.sagi

        super(Circle, self).area(radius, radius)

    def edgetest(self, radius=None):
        if not radius is None:
            radius = self.sagi
        super(Circle, self).edgetest(radius, 0)
