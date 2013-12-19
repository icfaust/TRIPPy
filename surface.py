import geometry
import scipy

edges = scipy.array(([-1,-1],
                     [-1, 1],
                     [ 1,-1],
                     [ 1, 1]))

class Surf(geometry.Origin):

    def __init__(self, x_hat, ref, area, Vec=[], err=scipy.array((0,0,0)), angle=[], flag=[]):

        super(Surf,self).__init__(x_hat, ref, Vec=Vec, err=err, angle=angle, flag=flag)
        self.sagi.s = scipy.atleast_1d(area[0])/2
        self.meri.s = scipy.atleast_1d(area[1])/2 
        # this utilizes an unused attribute of the geometry.Origin where the length
        # of the defining coordinate system unit vectors are used to define the
        # cross sectional area of the rectangle



class Rect(Surf):

    def area(self,sagi = [], meri = []):
        if not sagi:
            sagi = self.sagi.s
        if not meri:
            meri = self.meri.s

        return sagi*meri*4

    def edge(self):
        """ return points at the edge of rectangle """
        temp1 = self.sagi.x()
        temp2 = self.meri.x()
        return geometry.Point((self + geometry.Vecx(scipy.dot(edges,[temp1,temp2]).T)),self._origin)

    def edgetest(self,meri,sagi):

        if abs(meri) <= self.meri and abs(sagi) <= self.sagi:
            return True
        else:
            return False

    def grid(self):
        """ utilizes geometry.grid to change the rectangle into a generalized surface,
        it is specified with a single set of basis vectors to describe the meridonial,
        normal, and sagittal planes."""
        print('test')

"""
class Parabola(Surf):

class Cylinder(Surf):

class Sphere(Surf):
"""
class Ellipse(Surf):

    def area(self,sagi = None, meri = None):
        if not sagi is None:
            sagi = self.sagi.s
        if not meri is None:
            meri = self.meri.s

        return scipy.pi*sagi*meri

    def edgetest(self,meri,sagi):
        if (meri/self.meri)**2+(sagi/self.sagi)**2 <= 1:
            return True
        else:
            return False


class Circle(Ellipse):

    def area(self,radius = None):
        if not radius is None:
            radius = self.sagi

        super(Circle,self).area(radius,radius)

    def edgetest(self,radius = None):
        if not radius is None:
            radius = self.sagi
        super(Circle,self).edgetest(radius,0)
