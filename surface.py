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
        if self.sagi.flag:
            temp1 = self.sagi.c().x()
        else:
            temp1 = self.sagi.x()

        if self.meri.flag:
            temp2 = self.meri.c().x()
        else:
            temp2 = self.meri.x()

        return geometry.Point((self.vec + geometry.Vecx(scipy.dot(edges,[temp1,temp2]).T)).x().reshape(3,2,2),self._origin)


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

    def area(self,sagi = [], meri = []):
        if not sagi:
            sagi = self.sagi.s
        if not meri:
            meri = self.meri.s

        return scipy.pi*sagi*meri

class Circle(Ellipse):

    def area(radius):
        super(Circle,self).area(radius,radius)

