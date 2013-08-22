import geometry
import scipy

class Rect(geometry.Origin):

    def __init__(self, x_hat, ref, area, Vec=[], err=scipy.array((0,0,0)), angle=[], flag=[]):

        super(Rect,self).__init__(x_hat, ref, Vec=Vec, err=err, angle=angle, flag=flag)
        self.sagi.s = scipy.atleast_1d(area[0])/2
        self.meri.s = scipy.atleast_1d(area[1])/2 
        # this utilizes an unused attribute of the geometry.Origin where the length
        # of the defining coordinate system unit vectors are used to define the
        # cross sectional area of the rectangle

    def area(self):
        return self.sagi.s*self.meri.s*4

    def edge(self):
        """ return points at the edge of rectangle """        
        p00 = geometry.Point((self.vec - self.meri - self.sagi).x(),self._origin)
        p01 = geometry.Point((self.vec - self.meri + self.sagi).x(),self._origin)
        p10 = geometry.Point((self.vec + self.meri - self.sagi).x(),self._origin)
        p11 = geometry.Point((self.vec + self.meri + self.sagi).x(),self._origin)

        return ((p00,p01),(p10,p11))

    def grid(self):
        """ utilizes geometry.grid to change the rectangle into a generalized surface,
        it is specified with a single set of basis vectors to describe the meridonial,
        normal, and sagittal planes."""
        print('test')

    def plot(self):
        """ return coordinates of the edge values such that the rectangle can be plotted
        """
        temp = self.edge()
        return scipy.array((temp[0][0].x(),temp[0][1].x(),temp[1][1].x(),temp[1][0].x(),temp[0][0].x())).T
"""
class Parabola(Rect):

class Cylinder(Rect):

class Sphere(rect):

class Ellipse(Rect):
"""
