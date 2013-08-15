import geometry
import scipy

class Rect(geometry.Origin):

    def __init__(self, x_hat, ref, area, Vec=[], err=scipy.array((0,0,0)), angle=[], flag=[]):

        self.saggital = scipy.array(area[0])
        self.meridonial = scipy.array(area[1])
        super(Rect,self).__init__(x_hat, ref, Vec=Vec, err=err, angle=angle, flag=flag)

    def area(self):
        return self.saggital*self.meridonial

    def norm(self):
        """ return the vector normal to the surface in the coordinate system
        of the parent origin """
        if self._origin.flag:
            return geometry.Vecr(self._rot[2],s=1)
        else:
            return geometry.Vecx(self._rot[2],s=1)

    def edge(self):
        """ return points at the edge of rectangle """
        x0 = self.saggital/2
        x1 = self.meridonial/2
        
        p00 = geometry.Point(self.x()-x0*self._rot[0]-x1*self._rot[1])
        p01 = geometry.Point(self.x()-x0*self._rot[0]+x1*self._rot[1])
        p10 = geometry.Point(self.x()+x0*self._rot[0]-x1*self._rot[1])
        p11 = geometry.Point(self.x()+x0*self._rot[0]+x1*self._rot[1])

        return ((p00,p01),(p10,p11))


class Parabola(Rect):

class Cylinder(Rect):

class Sphere(rect):

class Ellipse(Rect):
