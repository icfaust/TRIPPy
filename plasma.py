import geometry,scipy,eqtools
import scipy.linalg

class Tokamak(geometry.Center):

    def __init__(self,equilib,flag=True):
        self.eq = equilib
        super(Tokamak,self).__init__(flag=flag)

    def genWireFrame(self,pts=250):
        theta = scipy.linspace(0,2*scipy.pi,pts)
        theta = scipy.append(theta,0)
        x = scipy.cos(theta)
        y = scipy.sin(theta)
        z = scipy.ones(theta.shape)
        return x,y,z

    def inVessel(self,ptin):
        Rlim,Zlim = self.eq.getMachineCrossSection()
        return bool(eqtools.inPolygon(Rlim,Zlim,ptin[0],ptin[2]))
            
    def getVessel(self,idx,pts=250):
        Rlim,Zlim = self.eq.getMachineCrossSection()
        if idx < Rlim.size:
            theta = scipy.linspace(0,2*scipy.pi,pts)
            return geometry.Point((Rlim[idx]*scipy.ones(theta.size),
                                   theta,
                                   Zlim[idx]*scipy.ones(theta.size)),self)

    # to be added, inversion functions etc! I'm glad I spent the day on this.

    def pnt2RhoTheta(self, point,t=0, method = 'psinorm', n=0, poloidal_plane=0):
        """ takes r,theta,z and the plasma, and map it to the toroidal position
        this will be replaced by a toroidal point generating system based of
        a seperate vector class"""
        try:
            if not point.flag:
                point = point.c()
        except AttributeError:
            if not point._origin.flag:
                point = point.c()

        r = self.eq.getMagR()
        z = self.eq.getMagZ()

        rho = self.eq.rz2rho(method,point[0],point[2],t)
        theta = (scipy.arctan2(point[2]-z,point[0]-r) + (n*point[1]) + scipy.pi)% (2*scipy.pi) - scipy.pi
        return rho,theta

    def _intersect2d(self,ray,s1,s2):
        Rlim,Zlim = self.eq.getMachineCrossSection()

        l1 = ray(s1)
        l = ray(s2) - ray(s1)
        lin = l.r()

        s = []
        for i in scipy.arange(len(Rlim)-1):
            mat = scipy.array([[lin[0],Rlim[i+1]-Rlim[i]],
                               [lin[2],Zlim[i+1]-Zlim[i]]])

            l2 = l1.r()

            temp = scipy.dot(scipy.linalg.inv(mat),scipy.array([[l2[0]-Rlim[i]],
                                                                [l2[2]-Zlim[i]]]))
            print(temp)

            if temp[1] < 1 and temp[1] > 0:
                s += [temp[0]]

                    
        print(s),
        print('here!!!')
        return scipy.array(s).min()
                                             



    def trace(self,ray,thresh=1e-3,s1=0,s2=None):
        
        if not ray._origin is self:
            ray.redefine(self)

        if s2 is None:
            s2 = ray.tangency(self,sigma=True)

        #if logical_xor(self.inVessel(ray(s2).r()),self.inVessel(ray(s1).r())):
        test = True
        sin = [0,s1,s2]
        op = 1

        while test:
            s = self._intersect2d(ray,sin[1],sin[2])
 
            if s > 1:
                sin[2] = sin[1]+(sin[2]-sin[1])*s
            elif s < 0:
                sin[1] = sin[1]+(sin[2]-sin[1])*s
            else:
                sin[op] = sin[1]+(sin[2]-sin[1])*s
                op = -op

            if abs(sin[2]-sin[0]) < thresh:
                test = False

        return sin[1]

        
