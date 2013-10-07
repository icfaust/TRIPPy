import geometry,scipy,eqtools

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
