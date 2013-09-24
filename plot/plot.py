import mayavi,scipy,eqtools


class Vessel(object):
    
    
    def __init__(self,plasma):
      self._eq = plasma
        
        
    def genWireFrame(self,pts=250):
      theta = scipy.linspace(0,2*scipy.pi,pts):
      theta = scipy.append(theta,0)
      x = scipy.cos(theta)
      y = scipy.sin(theta)
      z = scipy.ones(theta.shape)
      return x,y,z
      
    def plot(self,pts=250)
      c,s,zee = self.genWireFrame(pts=pts)
      for r,z in self._eq.getMachineCrossSection():
        mayavi.mlab.plot_3d(r*c,r*s,z*zee)
