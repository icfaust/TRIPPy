import mayavi.mlab as mlab,scipy,eqtools
import tvtk.api
import sys
sys.path.append('/home/ian/python/TRIPPy')
import geometry

class Vessel(object):
    
    
    def __init__(self,plasma):
        self._eq = plasma
        
        
    def genWireFrame(self,pts=250):
        theta = scipy.linspace(0,2*scipy.pi,pts)
        theta = scipy.append(theta,0)
        x = scipy.cos(theta)
        y = scipy.sin(theta)
        z = scipy.ones(theta.shape)
        return x,y,z
      
    def plot(self,pts=250,**kwargs):
        c,s,zee = self.genWireFrame(pts=pts)
        r,z = self._eq.getMachineCrossSection()
        for idx in range(r.size):
            mlab.plot3d(r[idx]*c,r[idx]*s,z[idx]*zee,zee,**kwargs)

    def plotVol(self,pts=15,**kwargs):
        fluxGrid = scipy.squeeze(self._eq.getFluxGrid()).T

        datain = scipy.zeros((fluxGrid.shape[0],
                              pts,
                              fluxGrid.shape[1]))

        for idx in range(pts):
            datain[:,idx,:] = fluxGrid



        temp = genCylGrid(self._eq.getRGrid(),
                          scipy.linspace(0,2*scipy.pi,pts),
                          self._eq.getZGrid())
        verticies = genVertsFromPixel(temp)

        hex_type = tvtk.api.tvtk.Hexahedron().cell_type
        temp = temp.reshape((temp.size/3,3))
        verticies = verticies.reshape((verticies.size/8,8))

        sg = tvtk.api.tvtk.UnstructuredGrid(points=temp)
    
        sg.set_cells(hex_type,verticies)
        sg.point_data.scalars = datain.flatten()
        sg.point_data.scalars.name = 'temp'
    
        psi_0 = self._eq.getFluxAxis()
        psi_LCFS = self._eq.getFluxLCFS()
        psi_min = fluxGrid.min()
        psi_max = fluxGrid.max()
        v1 = (psi_0 - psi_min)/(psi_max - psi_min)
        v2 = (psi_LCFS - psi_min)/(psi_max - psi_min)
        mlab.pipeline.volume(sg)

    def plotIso(self,pts=15,**kwargs):
        fluxGrid = scipy.squeeze(self._eq.getFluxGrid()).T

        datain = scipy.zeros((fluxGrid.shape[1],
                              pts,
                              fluxGrid.shape[0]))

        for idx in range(pts):
            datain[:,idx,:] = fluxGrid
            
        temp = genCylGrid(self._eq.getRGrid(),
                          scipy.linspace(0,2*scipy.pi,pts),
                          self._eq.getZGrid())
        temp = temp.reshape((temp.size/3,3))
        
        sgrid = tvtk.api.tvtk.StructuredGrid(dimensions=(self._eq.getRGrid().size,
                                                         pts,
                                                         self._eq.getZGrid().size))
        
        sgrid.points = temp
        sgrid.point_data.scalars = datain.ravel()
        sgrid.point_data.scalars.name = 'scalars'
        mlab.pipeline.iso_surface(sgrid,**kwargs)               

    def plot_sightlines(self,b,color=1,start=.001,end=4,points=100,**kwargs):
        cent = len(b)*[0]
        s = scipy.linspace(start,end,points)
        
        
        for i in range(len(cent)):
            cent[i] = geometry.pts2Vec(b[i],b[-1])
            centerline = geometry.Vecx((-1,0,0))
            cent[i].s = scipy.linspace(start,
                                       scipy.clip(scipy.absolute(b[i][0]/(centerline * cent[i])),
                                                  .001,
                                                  4),
                                       100)
            print(b[i].vec,cent[i])
            cent[i] = b[i].vec + cent[i]
            cent[i] = cent[i].c()
            print(cent[i].s[-1])
            mlab.plot3d(cent[i][0],cent[i][1],cent[i][2],scipy.ones(cent[i][0].shape),**kwargs)






def generate_annulus(r, theta, z):
    """ Generate points for structured grid for a cylindrical annular
    volume.  This method is useful for generating a unstructured
    cylindrical mesh for VTK.
    """
    # Find the x values and y values for each plane.
    x_plane = (scipy.cos(theta)*r[:,None]).ravel()
    y_plane = (scipy.sin(theta)*r[:,None]).ravel()
    
    # Allocate an array for all the points.  We'll have len(x_plane)
    # points on each plane, and we have a plane for each z value, so
    # we need len(x_plane)*len(z) points.
    points = scipy.empty([len(x_plane)*len(z),3])
    
    # Loop through the points for each plane and fill them with the
    # correct x,y,z values.
    start = 0
    for z_plane in z:
        end = start+len(x_plane)
        # slice out a plane of the output points and fill it
        # with the x,y, and z values for this plane.  The x,y
        # values are the same for every plane.  The z value
        # is set to the current z
        plane_points = points[start:end]
        plane_points[:,0] = x_plane
        plane_points[:,1] = y_plane
        plane_points[:,2] = z_plane
        start = end
        
    return points
        
def genCartGrid(x0, x1, x2, edges = False):
    if edges:
        for i in (x0,x1,x2):
            i = scipy.insert(i,0,2*i[1]-i[2])
            i = scipy.append(i,2*i[-1]-i[-2])
            i = (i[1:]+i[:-1])/2
        
    pnts = scipy.empty((x0.size, x1.size, x2.size,3))
    x0in,x1in,x2in = scipy.meshgrid(x0, x1, x2, indexing='ij')
    pnts[:,:,:,0] = x0in
    pnts[:,:,:,1] = x1in
    pnts[:,:,:,2] = x2in
    return pnts

def genCylGrid(x0,x1,x2,edges=False):

    if edges:
        for i in (x0,x1,x2):
            i = scipy.insert(i,0,2*i[1]-i[2])
            i = scipy.append(i,2*i[-1]-i[-2])
            i = (i[1:]+i[:-1])/2

    pnts = scipy.empty((x0.size, x1.size, x2.size,3))
    xin = scipy.dot(scipy.atleast_2d(x0).T, scipy.atleast_2d(scipy.cos(x1)))
    yin = scipy.dot(scipy.atleast_2d(x0).T, scipy.atleast_2d(scipy.sin(x1)))
    zee = scipy.ones(yin.shape)
    for i in range(x2.size):
        pnts[:,:,i,0] = xin
        pnts[:,:,i,1] = yin
        pnts[:,:,i,2] = x2[i]*zee
    return pnts

def genVertsFromPixel(grid):
    """reduces the lengths of the dimensions by 1"""
    output = scipy.empty(scipy.append(scipy.array(grid.shape[:-1]) - 1,8),dtype=int)
    shape = grid.shape
    #INDEXING DEPENDENT
    idx,jdx,kdx = scipy.mgrid[0:shape[0] - 1,
                              0:shape[1] - 1,
                              0:shape[2] - 1]
    output[...,0] = idx + shape[0]*(jdx + shape[1]*kdx)
    output[...,1] = idx + 1 + shape[0]*(jdx + shape[1]*kdx)
    output[...,2] = idx + 1 + shape[0]*(jdx + 1 + shape[1]*kdx)
    output[...,3] = idx + shape[0]*(jdx + 1 + shape[1]*kdx)
    output[...,4] = idx + shape[0]*(jdx + shape[1]*(kdx + 1))
    output[...,5] = idx + 1 + shape[0]*(jdx + shape[1]*(kdx + 1))
    output[...,6] = idx + 1 + shape[0]*(jdx + 1 + shape[1]*(kdx + 1))
    output[...,7] = idx + shape[0]*(jdx + 1 + shape[1]*(kdx + 1))
    return output
