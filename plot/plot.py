import mayavi.mlab as mlab,scipy,eqtools
import tvtk.api
from tvtk.util.ctf import PiecewiseFunction
import sys
sys.path.append('/home/ian/python/TRIPPy')
import geometry


def plotLine(vector,close = False, **kwargs):
    temp = vector.x()
    temp = temp.reshape(3,temp.size/3)
    print(temp)
    if close:
        temp = scipy.concatenate((temp,
                                  scipy.atleast_2d(temp[:,0]).T),
                                 axis = 1)

    mlab.plot3d(temp[0],temp[1],temp[2],scipy.ones((temp.size/3,)),**kwargs)
    
def plotView(rays,pts=None, **kwargs):
      

           
    if not pts is None:
        x = scipy.zeros((len(rays)+1,pts))
        y = scipy.zeros(x.shape)
        z = scipy.zeros(x.shape)
        for i in rays:
            i.norm.s = scipy.linspace(i.norm.s[0],i.norm.s[-1],pts)
    else:
        x = scipy.zeros((len(rays)+1,len(rays[0].norm.s)))
        y = scipy.zeros(x.shape)
        z = scipy.zeros(x.shape)
    

    for i in xrange(len(rays)):
        if rays[i]._origin.flag:
            rays[i] = rays[i].c()
        x[i] = rays[i].x()[0]
        y[i] = rays[i].x()[1]
        z[i] = rays[i].x()[2]

    x[-1] = rays[0].x()[0]
    y[-1] = rays[0].x()[1]
    z[-1] = rays[0].x()[2]

    mlab.mesh(x,y,z,**kwargs)

def plotVol(volume,pts=15,**kwargs):
    fluxGrid = scipy.squeeze(volume.getFluxGrid()).T

    datain = scipy.zeros((fluxGrid.shape[0],
                          pts,
                          fluxGrid.shape[1]))

    for idx in range(pts):
        datain[:,idx,:] = fluxGrid

        temp = genCylGrid(plasma.eq.getRGrid(),
                      scipy.linspace(0,2*scipy.pi,pts),
                      plasma.eq.getZGrid())
    verticies = genVertsFromPixel(temp)

    hex_type = tvtk.api.tvtk.Hexahedron().cell_type
    temp = temp.reshape((temp.size/3,3))
    verticies = verticies.reshape((verticies.size/8,8))

    sg = tvtk.api.tvtk.UnstructuredGrid(points=temp)
    
    sg.set_cells(hex_type,verticies)
    sg.point_data.scalars = datain.flatten()
    sg.point_data.scalars.name = 'temp'
    
    psi_0 = plasma.eq.getFluxAxis()
    psi_LCFS = plasma.eq.getFluxLCFS()
    psi_min = fluxGrid.min()
    psi_max = fluxGrid.max()
    v1 = (psi_0 - psi_min)/(psi_max - psi_min)
    v2 = (psi_LCFS - psi_min)/(psi_max - psi_min)
    mlab.pipeline.volume(sg)


def plotSymIso(plasma,pts=15,**kwargs):
    fluxGrid = scipy.squeeze(plasma.eq.getFluxGrid()).T
    
    datain = scipy.zeros((fluxGrid.shape[1],
                          pts,
                          fluxGrid.shape[0]))
    
    for idx in range(pts):
        datain[:,idx,:] = fluxGrid
        
    temp = genCylGrid(plasma.eq.getRGrid(),
                      scipy.linspace(0,2*scipy.pi,pts),
                      plasma.eq.getZGrid())
    temp = temp.reshape((temp.size/3,3))
        
    sgrid = tvtk.api.tvtk.StructuredGrid(dimensions=(plasma.eq.getRGrid().size,
                                                     pts,
                                                     plasma.eq.getZGrid().size))
        
    sgrid.points = temp
    sgrid.point_data.scalars = datain.ravel()
    sgrid.point_data.scalars.name = 'scalars'
    mlab.pipeline.iso_surface(sgrid,**kwargs)   


def plotVol2(plasma,pts=None,lim=False,**kwargs):

    if pts is None:
        pts = plasma.eq.getRGrid().size

    zmin = plasma.eq.getMachineCrossSection()[1].min()
    zmax = plasma.eq.getMachineCrossSection()[1].max()
    rmax = plasma.eq.getMachineCrossSection()[0].max()
    x,y,z = scipy.mgrid[-rmax:rmax:(2*rmax)/pts,-rmax:rmax:(2*rmax)/pts,zmin:zmax:(zmax-zmin)/pts]
    r = (x**2 + y**2)**.5
    psi = plasma.eq.rz2psi(r,z)
    if lim:
        vmin = plasma.eq.getFluxAxis()*-1
        vmax = plasma.eq.getFluxLCFS()
        mmax = plasma.eq.getFluxGrid().max()
        vmin = 0
        vmax = (vmax - vmin)/(mmax-vmin)
        print(vmax)
        vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(x,y,z,psi),vmin=0,vmax=vmax,**kwargs)
#        otf = PiecewiseFunction()
#        otf.add_point(-vmin, .8)
#        vol._otf = otf
#        vol._volume_property.set_scalar_opacity(otf)
    else:
        mlab.pipeline.volume(mlab.pipeline.scalar_field(x,y,z,psi),**kwargs)

        
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

def plotTangency(plasma,beam):
    """ nonlinear minimization of r^2 between the plasma center (in R,Z) versus defined beam s
    vector """
    print('not implemented yet')
