import warnings
import sys
import scipy

try:
    import mayavi.mlab as mlab
    import tvtk.api
    from tvtk.util.ctf import PiecewiseFunction
except ImportError:
    warnings.warn("mayavi and tvtk modules could not be loaded"
                  "mayavi plotting is unavailable")


def plotLine(vector,val=1.0, close=False, tube_radius=None, index=None, **kwargs):
    """
    PlotLine creates a single plot object from a singular vector or from a n-dimensional
    tuple or list.
    """
    plot = False
    try:
        x = vector.x()
        temp0 = x[0]
        temp1 = x[1]
        temp2 = x[2]
        s =  val*scipy.ones(temp0.shape)

            # For surface objects, this keyword allows for the last corner to connect with the first
        if close:
            temp0 = scipy.concatenate((temp0,scipy.atleast_1d(temp0[0])))
            temp1 = scipy.concatenate((temp1,scipy.atleast_1d(temp1[0])))
            temp2 = scipy.concatenate((temp2,scipy.atleast_1d(temp2[0])))
            s = scipy.concatenate((s,scipy.atleast_1d(s[0])))

        if not index is None:
            N = len(temp0)
            connect = scipy.vstack([scipy.arange(index,   index + N - 1.5),
                                    scipy.arange(index + 1, index + N - .5)]
                                   ).T # I want to rewrite this...
            index += N

    except AttributeError:
             
        temp0 = []
        temp1 = []
        temp2 = []
        s = []
        connect = []
       
        # if it is not some sort of vector or vector-derived class, iterate through and make a surface object 
        if index is None:
            index = 0
            plot = True
                    
        for i in vector:
            output = plotLine(i, close=close, index=index, **kwargs)
            temp0 += [output[0]]
            temp1 += [output[1]]
            temp2 += [output[2]]
            s += [output[3]]
            connect += [output[4]]
            index = output[5]

        #turn to arrays here so I don't accidentally nest lists or tuples
        temp0 = scipy.hstack(temp0)
        temp1 = scipy.hstack(temp1)
        temp2 = scipy.hstack(temp2)
        s = scipy.hstack(s)
        connect = scipy.vstack(connect)

    if index is None:

        try:
            mlab.plot3d(temp0, 
                        temp1,
                        temp2, 
                        s,
                        vmin=0.,
                        vmax=1.,
                        tube_radius=tube_radius,
                        **kwargs)
        except ValueError:
            mlab.plot3d(temp0.flatten(), 
                        temp1.flatten(),
                        temp2.flatten(), 
                        s.flatten(),
                        vmin=0.,
                        vmax=1.,
                        tube_radius=tube_radius,
                        **kwargs)
 
    else:
        if plot:
            # follows http://docs.enthought.com/mayavi/mayavi/auto/example_plotting_many_lines.html#example-plotting-many-lines
   
            src = mlab.pipeline.scalar_scatter(temp0, temp1, temp2, s)
            src.mlab_source.dataset.lines = connect      
            lines = mlab.pipeline.stripper(src)
            mlab.pipeline.surface(lines, **kwargs)
 
        else:
            return (temp0,temp1,temp2,s,connect,index)

def plotTokamak(tokamak, angle=[0,scipy.pi*2], pts=250, section=None, **kwargs):
    temp = []
    for i in xrange(tokamak.norm.s.size):
        temp += [tokamak.getVessel(i).edge(angle=angle, pts=pts)]

    if not section is None:
        
        outline = tokamak.getMachineCrossSection()

        for i in scipy.linspace(angle[0],angle[1],section+1):
            temp += [outline.copy()]
            temp[-1].spin(i)
        plotLine(temp,**kwargs)
    else:
        plotLine(temp[1:],**kwargs)


def plotSurf(surf):

    k = surf.edge().x()
    mlab.mesh(k[0].reshape((2,2)),
              k[1].reshape((2,2)),
              k[2].reshape((2,2)))

    
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
