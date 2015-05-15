import surface, geometry, beam, scipy

def HXR(origin, place = (1.87,-.03,0.), angle=(0.,0.,0.)):
    pos = geometry.Origin(place, origin, angle=angle)
    area = [5e-3,5e-3]
    vec = [geometry.Vecx((0.,0.,1.)), geometry.Vecx((1.,0.,0.))]
    vecy = geometry.Vecx((0.,1.,0.))
    aperture = surface.Rect((0.,0.,0.),pos,area,vec=vec,flag=False)
    output = []

    rdet=scipy.array([[3.908e-1, 0.,1.015e-1],
                      [3.924e-1, 0., 9.51e-2],
                      [3.939e-1, 0., 8.86e-2],
                      [3.953e-1, 0., 8.22e-2],
                      [3.966e-1, 0., 7.57e-2],
                      [3.978e-1, 0., 6.91e-2],
                      [3.988e-1, 0., 6.26e-2],
                      [3.998e-1, 0., 5.61e-2],
                      [4.007e-1, 0., 4.95e-2],
                      [4.014e-1, 0., 4.29e-2],
                      [4.021e-1, 0., 3.63e-2],
                      [4.026e-1, 0., 2.97e-2],
                      [4.030e-1, 0., 2.31e-2],
                      [4.034e-1, 0., 1.65e-2],
                      [4.036e-1, 0., 9.8e-3],
                      [4.037e-1, 0., 3.2e-3],
                      [4.037e-1, 0., -3.4e-3],
                      [4.036e-1, 0., -1.00e-2],
                      [4.034e-1, 0., -1.67e-2],
                      [4.030e-1, 0., -2.33e-2],
                      [4.026e-1, 0., -2.97e-2],
                      [4.021e-1, 0., -3.65e-2],
                      [4.014e-1, 0., -4.31e-2],
                      [4.006e-1, 0., -4.97e-2],
                      [3.997e-1, 0., -5.62e-2],
                      [3.988e-1, 0., -6.28e-2],
                      [3.977e-1, 0., -6.93e-2],
                      [3.965e-1, 0., -7.59e-2],
                      [3.952e-1, 0., -8.23e-2],
                      [3.938e-1, 0., -8.88e-2],
                      [3.923e-1, 0., -9.53e-2],
                      [3.907e-1, 0.,-1.017e-1]])

    for i in rdet:
        vec2 = geometry.Vecx(i)
        veci = [geometry.cross(vec2,vecy),vec2]
        output += [surface.Rect(vec2, pos, area, vec=veci)]
        output[-1].redefine(origin)

    output.append(aperture)
    output[-1].redefine(origin)

    return output

def HXRbeam(alcator):
    temp = HXR(alcator)
    beams = beam.multiBeam(temp[:-1],temp[-1])
    alcator.trace(beams)
    return beams

