import scipy
import surface


def pilatus(pt0, ref, vec=None, angle=None, flag=None):

    dim = 1.72e-4
    pix = scipy.array([487,195])
    temp = surface.Rect(pt0,
                        ref,
                        dim*pix,
                        vec=vec,
                        angle=angle,
                        flag=flag)
    return temp.split(pix[0],pix[1])

def pilatus2(pt0, ref, vec=None, angle=None, flag=None):

    dim = 1.72e-4
    pix = scipy.array([487,195])
    temp = surface.Rect(pt0,
                        ref,
                        dim*pix,
                        vec=vec,
                        angle=angle,
                        flag=flag)
    return temp
