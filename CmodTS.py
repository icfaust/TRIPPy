import beam,geometry,scipy

def Chords():
    output = 16*[]
    r1 = scipy.array(16*[.69])
    r2 = scipy.array(16*[.8565])
    z1 = scipy.array([.166,.044,-.027,-.197,-.082,-.139,.108,-.254,.196,.138,0.,.082,-.054,-.112,.-.168,-.227])
    z2 = 1e-3*scipy.array([12.5197,
                           3.3228,
                           -2.0386,
                           -14.8510,
                           -6.1898,
                           -10.4870,
                           8.1509,
                           -19.1277,
                           14.7758,
                           10.4119,
                           0.,
                           6.1901,
                           -4.0768,
                           -8.4525,
                           -12.6703,
                           -17.1035])
    for i in xrange(16):
        output[i] = beam.Ray(geometry.Point((r1[i],0.,z1[i])),
                             geometry.Point((r2[i],0.,z2[i])))
    return output
    
