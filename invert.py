import plasma
import beam
import scipy


def sens(beams,plasmameth,time,points):
    """ optimal to give multiple times """
    time = scipy.atleast_1d(time)
    # initialize output array of sensitivities
    output = scipy.zeros((len(time),len(beams),len(points)))
    

    for i in xrange(len(time)):
        
