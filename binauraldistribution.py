'''
Binaural distribution needs the following characteristics:

+ Set of triples (best_frequency, delay_left, delay_right)
'''

from shared import *

def delays_from_bd(bd):
    dl = -bd/2
    dr = bd/2
    return dl, dr

class BinauralDistribution(object):
    '''
    Base class for distributions of binaural properties, including best
    frequency and ITD. All derived classes must have the best_frequencies,
    delay_left and delay_right attributes, and should probably just consist
    of an __init__ method that calls this one.
    '''
    def __init__(self, params, cfrepeat=1):
        bf, dl, dr = zip(*params)
        self.best_frequencies = bf
        self.delay_left = array(dl)
        self.delay_right = array(dr)
        self.best_delay = array(dr)-array(dl)
        self.cfrepeat = cfrepeat
        
    def _make_subset(self, indices, cfrepeat=None):
        if cfrepeat is None:
            cfrepeat = self.cfrepeat
        bf = array(self.best_frequencies)[indices]
        dl = array(self.delay_left)[indices]
        dr = array(self.delay_right)[indices]
        return BinauralDistribution(zip(bf, dl, dr), cfrepeat=cfrepeat)
    
    def filterbank(self, source):
        nchannels = source.nchannels
        delays = hstack((self.delay_left, self.delay_right))
        return FractionalDelay(source, delays)
    
    def __str__(self):
        props = [self.best_frequencies, self.delay_left, self.delay_right,
                 self.cfrepeat]
        strhash = str(hash(tuple(map(str, props))))
        return '%s%s'%(self.__class__.__name__, strhash)
