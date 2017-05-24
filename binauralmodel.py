'''
Binaural model assumes that the delays have already been applied by the
binaural distribution, and just does the cross-correlation model. The class
takes an input of the form (L1,L2,...R1,R2,...) and returns (B1,B2...)

TODO: accelerated versions using weave/GPU
'''

from shared import *

class BinauralModel(object):
    '''
    Base class for binaural models, handles the routine tasks. Binaural
    models should implement their own apply method, making use of the
    parse method if they wish.
    '''
    def reinit(self):
        pass
        
    def parse(self, x):
        # shape of x is (bufsize, nchannels)
        # channels of x are (L1,L2,...,R1,R2,...)
        L = x[:, :x.shape[1]/2]
        R = x[:, x.shape[1]/2:]
        return L, R

    def apply(self, source, duration, buffersize=32):
        return zeros(source.nchannels)
    
    def __str__(self):
        return 'BinauralModel()'
    __repr__ = __str__
        
class FunctionBinauralModel(BinauralModel):
    '''
    Models of the form response(L, R) = int f(L(t), R(t)) dt
    '''
    def __init__(self, func):
        self.func = func
        self.dt = 1/samplerate
        self.source = None
        self.filterbank = None
    
    def apply(self, source, duration, buffersize=32):
        if self.source is source:
            filterbank = self.filterbank
        else:
            self.source = source
            func = lambda x: self.dt*self.func(*self.parse(x))
            filterbank = FunctionFilterbank(source, func,
                                            nchannels=source.nchannels/2)
            self.filterbank = filterbank 
        if not isinstance(duration, int):
            duration = int(duration*samplerate)
        endpoints = hstack((arange(0, duration, buffersize), duration))
        response = zeros(filterbank.nchannels)
        filterbank.buffer_init()
        for start, end in zip(endpoints[:-1], endpoints[1:]):
            R = filterbank.buffer_fetch(start, end)
            response += sum(R, axis=0)
        return response

    def __str__(self):
        return 'FunctionBinauralModel()'
    __repr__ = __str__

class BinauralPowerFunction(object):
    def __init__(self, a, b, c, power=2, rectify=False):
        self.a, self.b, self.c = a, b, c
        self.power = power
        self.rectify = rectify
    def __call__(self, L, R):
        if self.rectify:
            return self.a*clip(L+R+self.b, 0, Inf)**self.power+self.c
        else:
            return clip(self.a*(L+R+self.b)**self.power+self.c, 0, Inf)
    
class PowerBinauralModel(FunctionBinauralModel):
    '''
    This is the model from Brian Fischer's owl paper
    '''
    def __init__(self, a, b, c, power=2, rectify=False):
        self.a, self.b, self.c = a, b, c
        self.power = power
        self.rectify = rectify
        f = BinauralPowerFunction(a, b, c, power=power, rectify=rectify)
        FunctionBinauralModel.__init__(self, f)
    def __str__(self):
        return 'PowerBinauralModel(%g,%g,%g,%d,%s)'%(self.a, self.b, self.c,
                                                     self.power,
                                                     str(self.rectify))
    __repr__ = __str__

QuadraticBinauralModel = PowerBinauralModel

class BinauralProductFunction(object):
    def __init__(self, rectification=True):
        self.rectification = rectification
    def __call__(self, L, R):
        if self.rectification:
            return clip(L, 0, inf)*clip(R, 0, inf)
        else:
            return L*R
        
class ProductBinauralModel(FunctionBinauralModel):
    '''
    This is just the function f(L,R)=L*R, with L, R half-wave rectified if
    rectification=True (the default).
    '''
    def __init__(self, rectification=True):
        self.rectification = rectification
        f = BinauralProductFunction(rectification=rectification)
        FunctionBinauralModel.__init__(self, f)

    def __str__(self):
        return 'ProductBinauralModel(%s)'%str(self.rectification)
    __repr__ = __str__
