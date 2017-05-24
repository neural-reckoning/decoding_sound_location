from shared import *

def identity_compression(input):
    return input

class GainCompression(object):
    def __init__(self, const):
        self.const = const
    def __call__(self, input):
        return input*self.const
    def __str__(self):
        return 'GainCompression_'+'%.2f'%self.const

class PowerLawRectificationCompression(object):
    def __init__(self, const, power):
        self.const = const
        self.power = power
    def __call__(self, input):
        return self.const*clip(input, 0, Inf)**self.power
    def __str__(self):
        return 'PowerLawRectificationCompression_%.2f_%.2f'%(self.const, self.power)
