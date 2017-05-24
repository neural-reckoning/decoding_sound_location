from shared import *

class AcousticNoiseModel(object):
    '''
    Base class for different types of noise. The generate method should return
    a pair L, R of noises.
    '''
    def generate(self, duration, meta):
        pass

class NoAcousticNoise(object):
    def generate(self, duration, meta):
        return silence(duration, nchannels=2)
    def __str__(self):
        return 'NoAcousticNoise'
    __repr__ = __str__

class IndependentWhiteAcousticNoise(AcousticNoiseModel):
    def __init__(self, level_dB):
        if not isinstance(level_dB, (list, tuple)):
            level_dB = (level_dB, level_dB)
        self.level_dB = level_dB
        
    def generate(self, duration, meta):
        level_min, level_max = self.level_dB
        level = rand()*(level_max-level_min)+level_min
        meta['acoustic_noise_level'] = level
        return whitenoise(duration, nchannels=2).atlevel(level)

    def __str__(self):
        return 'IndependentWhiteAcousticNoise((%d,%d))'%tuple(map(int, self.level_dB))
    __repr__ = __str__
