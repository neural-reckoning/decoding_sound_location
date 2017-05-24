'''
Sets of stimuli should just be iterable objects of pairs:

  (sound, location, meta)

Where sound is a single channel Sound object, location is an index into
an HRTFSet of locations, and meta is metadata about the stimulus (can be
None if there is no metadata to store). Metadata could include, for example,
the name of the sound file, etc.

Stimuli should in most cases try to include the following meta data:

    time=time.gmtime()
    type: a string identifying the stimuli generation algorithm
    keywords specific to type that give the arguments to the stimuli generation
    duration
    level
    tag: an optional value that can be passed to the stimuli generator
'''

from shared import *

class Stimuli(object):
    '''
    Use this to pass stimuli to FixedModel.parallel_apply_stimuli because
    iterable objects are not picklable.
    
    For example:
    
        stimuli = Stimuli(whitenoise_stimuli, hrtfset, duration=50*ms)
    '''
    def __init__(self, func, *args, **kwds):
        self.func = func
        self.args = args
        self.kwds = kwds
    def __call__(self):
        return self.func(*self.args, **self.kwds)

def whitenoise_stimuli(space, duration=100*ms, level=70*dB, tag=None):
    '''
    Note that level can be a single number or a tuple (min, max) in which case
    levels will be selected uniformly between the two.
    
    The meta information is a dictionary with keys duration, level, tag,
    type='whitenoise', level_range and time=time.gmtime().
    '''
    if isinstance(level, tuple):
        level_range = level
    else:
        level_range = (level, level)
    while True:
        location = space.get_random_itd()
        level = rand()*(level_range[1]-level_range[0])+level_range[0]
        sound = Sound.whitenoise(duration, samplerate).atlevel(level)
        meta = {'level':level, 'tag':tag, 'duration':duration,
                'type':'whitenoise', 'level_range':level_range,
                'time':time.gmtime()}
        yield sound, location, meta

def powerlawnoise_stimuli(space, duration=100*ms, level=70*dB, alpha=(0, 2),
                          tag=None):
    '''
    Note that level and alpha can be a single number or a tuple (min, max) in which case
    levels/alphas will be selected uniformly between the two.
    
    The meta information is a dictionary with keys duration, alpha, level, tag,
    type='powerlawnoise', level_range, alpha_range and
    time=time.gmtime().
    '''
    if isinstance(level, tuple):
        level_range = level
    else:
        level_range = (level, level)
    if isinstance(alpha, tuple):
        alpha_range = alpha
    else:
        alpha_range = (alpha, alpha)
    while True:
        location = space.get_random_itd()
        level = rand()*(level_range[1]-level_range[0])+level_range[0]
        alpha = rand()*(alpha_range[1]-alpha_range[0])+alpha_range[0]
        sound = powerlawnoise(duration, alpha, samplerate).atlevel(level)
        meta = {'level':level, 'tag':tag, 'duration':duration, 'alpha':alpha,
                'type':'powerlawnoise', 'level_range':level_range,
                'alpha_range':alpha_range,
                'time':time.gmtime()}
        yield sound, location, meta

class WhiteNoiseStimuli(Stimuli):
    __doc__ = whitenoise_stimuli.__doc__
    def __init__(self, *args, **kwds):
        Stimuli.__init__(self, whitenoise_stimuli, *args, **kwds)

class PowerLawNoiseStimuli(Stimuli):
    __doc__ = powerlawnoise_stimuli.__doc__
    def __init__(self, *args, **kwds):
        Stimuli.__init__(self, powerlawnoise_stimuli, *args, **kwds)

natural_sound_directories = [
    #'../../sounds/pittsburgh_natural_sounds/*.wav',
    '../../sounds/neotropical_rainforest_mammals/CD 1/*.wav',
    '../../sounds/neotropical_rainforest_mammals/CD 2/*.wav',
    ]

def natural_sounds_stimuli(space, duration=100*ms, level=70*dB,
                           noiselevel=None,
                           tag=None,
                           directories=natural_sound_directories):
    '''
    Note that level can be a single number or a tuple (min, max) in which case
    levels will be selected uniformly between the two.
    
    The meta information is a dictionary with keys duration, level, tag,
    type='natural', filename, channel, start, level_range and time=time.gmtime().
    '''
    if isinstance(level, tuple):
        level_range = level
    else:
        level_range = (level, level)
    files = []
    for directory in natural_sound_directories:
        directory = os.path.join(os.path.split(__file__)[0], directory)
        files.extend(glob.glob(directory))
    while True:
        file = files[randint(len(files))]
        sound = loadsound(file)
        if sound.duration<duration:
            continue
        if int(sound.samplerate)!=int(samplerate):
            continue
        start = rand()*(sound.duration-duration)
        channel = randint(sound.nchannels)
        sound = sound.channel(channel)
        sound = sound[start:start+duration]
        sound = Sound(array(asarray(sound), copy=True))
        gc.collect()
        location = space.get_random_itd()
        level = rand()*(level_range[1]-level_range[0])+level_range[0]
        try:
            sound = sound.atlevel(level)
        except OverflowError:
            continue
        if noiselevel is not None:
            sound = sound+whitenoise(sound.duration).atlevel(noiselevel)
        meta = {'level':level, 'tag':tag, 'duration':duration,
                'noiselevel':noiselevel,
                'type':'natural', 'filename':file, 'start':start,
                'channel':'channel',
                'level_range':level_range,
                'time':time.gmtime()}
        yield sound, location, meta

class NaturalSoundsStimuli(Stimuli):
    __doc__ = natural_sounds_stimuli.__doc__
    def __init__(self, *args, **kwds):
        Stimuli.__init__(self, natural_sounds_stimuli, *args, **kwds)

def tone_stimuli(space, freqrange, duration=100*ms, level=70*dB,
                 tag=None):
    '''
    Note that level can be a single number or a tuple (min, max) in which case
    levels will be selected uniformly between the two.
    
    freqrange should be a sequence of at least two elements. Tone frequencies
    are chosen by first randomly selecting an interval defined by two successive
    elements on freqrange, and then uniformly randomly selecting a frequency
    in this interval. By passing the sorted cf array, you will get tones
    approximately distributed as the cfs.
    
    The meta information is a dictionary with keys duration, level, tonefreq, tag,
    type='tone', level_range, freqrange and
    time=time.gmtime().
    '''
    if isinstance(level, tuple):
        level_range = level
    else:
        level_range = (level, level)
    freqrange = array(freqrange)
    while True:
        location = space.get_random_itd()
        level = rand()*(level_range[1]-level_range[0])+level_range[0]
        i = randint(len(freqrange)-1)
        flow, fhigh = freqrange[i:i+2]
        tonefreq = (rand()*(fhigh-flow)+flow)*Hz
        sound = tone(tonefreq, duration, samplerate).atlevel(level)
        meta = {'level':level, 'tag':tag, 'duration':duration,
                'tonefreq':tonefreq,
                'type':'tone', 'level_range':level_range,
                'freqrange':freqrange,
                'time':time.gmtime()}
        yield sound, location, meta

class ToneStimuli(Stimuli):
    __doc__ = tone_stimuli.__doc__
    def __init__(self, *args, **kwds):
        Stimuli.__init__(self, tone_stimuli, *args, **kwds)


def bandpassnoise_stimuli(space, freqrange, duration=100*ms, level=70*dB,
                 tag=None):
    '''
    Note that level can be a single number or a tuple (min, max) in which case
    levels will be selected uniformly between the two.
    
    freqrange should be a sequence of at least two elements. Centre frequencies
    are chosen by first randomly selecting an interval defined by two successive
    elements on freqrange, and then uniformly randomly selecting a frequency
    in this interval. By passing the sorted cf array, you will get bandpass CFs
    approximately distributed as the cfs.
    
    The meta information is a dictionary with keys duration, level, bpcfreq, tag,
    type='bandpassnoise', level_range, freqrange and
    time=time.gmtime().
    '''
    if isinstance(level, tuple):
        level_range = level
    else:
        level_range = (level, level)
    freqrange = array(freqrange)
    while True:
        location = space.get_random_itd()
        level = rand()*(level_range[1]-level_range[0])+level_range[0]
        i = randint(len(freqrange)-1)
        flow, fhigh = freqrange[i:i+2]
        tonefreq = (rand()*(fhigh-flow)+flow)*Hz
        sound = Sound(Gammatone(whitenoise(duration, samplerate),
                                [tonefreq]).process(), samplerate).atlevel(level)
        meta = {'level':level, 'tag':tag, 'duration':duration,
                'bpcfreq':tonefreq,
                'type':'bandpassnoise', 'level_range':level_range,
                'freqrange':freqrange,
                'time':time.gmtime()}
        yield sound, location, meta

class BandPassNoiseStimuli(Stimuli):
    __doc__ = tone_stimuli.__doc__
    def __init__(self, *args, **kwds):
        Stimuli.__init__(self, bandpassnoise_stimuli, *args, **kwds)

if __name__=='__main__':
    from spatial import *
    space = Spatial(itd_max=100*usecond)
    cf = erbspace(100*Hz, 1500*Hz, 480)
    f = []
    for sound, location, meta in bandpassnoise_stimuli(space, cf):
#        f.append(meta['tonefreq'])
#        if len(f)>=10000:
#            break
        for k, v in meta.iteritems():
            print k, ':', v
        plot(sound)
        show()
        exit()
    hist(f)
    show()
    