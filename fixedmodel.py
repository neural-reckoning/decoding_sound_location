'''
Putting together the basic pieces to make the model
'''

from shared import *
from stimuli import *
import itertools
from brian.utils.progressreporting import ProgressReporter
from random import shuffle
from brian.tools.taskfarm import run_tasks

# identity_function and gain_function are both just used so that
# we don't have to use lambda functions which aren't picklable
# and so can't be used with multiprocessing or playdoh
def identity_function(x):
    return x

class GainFunction(object):
    def __init__(self, gains):
        self.gains = gains
    def __call__(self, input):
        return self.gains*input

class FixedModel(object):
    '''
    Parameters are:
    
    ``space``
        The Spatial object giving the range of ITDs to use.
    ``acousticnoisemodel``
        The model of acoustical noise to use.
    ``binauralmodel``
        The model of binaural neurons.
    ``responsenoisemodel``
        The model of response noise to use
    ``binauraldistribution``
        The distribution of binaural neurons (cf, itd)
    ``dataman``
        The DataManager used to store results, either a string or DataManager
        object. If it's a string, the path 'optimal_localisation_decoding/' will
        be prepended and a DataManager object created for it. Responses are
        saved with a randomly generated string key in the format
        ``(response, location, meta, sound)``.
    ``store_sounds``
        Set to ``True`` to save the sounds along with the response, location and
        metadata.
        
    **Methods:**
    
    ``isaved_responses()``, ``saved_responses()``
        Iterator/list of saved responses.
    ``apply(sound, location, meta=None, store=False)``
        Applies the model to the given sound at the given location. This
        method is normally used for testing purposes, so the results are
        not stored, but you can do this using store=True, in which case you
        may also want to store some metadata in meta. Returns data in the
        form ``(response, location, meta, sound)``.
    ``iapply_stimuli(stimuli, max_stimuli=None, store=True)``
        Applies the model to a set of stimuli (use at most max_stimuli if
        provided). Returns an iterator of elements of the form
        returned by apply. Note that at the moment it just repeatedly calls
        apply, but this might be vectorised/optimised in the future.
    ``apply_stimuli(stimuli, max_stimuli=None, store=True, return_vals=False)``
        Applies the model to a set of stimuli (use at most max_stimuli if
        provided). If ``return_vals=True`, returns a list of elements of the form
        returned by apply. Note that at the moment it just repeatedly calls
        apply, but this might be vectorised/optimised in the future.
    ``parallel_apply_stimuli(stimuli, max_stimuli, report='stderr', cpu=0)``
        Applies the model to the stimuli using multiple CPUs in parallel.
        You must provide a maximum number of stimuli, and you should note that
        on linux the process forks so that randomisation will produce the same
        results on multiple CPUs
    '''
    def __init__(self, space, acousticnoisemodel, binauralmodel,
                 responsenoisemodel, binauraldistribution,
                 dataman, 
                 compression=None,
                 gains=None, normalisation=None,
                 rms_normalisation=None, rms_normalisation_power=None,
                 store_sounds=False, gammatone_params=None,
                 use_ram_disk=None):
        self.space = space
        self.acousticnoisemodel = acousticnoisemodel
        self.responsenoisemodel = responsenoisemodel
        self.gammatone_params = gammatone_params
        self.binauralmodel = binauralmodel
        self.binauraldistribution = binauraldistribution
        if not isinstance(dataman, DataManager):
            datamanname = 'optimal_localisation_decoding/'+dataman
            if use_ram_disk is None:
                dataman = DataManager(datamanname)
            else:
                dataman = datamanager.DataManager(datamanname, use_ram_disk)
        self.dataman = dataman
        self.store_sounds = store_sounds
        cf = centre_frequencies = binauraldistribution.best_frequencies
        self.cfN = len(cf)
        centre_frequencies = hstack((cf, cf))
        self.centre_frequencies = centre_frequencies
        sound = silence(1*ms, nchannels=2)
        self.soundinput = DoNothingFilterbank(sound)
        # 3. Application of cochlear filtering
        if gammatone_params is None:
            gammatone_params = {}
        self.gammatone_filterbank = Gammatone(
                Repeat(self.soundinput, self.cfN),
                centre_frequencies, **gammatone_params)
        if compression is None:
            compression = identity_function
        self.compression_filterbank = FunctionFilterbank(
                self.gammatone_filterbank,
                compression)
        if gains is None:
            gains = ones(self.cfN*2)
        self.gains = gains
        self.gains_filterbank = FunctionFilterbank(self.compression_filterbank,
                                                   GainFunction(self.gains))
        # 4. Application of BDs
        self.binauraldistribution_filterbank = self.binauraldistribution.filterbank(self.gains_filterbank)
        
        if normalisation is None:
            normalisation = ones(self.cfN)
        self.normalisation = normalisation
        
        self.rms_normalisation = rms_normalisation
        if rms_normalisation_power is None:
            rms_normalisation_power = binauralmodel.power
        self.rms_normalisation_power = rms_normalisation_power
        
        self.use_ideal_responses = False
        
    def isaved_responses(self):
        dataman = self.dataman
        for _, response in self.dataman.iteritems():
            R, location, meta, sound = response
            if self.use_ideal_responses:
                R = meta['noiseless_response']
            if hasattr(self, 'redo_response_noise'):
                R = R*self.normalisation
                R = self.responsenoisemodel(R)
#                print R
                yield (R*1.0, location, meta, sound)
            else:
                yield (R*self.normalisation, location, meta, sound)
    
    def saved_responses(self):
        return list(self.isaved_responses())
    
    def num_saved_responses(self):
        return self.dataman.itemcount()
    
    def normalise_to_saved_responses(self, max_response=1.):
        self.normalisation = ones(self.cfN)
        Rmax = zeros(self.cfN)
        for R, _, _, _ in self.isaved_responses():
            Rmax = maximum(Rmax, R)
        self.normalisation *= max_response/Rmax

    def apply(self, sound, location, meta=None, store=False):
        gc.collect()
        if meta is None:
            meta = {}
        duration = sound.duration
        nsamples = sound.nsamples
        delay_offset = self.binauraldistribution_filterbank.delay_offset
        signal_duration = duration+global_max_delay
        decay_duration = signal_duration+max_decay_time
        bd_end = delay_offset+decay_duration+global_max_delay
        # 1. Application of HRTF
        sound = self.space.apply(sound, location)[:bd_end]
        # 2. Application of noise model
        noise = self.acousticnoisemodel.generate(signal_duration, meta)
        sound = sound+noise
        self.soundinput.source = sound
        # If RMS normalisation is used, we need to do two passes
        if self.rms_normalisation is not None:
            K = self.rms_normalisation
            self.gains[:] = 1
            def sum_of_powers(input, running):
                return running+sum(input**self.rms_normalisation_power, axis=0)            
            rms = self.gains_filterbank.process(sum_of_powers, decay_duration)/int(samplerate*signal_duration)
            rms = rms**(1.0/self.rms_normalisation_power)
            self.gains[:] = K/rms
        # Application of (3, 4, 5)
        R = self.binauralmodel.apply(self.binauraldistribution_filterbank,
                                     sound.duration)
        noiseless_response = R.copy()
        meta['noiseless_response'] = noiseless_response
        R = R*self.normalisation
        R = self.responsenoisemodel(R)
        if not self.store_sounds:
            sound = None
        if store:
            session = self.dataman.computer_session()
            session[self.dataman.make_unique_key()] = (response, location, meta, sound) 
        if self.use_ideal_responses:
            R = meta['noiseless_response']
        return (R, location, meta, sound)

    def iapply_stimuli(self, stimuli, max_stimuli=None, store=True,
                       report=None):
        if isinstance(stimuli, Stimuli):
            stimuli = stimuli()
        if max_stimuli is None:
            report = None
        if isinstance(report, str):
            report = ProgressReporter(report)
        if report is not None: report.start()
        for i, (sound, location, meta) in enumerate(itertools.islice(stimuli, max_stimuli)):
            yield self.apply(sound, location, meta, store)
            if report is not None: report.update((i+1.)/max_stimuli)
        if report is not None: report.finish()

    def apply_stimuli(self, stimuli, max_stimuli=None, store=True,
                      report=None, return_vals=False):
        iter = self.iapply_stimuli(stimuli, max_stimuli, store, report=report)
        if return_vals:
            return list(iter)
        else:
            for _ in iter:
                pass
            
    def parallel_apply_stimuli(self, stimuli, max_stimuli=None, gui=True,
                               cpu=0):
        #if cpu==1:
        #    self.apply_stimuli(stimuli, max_stimuli=max_stimuli, store=True,
        #                       return_vals=False)
        #    return
        seedstart = rand()
        items = itertools.repeat(None, max_stimuli)
        run_tasks(self.dataman, ApplyFixedModelToStimuliTask, items, gui=gui,
                  poolsize=cpu, initargs=(self, stimuli, seedstart),
                  numitems=max_stimuli)
        
    def __str__(self):
        return 'FixedModel()'
    __repr__ = __str__
        

class ApplyFixedModelToStimuliTask(object):
    '''
    We use this for multiple runs of the model in parallel on a single
    machine using run_tasks.
    '''
    def __init__(self, fm, stimuli, seedstart):
        random.seed(seedstart+os.getpid())
        numpy.random.seed(random.randint(0, sys.maxint-1)+1)
        self.fm = fm
        self.stimuli = stimuli()
    
    def __call__(self, dummy):
        return self.fm.apply(*self.stimuli.next())
    

class SubsetFixedModel(FixedModel):
    def __init__(self, fm, indices=None, cfrepeat=1):
        if isinstance(fm, SubsetFixedModel):
            baseindices = fm.indices
        copy_params = ['space', 'acousticnoisemodel', 'responsenoisemodel',
                       'binauralmodel', 'dataman', 'use_ideal_responses',
                       'gammatone_params',
                       'redo_response_noise',
                       ]
        self.original_fixedmodel = fm
        for p in copy_params:
            if hasattr(fm, p):
                setattr(self, p, getattr(fm, p))
        if indices is None:
            indices = arange(fm.cfN)
            num_indices = fm.cfN
        elif isinstance(indices, int):
            num_indices = indices
            indices = arange(fm.cfN)
            if cfrepeat>1:
                repindices = arange(fm.cfN/cfrepeat)
                shuffle(repindices)
                repindices = repindices[:num_indices/cfrepeat]
                repindices = repeat(repindices*cfrepeat, cfrepeat)+tile(arange(cfrepeat), num_indices/cfrepeat)
                indices = indices[repindices]
            else:
                shuffle(indices)
                indices = indices[:num_indices]
        else:
            num_indices = len(indices)
        self.num_indices = num_indices
        if isinstance(fm, SubsetFixedModel):
            self.indices = baseindices[indices]
        else:
            self.indices = indices
        self.binauraldistribution = fm.binauraldistribution._make_subset(indices)
        self.cfN = self.num_indices
        self.centre_frequencies = fm.centre_frequencies[hstack((indices, indices+fm.cfN))]
        self.normalisation = fm.normalisation[indices]
    
    def isaved_responses(self):
        dataman = self.dataman
        for _, response in self.dataman.iteritems():
            R, location, meta, sound = response
            if 'noiseless_response' in meta:
                meta = meta.copy()
                meta['noiseless_response'] = meta['noiseless_response'][self.indices]
            R = R[self.indices]
            if hasattr(self, 'redo_response_noise'):
                R = R*self.normalisation
                R = self.responsenoisemodel(R)
#                print R
                yield (R*1.0, location, meta, sound)
            else:
                yield (R*self.normalisation, location, meta, sound)
            yield (R*self.normalisation, location, meta, sound)    
    
    def _raise_error(self, *args, **kwds):
        raise RuntimeError('Cannot do new computations with SubsetFixedModel')
    
    apply = _raise_error
    iapply_stimuli = _raise_error
    apply_stimuli = _raise_error
    parallel_apply_stimuli = _raise_error
