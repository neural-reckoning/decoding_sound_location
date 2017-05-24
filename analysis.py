from shared import *
from fixedmodel import *
import inspect
#from analyse_model import *
from copy import copy
import multiprocessing
import random as pyrandom

class Analysis(object):
    def __init__(self, analysis_settings, extra_analysis_settings):
        basename = analysis_settings['basename']
        extraname = analysis_settings['extraname']
        analysis_settings_hashable = dict((k, str(v)) for k, v in analysis_settings.items())
        analysis_settings_hash = hash(tuple(sorted(analysis_settings_hashable.items())))
        dataname = generate_extended_datamanager_name(basename, extraname)
        analysisdir = dataname[:-4] + 'analysis_' + str(analysis_settings_hash)
        analysisdir = os.path.normpath(os.path.join(datapath, 'optimal_localisation_decoding/', analysisdir))
        if not os.path.exists(analysisdir):
            os.mkdir(analysisdir)
            os.mkdir(os.path.join(analysisdir, 'data'))
            for k in ['limit_frequency',
                      'use_ideal_responses', 'num_shuffles', 'training_size',
                      'testing_size']:
                k = k + '_' + analysis_settings_hashable[k]
                open(os.path.join(analysisdir, k), 'w')
            f = open(os.path.join(analysisdir, 'filters.txt'), 'w')
            import pprint
            for k in ['training_filters', 'testing_filters', 'sample_filters']:
                f.write(k + ' = ' + pprint.pformat(analysis_settings[k]) + '\n')
            f.close()
        pickle.dump(analysis_settings, open(os.path.join(analysisdir, 'data', 'settings.pkl'), 'wb'), -1)
        self.analysis_settings = self.settings = analysis_settings
        self.moresettings = extra_analysis_settings
        self.analysisdir = os.path.join(analysisdir, 'data')
        self.cached_responses_key = None
        self.cached_responses_fm = None
        self.fm_base = self.get_fixedmodel()
        fm = SubsetFixedModel(self.fm_base)
        fm.use_ideal_responses = analysis_settings['use_ideal_responses']
        limit_frequency = analysis_settings['limit_frequency']
        cf = analysis_settings['cf']
        if limit_frequency is not False:
            print 'CFs limited to range', limit_frequency
            fmin, fmax = limit_frequency
            indices, = logical_and(cf > fmin, cf < fmax).nonzero()
            fm = SubsetFixedModel(fm, indices)
            analysis_settings['Nbin'] = len(indices)
        fm = SubsetFixedModel(fm, analysis_settings['Nbin'])
        self.fm_orig = self.fm = fm
        # Set default Nbinuserange
        cfrepeat = analysis_settings['cfrepeat']
        maxNbin = analysis_settings['Nbin']
        if limit_frequency is not False:
            maxNbin = sum(logical_and(cf > fmin, cf < fmax))
        if cfrepeat == 1:
            #Nbinuserange = unique(hstack((arange(1, maxNbin, 40), [maxNbin])))
            Nbinuserange = unique(hstack((arange(40, maxNbin, 40), [maxNbin])))
        else:
            Nbinuserange = arange(cfrepeat, maxNbin + 1, cfrepeat)
        self.default_Nbinuserange = Nbinuserange
        poolsize = extra_analysis_settings['cpu']
        if poolsize <= 0:
            numprocesses = poolsize + multiprocessing.cpu_count()
        elif poolsize > 0:
            numprocesses = poolsize
        self.poolsize = poolsize

    @property
    def pool(self):
        if not hasattr(self, '_pool'):
            self._pool = multiprocessing.Pool(self.poolsize)
        return self._pool

    def get_fixedmodel(self):
        exec '''
fm = FixedModel(space, acousticnoisemodel, binauralmodel,
                responsenoisemodel, binauraldistribution,
                extended_name,
                gammatone_params=gammatone_params,
                compression=compression,
                rms_normalisation=rms_normalisation,
                use_ram_disk=use_ram_disk,
                )
''' in {'FixedModel':FixedModel,
        'use_ram_disk':self.moresettings['use_ram_disk'],
        }, self.settings
        fm = self.settings.pop('fm')
        self.fm_orig = self.fm = fm
        return fm

    def generate_data(self, stimuli, cpu=0, extra_size=0, forced_size=0, doexit=True):
        training_size = self.settings['training_size']
        testing_size = self.settings['testing_size']
        fm = self.fm_base
        starting_num_saved_responses = fm.num_saved_responses()
        if starting_num_saved_responses<training_size+testing_size+extra_size or forced_size>0:
            remaining_size = training_size + testing_size + extra_size - fm.num_saved_responses()
            if remaining_size < 0: remaining_size = 0
            remaining_size += forced_size
            print 'Already computed', fm.num_saved_responses(), 'responses,', remaining_size, 'remaining'
            fm.parallel_apply_stimuli(stimuli, remaining_size, cpu=cpu)
            if extra_size==training_size==testing_size==0:
                if fm.num_saved_responses()-starting_num_saved_responses<forced_size:
                    print 'Aborted before finishing.'
                    if doexit:
                        exit()
                    return False
            if fm.num_saved_responses() < training_size + testing_size + extra_size:
                print 'Aborted before finishing.'
                if doexit:
                    exit()
                return False
        return True

    def __call__(self, func, estimator, *args, **kwds):
        forcecompute = kwds.pop('forcecompute', False)
        if hasattr(self, 'force_forcecompute'):
            forcecompute = True
        kwds_hashable = dict((k, str(v)) for k, v in kwds.items())
        kwds_hash = hash(tuple(sorted(kwds.items())) + tuple(str(k) for k in args))
        extname = func.__name__ + '_' + str(estimator) + '_' + str(kwds_hash) + '.pkl'
        fname = os.path.join(self.analysisdir, extname)
        found_results = False
        if not forcecompute and os.path.exists(fname):
            try:
                results = pickle.load(open(fname, 'rb'))
                found_results = True
            except:
                print 'Pickling error.'
                found_results = False
        if not found_results:
            print 'Computing for:', extname
            results = func(estimator, kwds_hash, *args, **kwds)
            pickle.dump(results, open(fname, 'wb'))
        return results

    def updateprogress(self):
        pass

    def responses(self, kwd_hash=None):
        if kwd_hash is not None and self.cached_responses_key == kwd_hash and self.cached_responses_fm is self.fm:
            return self.cached_responses
        responses = self.fm.saved_responses()
        random.seed(kwd_hash)
        random.shuffle(responses)
        random.seed()
        self.cached_responses_key = kwd_hash
        self.cached_responses = responses
        self.cached_responses_fm = self.fm
        return responses

    def shuffled_responses(self, hash1, hash2):
        responses = list(self.responses(hash1)) + []
        random.seed(hash2)
        random.shuffle(responses)
        random.seed()
        return responses

    def response_subset(self, responses, filters):
        rejected_responses = []
        for filter in filters:
            newresponses = []
            for i, (R, location, meta, sound) in enumerate(responses):
                localns = meta.copy()
                localns['location'] = location
                localns['i'] = i
                if eval(filter, self.analysis_settings, localns):
                    newresponses.append((R, location, meta, sound))
                else:
                    rejected_responses.append((R, location, meta, sound))
            responses = newresponses
        return responses, rejected_responses

    def traintest_responses(self, kwd_hash=None, shuffle=None):
        if shuffle is not None:
            responses = self.shuffled_responses(kwd_hash, (kwd_hash, shuffle))
        else:
            responses = self.responses(kwd_hash)
        training_responses, rejected = self.response_subset(responses,
                                            self.settings['training_filters'])
        testing_responses, _ = self.response_subset(rejected,
                                            self.settings['testing_filters'])
        return training_responses, testing_responses

    def useful_results(self, results):
        trueval, guessedval, _, _ = zip(*results)
        trueval = array(trueval)
        guessedval = array(guessedval)
        return results, trueval, guessedval

    def traintest_results(self, estimator, kwds_hash, shufflenum=None):
        training_responses, testing_responses = self.traintest_responses(kwds_hash, shufflenum)
        estimator.train(responses=training_responses)
        results = estimator.test(testing_responses)
        results, trueval, guessedval = self.useful_results(results)
        return results, trueval, guessedval, training_responses, testing_responses

    def pre_multiprocessing(self):
        keep = self.updateprogress, self._pool
        self.updateprogress, self._pool = None, None
        return keep

    def post_multiprocessing(self, keep):
        self.updateprogress, self._pool = keep

    ### ANALYSIS RESULTS

    def summary_results(self, estimator, kwds_hash):
        self.updateprogress()
        results, trueval, guessedval, _, _ = self.traintest_results(estimator,
                                                                    kwds_hash)
        return results, trueval, guessedval

    def shuffled_results(self, estimator, kwds_hash, num_shuffles=None):
        if num_shuffles is None:
            num_shuffles = self.settings['num_shuffles']
        shuffled_results = []
        if True:#self.poolsize==1:
            for shufflenum in xrange(num_shuffles):
                self.updateprogress()
                results, trueval, guessedval, _, testing_responses = self.traintest_results(
                                                            estimator, kwds_hash, shufflenum)
                shuffled_results.append((results, trueval, guessedval, testing_responses))
        else:
            pool_args = [(self, estimator, kwds_hash, shufflenum) for shufflenum in xrange(num_shuffles)]
            pool = self.pool
            keep = self.pre_multiprocessing()
            pool_results = pool.map(traintest_results, pool_args, 1)
            self.post_multiprocessing(keep)
            for results, trueval, guessedval, _, testing_responses in pool_results:
                shuffled_results.append((results, trueval, guessedval, testing_responses))
                self.updateprogress()
        return shuffled_results

    def cells_results(self, estimatortype, kwds_hash, num_shuffles=None,
                      Nbinuserange=None, angles=False):
        if num_shuffles is None:
            num_shuffles = self.settings['num_shuffles']
        if Nbinuserange is None:
            Nbinuserange = self.default_Nbinuserange
        meanerror = dict((Nbinuse, []) for Nbinuse in Nbinuserange)
        for i, Nbinuse in enumerate(Nbinuserange):
            for shufflenum in xrange(num_shuffles):
                fm = SubsetFixedModel(self.fm_orig, Nbinuse, self.settings['cfrepeat'])
                self.fm = fm
                estimator = estimatortype(fm)
                self.updateprogress()
                try:
                    results, trueval, guessedval, _, _ = self.traintest_results(
                                                estimator, kwds_hash, shufflenum)
                    if angles:
                        trueval = arcsin(trueval/self.fm.space.itd_max)*180/pi
                        guessedval = clip(guessedval, -self.fm.space.itd_max,
                                         self.fm.space.itd_max)
                        guessedval = arcsin(guessedval/self.fm.space.itd_max)*180/pi
                    me = mean(abs(trueval - guessedval))
                except e:
                    me = nan
                meanerror[Nbinuse].append(me)
        return meanerror

    def cutoff_frequency_results(self, estimatortype, kwds_hash,
                                 num_shuffles=None,
                                 Nbinuserange=None, angles=False):
        if num_shuffles is None:
            num_shuffles = self.settings['num_shuffles']
        if Nbinuserange is None:
            Nbinuserange = self.default_Nbinuserange
        meanerror = dict((Nbinuse, []) for Nbinuse in Nbinuserange)
        cutoff = []
        for i, Nbinuse in enumerate(Nbinuserange):
            for shufflenum in xrange(num_shuffles):
                fm_orig = self.fm_orig
                cf = fm_orig.centre_frequencies[:fm_orig.cfN]
                cfN = fm_orig.cfN
                I = argsort(cf)
                fm = SubsetFixedModel(self.fm_orig, I[:Nbinuse], self.settings['cfrepeat'])
                self.fm = fm
                estimator = estimatortype(fm)
                self.updateprogress()
                try:
                    results, trueval, guessedval, _, _ = self.traintest_results(
                                                estimator, kwds_hash, shufflenum)
                    if angles:
                        trueval = arcsin(trueval/self.fm.space.itd_max)*180/pi
                        guessedval = clip(guessedval, -self.fm.space.itd_max,
                                         self.fm.space.itd_max)
                        guessedval = arcsin(guessedval/self.fm.space.itd_max)*180/pi
                    me = mean(abs(trueval - guessedval))
                except e:
                    me = nan
                meanerror[Nbinuse].append(me)
            cutoff.append(amax(cf[I][:Nbinuse]))
        return meanerror, cutoff
    
    def frequency_results(self, estimatortype, kwds_hash,
                          num_frequency_pools, frequency_pool_width,
                          num_frequency_points,
                          num_shuffles=None,
                          bias_fraction=1,
                          ):
        if num_shuffles is None:
            num_shuffles = self.settings['num_shuffles']
        cfrepeat = self.settings['cfrepeat']
        fm_orig = self.fm_orig
        cf = fm_orig.centre_frequencies[:fm_orig.cfN]
        cfN = fm_orig.cfN
        I = argsort(cf)
        cfsorted = cf[I]
        fmid = []
        allcurI = []
        allcfs = []
        if cfrepeat!=1:
            cfu = sorted(unique(cf))
            for fi in xrange(len(cfu) + 1 - num_frequency_pools):
                curI = I[fi * cfrepeat:(fi + num_frequency_pools) * cfrepeat]
                curcfu = unique(cf[curI])
                fmid.append(mean(curcfu))
                allcurI.append(curI)
                allcfs.append(cf[curI])
        else:
            for istart in array(linspace(0, cfN-frequency_pool_width,
                                         num_frequency_points), dtype=int):
                curI = I[istart:istart+frequency_pool_width]
                fmid.append(mean(cf[curI]))
                allcurI.append(curI)
                allcfs.append(cf[curI])
        error = defaultdict(list)
        bias = defaultdict(list)
        freq_results = defaultdict(dict)
        for fi, curI in enumerate(allcurI):
            for shufflenum in xrange(num_shuffles):
                fm = SubsetFixedModel(fm_orig, curI)
                self.fm = fm
                estimator = estimatortype(fm)
                self.updateprogress()
                results, trueval, guessedval, _, testing_responses = self.traintest_results(
                                            estimator, kwds_hash, shufflenum)
                me = mean(abs(trueval - guessedval))
                if bias_fraction!=1:
                    K = abs(trueval)<(amax(abs(trueval))*bias_fraction)
                    b = 100 * (1 - sum(trueval[K] * guessedval[K]) / sum(trueval[K] ** 2))
                else:
                    b = 100 * (1 - sum(trueval * guessedval) / sum(trueval ** 2))
                error[shufflenum].append(me)
                bias[shufflenum].append(b)
                freq_results[fi][shufflenum] = results, trueval, guessedval, testing_responses
        return allcfs, fmid, error, bias, freq_results

    def artificial_bandpass_results(self, estimatortype, kwds_hash,
                                    num_frequency_pools, frequency_pool_width,
                                    num_frequency_points,
                                    num_shuffles=None):
        if num_shuffles is None:
            num_shuffles = self.settings['num_shuffles']
        if isinstance(frequency_pool_width, int):
            frequency_pool_number = frequency_pool_width
        else:
            frequency_pool_width, frequency_pool_number = frequency_pool_width
        cfrepeat = self.settings['cfrepeat']
        fm_orig = self.fm_orig
        cf = fm_orig.centre_frequencies[:fm_orig.cfN]
        cfN = fm_orig.cfN
        I = argsort(cf)
        cfsorted = cf[I]
        fmid = []
        allcurI = []
        allcfs = []
        if cfrepeat!=1:
            cfu = sorted(unique(cf))
            for fi in xrange(len(cfu) + 1 - num_frequency_pools):
                curI = I[fi * cfrepeat:(fi + num_frequency_pools) * cfrepeat]
                curcfu = unique(cf[curI])
                fmid.append(mean(curcfu))
                allcurI.append(curI)
                allcfs.append(cf[curI])
        else:
            for istart in array(linspace(0, cfN-frequency_pool_width,
                                         num_frequency_points), dtype=int):
                curI = I[istart:istart+frequency_pool_width]
                fmid.append(mean(cf[curI]))
                allcurI.append(curI)
                allcfs.append(cf[curI])
        error = defaultdict(list)
        bias = defaultdict(list)
        freq_results = defaultdict(dict)
        fm = fm_orig
        self.fm = fm
        def zeroify(x, I):
            y = zeros_like(x)
            y[I] = x[I]
            return y
        for shufflenum in xrange(num_shuffles):
            estimator = estimatortype(fm)
            training_responses, testing_responses = self.traintest_responses(kwds_hash, shufflenum)
            estimator.train(responses=training_responses)
            for fi, curI in enumerate(allcurI):
                J = arange(len(curI))
                pyrandom.shuffle(J)
                curI = curI[J[:frequency_pool_number]]
                self.updateprogress()
                bandpassed_testing_responses = [(zeroify(R, curI), location, meta, sound) for (R, location, meta, sound) in testing_responses]
                results = estimator.test(bandpassed_testing_responses)
                results, trueval, guessedval = self.useful_results(results)
                me = mean(abs(trueval - guessedval))
                b = 100 * (1 - sum(trueval * guessedval) / sum(trueval ** 2))
                error[shufflenum].append(me)
                bias[shufflenum].append(b)
                freq_results[fi][shufflenum] = results, trueval, guessedval, testing_responses
        return allcfs, fmid, error, bias, freq_results


def get_analysis_from_namespace():
    import analyse_model
    frame = inspect.stack()[1][0]
    global_ns, local_ns = frame.f_globals, copy(frame.f_locals)
    for name, calc in analyse_model.default_analysis_value_calculations:
        if not name in global_ns and not name in local_ns:
            local_ns[name] = eval(calc, global_ns, local_ns)
    analysis_settings = dict((name, eval(name, global_ns, local_ns)) for name in [
        'basename', 'space', 'Nbin', 'cfrepeat', 'cf', 'bd', 'compression',
        'binauralmodel', 'max_firing_rate', 'rms_normalisation',
        'responsenoisemodel', 'binauraldistribution',
        'limit_frequency', 'use_ideal_responses',
        'num_shuffles', 'training_size', 'testing_size',
        'acousticnoisemodel',
        'training_filters', 'testing_filters', 'sample_filters',
        'extended_name', 'extraname', 'gammatone_params',
        ])
    extra_analysis_settings = dict((name, eval(name, global_ns, local_ns)) for name in [
        'level', 'sn_level_range', 'num_noise_bins',
        'alpha_bins', 'itd_bins', 'use_ram_disk', 'cpu',
        'tonefreq_bins',
        ])
    return Analysis(analysis_settings, extra_analysis_settings)

def traintest_results((analysis, estimator, kwds_hash, shufflenum)):
    return analysis.traintest_results(estimator, kwds_hash, shufflenum)
