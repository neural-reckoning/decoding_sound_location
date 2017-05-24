from shared import *

# Default values
num_noise_bins = 20
itdmax_extend = 1.1
alpha_bins = linspace(0, 2, 10)
limit_frequency = False
level = 70*dB
sn_level_range = (-15*dB, 30*dB)
cpu = 1
default_analysis_value_calculations = (
    ('itd_bins', 'linspace(-space.itd_max*itdmax_extend, space.itd_max*itdmax_extend, 40)'),
    ('sample_filters', 'training_filters'),
    ('extended_name', 'generate_extended_datamanager_name(basename, extraname)'),
    ('tonefreq_bins', 'linspace(amin(cf), amax(cf), 10)'),
    )
use_ram_disk = None
curtask = 0

def updateprogress():
    global curtask
    progress.equal_subtask(curtask, numtasks)
    progress.update(0.0)
    curtask += 1

def do_show_sample_results():
    ############ SHOW SAMPLE RESULTS #######################################
    if show_responses>0:
        responses = fm.saved_responses()
        random.shuffle(responses)
        for R, location, _, _ in analysis.response_subset(responses, sample_filters)[0][:show_responses]:
            figure()
            suptitle(basename)
            fmbd = fm.binauraldistribution.best_delay
            plot(fmbd, R, '.')
            axvline(float(location))
            title(str(location))

def do_show_summary_info():
    ############## SHOW SUMMARY INFO ########################################
    if show_summary_info:
        estimators = [(f(fm), name) for f, name in estimator_types]
        for estimator, name in estimators:
            if __name__=='__main__':
                updateprogress()
            results, trueval, guessedval = analysis(analysis.summary_results,
                                                    estimator,
                                                    forcecompute=force_summary_compute)
            print 'Number of test results:', len(results)
            print name, 'mean error (usecond):', mean(abs(trueval-guessedval))*second/usecond
            
    # Don't show any more sample figures after this point
    no_more_sample_figs()

def do_show_analysis():
    ######### SHOW PERFORMANCE AS A FUNCTION OF ACOUSTIC NOISE LEVEL #########
    if performance_with_acoustic_noise:
        fig_performance_with_acoustic_noise(analysis, estimator_types)

    ######## SHOW PERFORMANCE AND BIAS AS A RESULT OF POWER LAW NOISE ALPHA ###
    if performance_with_stimuli_alpha:
        fig_performance_with_stimuli_alpha(analysis, estimator_types)
            
    ################## SHOW CONFUSION MATRICES #############################
    if show_error_matrix:
        show_confusion_matrices(analysis, estimator_types)

    ########### SHOW PERFORMANCE AND BIAS AS A FUNCTION OF LOCATION ############
    # TODO: Improve this to use the guessedval/trueval arrays directly instead
    # of the confusion matrix, and improve the definition of the bias?
    if performance_with_location:
        show_performance_with_location(analysis, estimator_types)
                   
    ############## SHOW HOW PERFORMANCE IMPROVES WITH MORE CELLS ############
    if performance_with_cells:
        fig_performance_with_cells(analysis, estimator_types)
        
    ############## SHOW HOW PERFORMANCE DEPENDS ON FREQUENCY POOLS ############
    if performance_with_frequency:
        fig_performance_with_frequency(analysis, estimator_types,
                                       num_frequency_pools,
                                       frequency_pool_width,
                                       num_frequency_points)

if __name__=='__main__':
    
    # change this to change the model
    #from models.joris_cat import *
    #from models.joris_tollin_cat import *
    #from models.ircam_human_uniform_bipd import *
    from models.mcalpine_guinea_pig import *
    #from models.ircam_mcalpine_guinea_pig import *
        
    show_summary_info = True
    force_summary_compute = False
    show_error_matrix = True
    performance_with_location = False
    performance_with_cells = False
    performance_with_acoustic_noise = False
    performance_with_stimuli_alpha = False
    performance_with_frequency = False
    num_frequency_pools = 3
    frequency_pool_width = 40
    num_frequency_points = 10
    #limit_frequency = (500*Hz, Inf)
    limit_frequency = (0*Hz, 900*Hz)
    use_ideal_responses = False # Removes all response noise from the results
    show_responses = 0
    num_shuffles = 5
    training_size = 400
    testing_size = 200
    extra_size = 0#4400
    forced_size = 0#6400-training_size-testing_size
    cpu = 2 # number of CPUs, set to 0 to use all, or -1 to use all but one
    report = 'stderr'
    use_ram_disk = None
        
    ### STIMULI
    stimuli = WhiteNoiseStimuli(space, duration=100*ms, level=level)
    #stimuli = PowerLawNoiseStimuli(space, duration=100*ms, level=level)
    #stimuli = NatureSoundsStimuli(space, duration=100*ms, level=level, noiselevel=0*dB)
    
    ### ACOUSTIC NOISE
    acousticnoisemodel = NoAcousticNoise()
    #acousticnoisemodel = IndependentWhiteAcousticNoise(tuple(level-sn for sn in sn_level_range[::-1]))
    extraname['acousticnoisemodel'] = acousticnoisemodel
    
    ### FILTERS CONTROLLING WHAT TO USE FOR TRAINING AND TESTING DATA
    # Results are randomly shuffled and then each filter is applied in turn
    # to the whole results set, any data used for training will not be used
    # for testing
    training_filters = (
        'type=="whitenoise"',
        #'alpha<.75',
        #'"acoustic_noise_level" not in locals()',
        #'"acoustic_noise_level" in locals()',
        #'acoustic_noise_level<level-90*dB',
        #'rand()<.5',
        'i<training_size',
        )
    testing_filters = (
        'type=="whitenoise"',
        #'type=="powerlawnoise"',
        #'alpha>1.5',
        #'type=="natural"',
        #'"acoustic_noise_level" in locals()',
        #'acoustic_noise_level>level-40*dB',
        #'acoustic_noise_level>level',
        'i<testing_size',
        )
    sample_filters = training_filters
    #sample_filters = testing_filters

    estimator_types = (
    ##### JEFFRESS LIKE #################
        #(MakeEstimator(Jeffress), 'Naive Jeffress'),
        #(MakeEstimator(Jeffress, SmoothedMax(0.02*space.itd_max), phasemultiply=True), 'Smoothed Jeffress (phasemult)'),
        #(MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max), samplefigs=0), 'Smoothed Jeffress'),
        #(MakeEstimator(Jeffress, LLRSmoothedMax(0.3*space.itd_max), samplefigs=0), 'LLR Jeffress'),
        #(MakeEstimator(Jeffress, LQRSmoothedMax(0.6*space.itd_max), samplefigs=0), 'LQR Jeffress (long)'),
        #(MakeEstimator(Jeffress, LQRSmoothedMax(0.15*space.itd_max), samplefigs=0), 'LQR Jeffress (short)'),
        # Use for IRCAM IPD
        #(MakeEstimator(TrainedJeffress, DiscreteMax(window_len=13), bdmaximiser=DiscreteMax(13)), 'Jeffress (trained, smoothed, discrete)'),
        #(MakeEstimator(TrainedJeffress, LLRSmoothedMax(0.2*space.itd_max, discrete=True), bdmaximiser=LLRSmoothedMax(0.2*space.itd_max, discrete=True), samplefigs=5), 'Jeffress (trained, llr)'),
        #(MakeEstimator(TrainedJeffress, SmoothedMax(0.05*space.itd_max), bdmaximiser=FitSmoothedMax(0.05*space.itd_max, neighbourhood=0.075), samplefigs=2), 'Jeffress (trained)'),
        #(MakeEstimator(TrainedJeffress, SmoothedMax(0.05*space.itd_max), bdmaximiser=FitSmoothedMax(0.05*space.itd_max, neighbourhood=0.075), phasemultiply=True), 'Jeffress (trained, pm)'),
            # Use next TrainedJeffress for Tollin cat
        #(MakeEstimator(TrainedJeffress, SmoothedMax(0.05*space.itd_max), bdmaximiser=FitSmoothedMax(0.05*space.itd_max, neighbourhood=0.15)), 'Jeffress (trained)'),
    ##### DISTRIBUTION FITTING ###########
        #(MakeEstimator(FitDistributionMLE, SmoothedMax(0.15*space.itd_max)), 'Fit dist MLE (smoothing)'),
        #(MakeEstimator(FitDistributionMLE, FitSmoothedMax(0.15*space.itd_max)), 'Fit dist MLE (smoothing+fit)'),
    ##### PATTERN MATCHING ###############
        (MakeEstimator(PatternMatch), 'Pattern match'),
        #(MakeEstimator(PatternMatch, normalise=True, samplefigs=2), 'Pattern match (norm)'),
        #(MakeEstimator(PatternMatch, FitSmoothedMax(0.15*space.itd_max, degmax=3, neighbourhood=0.3), samplefigs=2), 'Pattern match (fit)'), # Better version
    ##### TWO CHANNEL ####################
        (MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
        #(MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend, difference=True), 'Two channel (diff)'),
        (MakeEstimator(TwoChannelCrossFrequency, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel (xf)'),
        (MakeEstimator(TwoChannelBanded, bandsize=40), 'Two channel (banded)'),
    ##### VECTOR/CENTROID ################
        #(MakeEstimator(Vector, PolyClosest(4)), 'Vector'), # BEST
    ##### SCIKITS.LEARN REGRESSION #######
      ### Linear Regression ###
        #(MakeEstimator(ScikitsLearnEstimator, LinearRegression()), 'Linear regression'), # OK
        #(MakeEstimator(ScikitsLearnEstimator, RidgeCV()), 'Ridge regression CV'), # GOOD
      ### Nearest neighbours ###
        #(MakeEstimator(ScikitsLearnEstimator, NeighborsRegressor()), 'Nearest neighb'), # Excellent
        )

    extended_name = generate_extended_datamanager_name(basename, extraname)
    
    analysis = get_analysis_from_namespace()
    
    ############################################################################
    ########## END OF SETTINGS, BEGINNING OF COMPUTATIONS ######################
    ############################################################################

    if acousticnoisemodel.__class__==NoAcousticNoise:
        performance_with_acoustic_noise = False

    print 'Model base name:', basename

    # TODO: remove this? Not quite so simple because we need to generate the
    # data sometimes...
    analysis.generate_data(stimuli, cpu, extra_size, forced_size)

    fm = analysis.fm
    fm_orig = analysis.fm_orig
    
    if testing_size==0:
        exit()

    # Compute number of tasks we will do
    numtasks = 0
    curtask = 0
    if show_summary_info:
        numtasks += len(estimator_types)
    if any([show_error_matrix, performance_with_location, performance_with_acoustic_noise,
            performance_with_stimuli_alpha]):
        numtasks += len(estimator_types)*num_shuffles
    if performance_with_cells:
        Nbinuserange = analysis.default_Nbinuserange
        numtasks += len(estimator_types)*num_shuffles*len(Nbinuserange)
    if performance_with_frequency:
        if cfrepeat>1:
            ocf = fm_orig.centre_frequencies[:fm_orig.cfN]
            I = argsort(ocf)
            cfsorted = ocf[I]
            cfu = sorted(unique(cf))
            numtasks += (len(cfu)+1-num_frequency_pools)*len(estimator_types)*num_shuffles
        else:
            numtasks += num_frequency_points*num_shuffles*len(estimator_types)
    progress = ProgressReporter(report)
    print 'Evaluating estimators, %d evaluations total.' % numtasks
    analysis.updateprogress = updateprogress

    from figures import fig_performance_with_acoustic_noise, \
                        fig_performance_with_stimuli_alpha, \
                        show_confusion_matrices, \
                        show_performance_with_location, \
                        fig_performance_with_cells, \
                        fig_performance_with_frequency
    
    do_show_sample_results()
    do_show_summary_info()
    do_show_analysis()
    progress.finish()
    
    show()
