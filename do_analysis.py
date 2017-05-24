from analyse_model import *
import analyse_model
   
def putinnamespace():
    for k, v in globals().iteritems():
        if not k.startswith('_'):
            setattr(analyse_model, k, v)
   
if __name__=='__main__':
    # defaults
    modelname = 'joris_cat'
    analysis_type = 'performance'
    #analysis_type = 'noise'
    #analysis_type = 'colour'
    training_size = 400
    testing_size = 200
    num_shuffles = 5

    num_frequency_pools = 3
    frequency_pool_width = 40
    num_frequency_points = 10
    use_ideal_responses = False # Removes all response noise from the results

    show_summary_info = True
    force_summary_compute = False
    show_error_matrix = False
    performance_with_location = False
    performance_with_cells = False
    performance_with_acoustic_noise = False
    performance_with_stimuli_alpha = False
    performance_with_frequency = False
    show_responses = 0
    report = 'stderr'
    use_ram_disk = None
    
    # overrides
    for arg in sys.argv[1:]:
        exec arg
        
    if analysis_type=='cells':
        performance_with_cells = True
        acoustic_noise_type = 'none'
    elif analysis_type=='frequency':
        performance_with_frequency = True
        acoustic_noise_type = 'none'
    elif analysis_type=='noise':
        performance_with_acoustic_noise = True
        acoustic_noise_type = 'independent'
    elif analysis_type=='colour':
        performance_with_stimuli_alpha = True
        acoustic_noise_type = 'none'
    else:
        raise ValueError('Unknown analysis type.')

    for arg in sys.argv[1:]:
        exec arg

    _limit_frequency = limit_frequency
    exec 'from models.'+modelname+' import *'
    limit_frequency = _limit_frequency        

    ### ACOUSTIC NOISE
    if acoustic_noise_type=='none':
        acousticnoisemodel = NoAcousticNoise()
    elif acoustic_noise_type=='independent':
        acousticnoisemodel = IndependentWhiteAcousticNoise(tuple(level-sn for sn in sn_level_range[::-1]))
    else:
        raise ValueError('acoustic_noise_type unknown')
    extraname['acousticnoisemodel'] = acousticnoisemodel
        
    sample_filters = []

    if modelname=='joris_cat':
        estimator_types = (
            (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed Jeffress'),
            (MakeEstimator(PatternMatch), 'Pattern match'),
            (MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
            (MakeEstimator(TwoChannelCrossFrequency, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel (xf)'),
            (MakeEstimator(Vector, PolyClosest(4)), 'Vector'),
            )
    elif modelname=='ircam_human_uniform_bipd':
        winlen = 13
#        if analysis_type=='frequency':
#            winlen = 5
        estimator_types = (
            (MakeEstimator(TrainedJeffress, DiscreteMax(window_len=winlen), bdmaximiser=DiscreteMax(winlen)), 'Jeffress'),
            (MakeEstimator(PatternMatch), 'Pattern match'),
            (MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
            )
    elif modelname=='ircam_mcalpine_guinea_pig':
        winlen = 13
#        if analysis_type=='frequency':
#            winlen = 3
        estimator_types = (
            (MakeEstimator(TrainedJeffress, DiscreteMax(window_len=winlen), bdmaximiser=DiscreteMax(winlen)), 'Jeffress'),
            (MakeEstimator(PatternMatch), 'Pattern match'),
            (MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
            )
    elif modelname=='joris_tollin_cat':
        estimator_types = (
            (MakeEstimator(TrainedJeffress, SmoothedMax(0.05*space.itd_max), bdmaximiser=FitSmoothedMax(0.05*space.itd_max, neighbourhood=0.15)), 'Jeffress'),
            (MakeEstimator(PatternMatch), 'Pattern match'),
            (MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
            )
    elif modelname=='mcalpine_guinea_pig':
        estimator_types = (
            (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed Jeffress'),
            (MakeEstimator(PatternMatch), 'Pattern match'),
            (MakeEstimator(PatternMatch, normalise_banded=40), 'Pattern match banded'),
            (MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
            (MakeEstimator(TwoChannelCrossFrequency, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel (xf)'),
            (MakeEstimator(Vector, PolyClosest(4)), 'Vector'),
            )
    else:
        raise ValueError('Unknown modelname for selecting estimators.')

    extended_name = generate_extended_datamanager_name(basename, extraname)
       
    ############################################################################
    ########## END OF SETTINGS, BEGINNING OF COMPUTATIONS ######################
    ############################################################################

    from figures import fig_performance_with_acoustic_noise, \
                        fig_performance_with_stimuli_alpha, \
                        show_confusion_matrices, \
                        show_performance_with_location, \
                        fig_performance_with_cells, \
                        fig_performance_with_frequency

    if analysis_type=='cells':
        training_filters = (
            'type=="whitenoise"',
            'i<training_size',
            )
        testing_filters = (
            'type=="whitenoise"',
            'i<testing_size',
            )
        analysis = get_analysis_from_namespace()
        fm = analysis.fm
        putinnamespace()
        do_show_sample_results()
        do_show_summary_info()
        do_show_analysis()
    elif analysis_type=='frequency':
        training_filters = (
            'type=="whitenoise"',
            'i<training_size',
            )
        testing_filters = (
            'type=="whitenoise"',
            'i<testing_size',
            )
        analysis = get_analysis_from_namespace()
        fm = analysis.fm
        putinnamespace()
        do_show_analysis()
    elif analysis_type=='noise':
        training_filters = (
            'type=="whitenoise"',
            '"acoustic_noise_level" not in locals()',
            'i<training_size',
            )
        testing_filters = (
            'type=="whitenoise"',
            '"acoustic_noise_level" in locals()',
            'i<testing_size',
            )
        analysis = get_analysis_from_namespace()
        putinnamespace()
        do_show_analysis()
        suptitle('training in quiet', x=.8)
#        training_filters = (
#            'type=="whitenoise"',
#            '"acoustic_noise_level" in locals()',
#            'i<training_size',
#            )
#        analysis = get_analysis_from_namespace()
#        putinnamespace()
#        do_show_analysis()
#        suptitle('training in noise', x=.8)
    elif analysis_type=='colour':
        training_filters = (
            'type=="whitenoise"',
            'i<training_size',
            )
        testing_filters = (
            'type=="powerlawnoise"',
            'i<testing_size',
            )
        use_ideal_responses = False
        analysis = get_analysis_from_namespace()
        putinnamespace()
        do_show_analysis()
        suptitle('(neural noise)', x=.8)
        use_ideal_responses = True
        analysis = get_analysis_from_namespace()
        putinnamespace()
        do_show_analysis()
        suptitle('(no neural noise)', x=.8)
    else:
        raise ValueError('Unknown analysis type.')
        
    show()
