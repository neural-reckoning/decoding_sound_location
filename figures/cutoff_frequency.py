#import os
#os.environ['SAMPLERATE'] = '97656' # cat
from base import *
from scipy import stats
from external_data.number_of_mso_neurons import *

def analyse_performance_with_cutoff_frequency(analysis, estimator_type, angles=False):
    Nbinuserange = analysis.default_Nbinuserange
    meanerror, cutoff = analysis(analysis.cutoff_frequency_results, estimator_type, angles=angles)
    error = []
    error_std = []
    for Nbinuse in Nbinuserange:
        me = meanerror[Nbinuse]
        me = [e for e in me if not isnan(e)]
        if len(me):
            mestd = std(me)
            me = mean(me)
        else:
            mestd = nan
            me = nan
        error.append(me)
        error_std.append(mestd)
    return Nbinuserange, array(error), array(error_std), array(cutoff)

def fig_performance_with_cutoff_frequency(analysis, estimator_types,
                                          angles=False,
                                          formatting=None,
                                          ):
    if formatting is None: formatting = dict((name, {}) for _, name in estimator_types)
    if angles:
        errorunit = 1
    else:
        errorunit = usecond
    if angles:
        axhline(60, ls='--', color='k')
    else:
        axhline((2./3)*float(analysis.settings['space'].itd_max/usecond),
                ls='--', color='k')
    for estimator_type, name in estimator_types:
        Nbinuserange, error, error_std, cutoff = analyse_performance_with_cutoff_frequency(analysis,
                                                estimator_type, angles=angles)
        errorbar(cutoff/kHz, error/errorunit, error_std/errorunit, label=name,
                 **formatting[name])
    #legend(loc='upper right', ncol=2)
    xlabel('Cutoff frequency (kHz)')
    if angles:
        ylabel('Mean error (deg)')
    else:
        ylabel(r'Mean error ($\mu$s)')
    ylim(ymin=0)

if __name__=='__main__':
    from analyse_model import * # load default values of some parameters
    # change this to change the model
    #from models.joris_cat import *
    #from models.joris_tollin_cat import *
    #from models.ircam_human_uniform_bitd import *
    #from models.ircam_human_uniform_bipd import *
    from models.mcalpine_guinea_pig import *
    use_ideal_responses = False # Removes all response noise from the results
    num_shuffles = 5
    training_size = 400
    testing_size = 200
    acousticnoisemodel = NoAcousticNoise()
    extraname['acousticnoisemodel'] = acousticnoisemodel
    training_filters = (
        'type=="whitenoise"',
        'i<training_size',
        )
    testing_filters = (
        'type=="whitenoise"',
        'i<testing_size',
        )
    estimator_types = (
        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed peak'),
        (MakeEstimator(PatternMatch), 'Pattern match'),
        (MakeEstimator(TwoChannel, PolyClosest(6), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
        )

    analysis = get_analysis_from_namespace()
    
    fig_performance_with_cutoff_frequency(analysis, estimator_types)

    show()
    