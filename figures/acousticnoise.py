from base import *

def analyse_performance_with_acoustic_noise(analysis, estimator):
    sn_level_range = analysis.moresettings['sn_level_range']
    num_noise_bins = analysis.moresettings['num_noise_bins']
    num_shuffles = analysis.settings['num_shuffles']
    level = analysis.moresettings['level']
    bins = linspace(*(sn_level_range+(num_noise_bins,)))
    binmids = 0.5*(bins[1:]+bins[:-1])
    shuffled_results = analysis(analysis.shuffled_results, estimator, num_shuffles)
    summeane, sum2meane = zeros(len(binmids)), zeros(len(binmids))
    biases = []
    for shufflenum in xrange(num_shuffles):
        results, trueval, guessedval, testing_responses = shuffled_results[shufflenum]
        acoustic_noise_levels = array([meta['acoustic_noise_level'] for _, _, meta, _ in testing_responses])
        sn = float(level)-acoustic_noise_levels
        errors = abs(trueval-guessedval)*second/usecond
        sume, _ = histogram(sn, bins, weights=errors)
        ne, _ = histogram(sn, bins)
        ne[ne==0] = 1
        meane = sume/ne
        summeane += meane
        sum2meane += meane**2
        sumXY, _ = histogram(sn, bins, weights=trueval*guessedval)
        sumXX, _ = histogram(sn, bins, weights=trueval**2)
        bias = 100*(1-sumXY/sumXX)
        biases.append(bias)
    meane = summeane/num_shuffles
    stde = sqrt(sum2meane/num_shuffles-meane**2)
    bias = array(biases)
    return binmids, meane, stde, mean(bias, axis=0), std(bias, axis=0)    

def fig_performance_with_acoustic_noise(analysis, estimator_types,
                                        do_figure=True,
                                        do_suptitle=True,
                                        show_mean_error=True,
                                        show_bias=True,
                                        show_legend=True,
                                        show_xlabel=True,
                                        show_ylabel=True,
                                        formatting=None,
                                        log_y_axis=False,
                                        tight=True,
                                        min_sn=-5*dB,
                                        ):
    if formatting is None: formatting = dict((name, {}) for _, name in estimator_types)
    basename = analysis.settings['basename']
    if do_figure: figure()
    if do_suptitle: suptitle(basename)
    for f, name in estimator_types:
        estimator = f(analysis.fm_orig)
        binmids, meane, stde, meanb, stdb = analyse_performance_with_acoustic_noise(analysis, estimator)
        I = binmids>min_sn
        binmids = binmids[I]
        meane = meane[I]
        stde = stde[I]
        meanb = meanb[I]
        stdb = stdb[I]
        if show_mean_error and show_bias: subplot(211)
        if show_mean_error:
            errorbar(binmids, meane, yerr=stde, label=name, **formatting[name])
            if log_y_axis:
                gca().set_yscale('log')
            if tight:
                xlim(float(min_sn), amax(binmids))
                ylim(amin(meane-stde), amax(meane+stde))
        if show_mean_error and show_bias: subplot(212)
        if show_bias:
            errorbar(binmids, meanb, yerr=stdb, label=name, **formatting[name])
            if tight:
                xlim(float(min_sn), amax(binmids))
                ylim(amin(meanb-stdb), amax(meanb+stdb))
    if show_mean_error and show_bias: subplot(211)
    if show_mean_error:
        if show_legend: legend(loc='upper right', ncol=2)
        if show_xlabel: xlabel('S/N (dB)')
        if show_ylabel: ylabel(r'Mean error ($\mu$s)')
    if show_mean_error and show_bias: subplot(212)
    if show_bias:
        axhline(0, ls='--', c='k')
        if show_xlabel: xlabel('S/N (dB)')
        if show_ylabel: ylabel('Bias to centre (%)')

if __name__=='__main__':
    from analyse_model import * # load default values of some parameters
    # change this to change the model
    #from models.artificial_twochan_bitd import *
    from models.joris_cat import *
    #from models.joris_cat_low_frequency import *
    #from models.joris_cat_uniform_bitd import *
    #from models.joris_cat_uniform_bipd import *
    #from models.ircam_human_uniform_bitd import *
    #from models.mcalpine_guinea_pig import *
    use_ideal_responses = True # Removes all response noise from the results
    num_shuffles = 5
    training_size = 800
    testing_size = 400
    acousticnoisemodel = IndependentWhiteAcousticNoise(tuple(level-sn for sn in sn_level_range[::-1]))
    extraname['acousticnoisemodel'] = acousticnoisemodel
    training_filters = (
        'type=="whitenoise"',
        #'"acoustic_noise_level" not in locals()',
        '"acoustic_noise_level" in locals()',
        'i<training_size',
        )
    testing_filters = (
        'type=="whitenoise"',
        '"acoustic_noise_level" in locals()',
        'i<testing_size',
        )
    estimator_types = (
    ##### JEFFRESS LIKE #################
        (MakeEstimator(Jeffress, SmoothedMax(0.02*space.itd_max), phasemultiply=True), 'Smoothed Jeffress (phasemult)'),
        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed Jeffress'),
        #(MakeEstimator(TrainedJeffress, SmoothedMax(0.05*space.itd_max), bdmaximiser=FitSmoothedMax(0.05*space.itd_max, neighbourhood=0.075)), 'TJeff'),
    ##### PATTERN MATCHING ###############
        (MakeEstimator(PatternMatch), 'Pattern match'),
    ##### TWO CHANNEL ####################
        #(MakeEstimator(TwoChannel, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
    ##### SCIKITS.LEARN REGRESSION #######
      ### Linear Regression ###
        (MakeEstimator(ScikitsLearnEstimator, LinearRegression()), 'Linear regression'), # OK
        (MakeEstimator(ScikitsLearnEstimator, RidgeCV()), 'Ridge regression CV'), # GOOD
      ### Nearest neighbours ###
        (MakeEstimator(ScikitsLearnEstimator, NeighborsRegressor()), 'Nearest neighb'), # Excellent
        )

    analysis = get_analysis_from_namespace()
    
    fig_performance_with_acoustic_noise(analysis, estimator_types)
    show()
    