from base import *

def analyse_performance_with_stimuli_tonefreq(analysis, estimator,
                                              tonefreqkey='tonefreq',
                                              fraction=1.0):
    sn_level_range = analysis.moresettings['sn_level_range']
    bins = analysis.moresettings['tonefreq_bins']
    num_shuffles = analysis.settings['num_shuffles']
    shuffled_results = analysis(analysis.shuffled_results, estimator, num_shuffles)
    binmids = 0.5*(bins[1:]+bins[:-1])
    summeane = zeros(len(binmids))
    sum2meane = zeros(len(binmids))
    summeanse = zeros(len(binmids))
    sum2meanse = zeros(len(binmids))
    biases = []
    for shufflenum in xrange(num_shuffles):
        results, trueval, guessedval, testing_responses = shuffled_results[shufflenum]
        I = abs(trueval)<amax(abs(trueval))*fraction
        trueval = trueval[I]
        guessedval = guessedval[I]
        testing_responses = [testing_responses[i] for i in xrange(len(I)) if I[i]]
        tonefreq = array([meta[tonefreqkey] for _, _, meta, _ in testing_responses])
#        figure()        
#        for i, (b1, b2) in enumerate(zip(bins[:-1], bins[1:])):
#            subplot(3, 3, i+1)
#            J = (b1<tonefreq)&(tonefreq<b2)
#            X = trueval[J]
#            Y = guessedval[J]
#            plot(X, Y, '.')
#            plot([amin(X), amax(X)], [amin(X), amax(X)], '-k')
#            g = sum(X*Y)/sum(X*X)
#            plot([amin(X), amax(X)], [g*amin(X), g*amax(X)], '-r')
#            xticks([])
#            yticks([])
#        show()
#        exit()
        serrors = (abs(trueval)-abs(guessedval))*second/usecond
        errors = abs(trueval-guessedval)*second/usecond
        sume, _ = histogram(tonefreq, bins, weights=errors)
        sumse, _ = histogram(tonefreq, bins, weights=serrors)
        ne, _ = histogram(tonefreq, bins)
        ne[ne==0] = 1
        meane = sume/ne
        meanse = sumse/ne
        summeane += meane
        sum2meane += meane**2
        summeanse += meanse
        sum2meanse += meanse**2
        sumXY, _ = histogram(tonefreq, bins, weights=trueval*guessedval)
        sumXX, _ = histogram(tonefreq, bins, weights=trueval**2)    
        bias = 100*(1-sumXY/sumXX)
        biases.append(bias)
    meane = summeane/num_shuffles
    stde = sqrt(sum2meane/num_shuffles-meane**2)
    bias = array(biases)
    return binmids, meane, stde, mean(bias, axis=0), std(bias, axis=0) 

def fig_performance_with_stimuli_tonefreq(analysis, estimator_types,
               tonefreqkey='tonefreq',
               fraction=1.0,
               do_figure=True,
               do_suptitle=True,
               show_mean_error=True,
               show_bias=True,
               show_legend=True,
               show_xlabel=True,
               show_ylabel=True,
               frequency_unit=Hz,
               formatting=None,
               axhline_formatting=None,
               fill_between_alpha=0,
               fill_between_pairs=[],
               tight=True,
               ):
    basename = analysis.settings['basename']
    if do_figure:
        figure()
    if do_suptitle:
        suptitle(basename)
    if show_bias:
        if show_mean_error:
            subplot(211)
        if axhline_formatting is None:
            axhline_formatting = {'ls':'--', 'c':'k', 'lw':0.5}
        axhline(0, **axhline_formatting)
    all_binmids = []
    all_meane = []
    all_meanb = []
    all_c = []
    yminb, ymaxb = 1e10, -1e10
    ymine, ymaxe = 1e10, -1e10
    for f, name in estimator_types:
        if formatting is None:
            F = {}
        else:
            F = formatting[name]
        estimator = f(analysis.fm_orig)
        binmids, meane, stde, meanb, stdb = analyse_performance_with_stimuli_tonefreq(analysis, estimator, tonefreqkey, fraction=fraction)
        binmids = array(binmids)/frequency_unit
        all_binmids.append(binmids)
        if 'color' in F:
            all_c.append(F['color'])
        if show_mean_error and show_bias:
            subplot(211)
        if show_mean_error:
            errorbar(binmids, meane, yerr=stde, label=name, **F)
            all_meane.append(meane)
            if amin(meane-stde)<ymine:
                ymine = amin(meane-stde)
            if amax(meane+stde)>ymaxe:
                ymaxe = amax(meane+stde)
        if show_mean_error and show_bias:
            subplot(212)
        if show_bias:
            errorbar(binmids, meanb, yerr=stdb, label=name, **F)
            all_meanb.append(meanb)
            if amin(meanb-stdb)<yminb:
                yminb = amin(meanb-stdb)
            if amax(meanb+stdb)>ymaxb:
                ymaxb = amax(meanb+stdb)
    if len(all_c)==0:
        all_c = [(0.7, 0.7, 0.7)]*len(all_binmids)
    if fill_between_alpha and fill_between_pairs:
        for i, j in fill_between_pairs:
            c = all_c[i]+(fill_between_alpha,)
            if show_mean_error and show_bias:
                subplot(211)
            if show_mean_error:
                fill_between(binmids, all_meane[i], all_meane[j], color=c)
            if show_mean_error and show_bias:
                subplot(212)
            if show_bias:
                fill_between(binmids, all_meanb[i], all_meanb[j], color=c)

    if show_mean_error and show_bias:
        subplot(211)
    if show_mean_error:
        if tight:
            xlim(amin(binmids), amax(binmids))
            ylim(ymine, ymaxe)
        if show_legend: legend(loc='upper right', ncol=2)
        if show_xlabel: xlabel('Frequency (%s)'%frequency_unit)
        if show_ylabel: ylabel('Mean error ($\\mu$s)')
    if show_mean_error and show_bias:
        subplot(212)
    if show_bias:
        if tight:
            xlim(amin(binmids), amax(binmids))
            ylim(yminb, ymaxb)
        if show_xlabel: xlabel('Frequency (%s)'%frequency_unit)
        if show_ylabel: ylabel('Bias to centre (%)')
    return ymine, ymaxe, yminb, ymaxb

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
    training_size = 800
    testing_size = 800
#    limit_frequency = (0*Hz, 900*Hz)
    acousticnoisemodel = IndependentWhiteAcousticNoise((level, level))
    extraname['acousticnoisemodel'] = acousticnoisemodel
    training_filters = (
        'type=="whitenoise"',
        'i<training_size',
        )
    testing_filters = (
        'type=="tone"',
        #'type=="bandpassnoise"',
        'i<testing_size',
        )
    estimator_types = (
        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed Jeffress'),
        (MakeEstimator(PatternMatch), 'Pattern match'),
        #(MakeEstimator(PatternMatch, normalise_banded=1), 'Pattern match (band=1)'),
        #(MakeEstimator(PatternMatch, normalise_banded=2), 'Pattern match (band=2)'),
        #(MakeEstimator(PatternMatch, normalise_banded=5), 'Pattern match (band=5)'),
        #(MakeEstimator(PatternMatch, normalise_banded=10), 'Pattern match (band=10)'),
        #(MakeEstimator(PatternMatch, normalise_banded=40), 'Pattern match (band=40)'),
        #(MakeEstimator(PatternMatch, normalise_banded=60), 'Pattern match (band=80)'),
        (MakeEstimator(TwoChannel, PolyClosest(6), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
#        (MakeEstimator(TwoChannelCrossFrequency, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel (xf)'),
#        (MakeEstimator(Vector, PolyClosest(4)), 'Vector'),
        )

    analysis = get_analysis_from_namespace()
    
    fig_performance_with_stimuli_tonefreq(analysis, estimator_types,
                                          #'bpcfreq',
                                          fraction=0.5,
                                          )
    show()
    