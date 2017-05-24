from base import *

def analyse_performance_with_stimuli_alpha(analysis, estimator):
    sn_level_range = analysis.moresettings['sn_level_range']
    bins = analysis.moresettings['alpha_bins']
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
        alpha = array([meta['alpha'] for _, _, meta, _ in testing_responses])
        serrors = (abs(trueval)-abs(guessedval))*second/usecond
        errors = abs(trueval-guessedval)*second/usecond
        sume, _ = histogram(alpha, bins, weights=errors)
        sumse, _ = histogram(alpha, bins, weights=serrors)
        ne, _ = histogram(alpha, bins)
        ne[ne==0] = 1
        meane = sume/ne
        meanse = sumse/ne
        summeane += meane
        sum2meane += meane**2
        summeanse += meanse
        sum2meanse += meanse**2
        sumXY, _ = histogram(alpha, bins, weights=trueval*guessedval)
        sumXX, _ = histogram(alpha, bins, weights=trueval**2)
        bias = 100*(1-sumXY/sumXX)
        biases.append(bias)
    meane = summeane/num_shuffles
    stde = sqrt(sum2meane/num_shuffles-meane**2)
    bias = array(biases)
    return binmids, meane, stde, mean(bias, axis=0), std(bias, axis=0) 

def fig_performance_with_stimuli_alpha(analysis, estimator_types,
                do_figure=True, do_suptitle=True, do_legend=True,
                show_mean_error=True, show_bias=True,
                show_xlabel=True, show_ylabel=True,
                formatting=None,
                tight=True,
                ):
    if formatting is None: formatting = dict((name, {}) for _, name in estimator_types)
    basename = analysis.settings['basename']
    if do_figure: figure()
    if do_suptitle: suptitle(basename)
    yminb, ymaxb = 1e10, -1e10
    ymine, ymaxe = 1e10, -1e10
    for f, name in estimator_types:
        estimator = f(analysis.fm_orig)
        binmids, meane, stde, meanb, stdb = analyse_performance_with_stimuli_alpha(analysis, estimator)
        if show_mean_error and show_bias:
            subplot(211)
        if show_mean_error:
            errorbar(binmids, meane, yerr=stde, label=name, **formatting[name])
            if amin(meane-stde)<ymine:
                ymine = amin(meane-stde)
            if amax(meane+stde)>ymaxe:
                ymaxe = amax(meane+stde)
        if show_mean_error and show_bias:
            subplot(212)
        if show_bias:
            errorbar(binmids, meanb, yerr=stdb, label=name, **formatting[name])
            if amin(meanb-stdb)<yminb:
                yminb = amin(meanb-stdb)
            if amax(meanb+stdb)>ymaxb:
                ymaxb = amax(meanb+stdb)
    if show_mean_error and show_bias:
        subplot(211)
        if do_legend:
            legend(loc='upper right', ncol=2)
    if show_mean_error:
        if tight:
            xlim(amin(binmids), amax(binmids))
            ylim(ymine, ymaxe)
        if show_xlabel: xlabel(r'Noise colour $\alpha$')
        if show_ylabel: ylabel('Mean error ($\\mu$s)')
    if show_mean_error and show_bias:
        subplot(212)
    if show_bias:
        if tight:
            xlim(amin(binmids), amax(binmids))
            ylim(yminb, ymaxb)
        axhline(0, ls='--', c='k')
        if show_xlabel: xlabel(r'Noise colour $\alpha$')
        if show_ylabel: ylabel('Bias to centre (%)')
