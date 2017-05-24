from base import *
from confusionmatrices import *

def analyse_performance_with_frequency(analysis, estimator_type,
                                       num_frequency_pools,
                                       frequency_pool_width,
                                       num_frequency_points,
                                       bandpass=False,
                                       bias_fraction=1,
                                       ):
    num_shuffles = analysis.settings['num_shuffles']
    if bandpass:
        f = analysis.artificial_bandpass_results
    else:
        f = analysis.frequency_results
    allcfs, fmid, E, B, FR = analysis(f, estimator_type,
                              num_frequency_pools, frequency_pool_width,
                              num_frequency_points, bias_fraction=bias_fraction)
    E = array([E[shufflenum] for shufflenum in xrange(num_shuffles)])
    B = array([B[shufflenum] for shufflenum in xrange(num_shuffles)])
    return allcfs, fmid, mean(E, axis=0), std(E, axis=0), mean(B, axis=0), std(B, axis=0)

def fig_performance_with_frequency(analysis, estimator_types,
                                   num_frequency_pools, frequency_pool_width,
                                   num_frequency_points,
                                   bias_fraction=1,
                                   bandpass=False,
                                   do_figure=True,
                                   do_suptitle=True,
                                   show_mean_error=True,
                                   show_bias=True,
                                   show_legend=True,
                                   show_freq_bars=True,
                                   formatting=None,
                                   frequency_unit=Hz,
                                   ):
    if formatting is None: formatting = dict((name, {}) for _, name in estimator_types)
    basename = analysis.settings['basename']
    num_shuffles = analysis.settings['num_shuffles']
    if do_figure: figure()
    if do_suptitle: suptitle(basename)#+' (%d pools)'%num_frequency_pools)
    plotlines = {}
    for estimator_type, name in estimator_types:
        allcfs, fmid, me, se, mb, sb = analyse_performance_with_frequency(analysis,
                                        estimator_type, num_frequency_pools,
                                        frequency_pool_width,
                                        num_frequency_points, bandpass=bandpass,
                                        bias_fraction=bias_fraction)
        fmid = []
        fwid = []
        for cf in allcfs:
            fmid.append((amax(cf)+amin(cf))/2)
            fwid.append((amax(cf)-amin(cf))/2)
        fmid = array(fmid)/frequency_unit
        plotlines[name] = (fmid, me, se, fwid)
        if show_mean_error and show_bias: subplot(211)
        if show_mean_error:
            if show_freq_bars:
                errorbar(fmid, me/usecond, se/usecond, fwid, label=name,
                         **formatting[name])
            else:
                errorbar(fmid, me/usecond, se/usecond, label=name,
                         **formatting[name])
        if show_mean_error and show_bias: subplot(212)
        if show_bias:
            if show_freq_bars:
                errorbar(fmid, mb, sb, fwid, label=name, **formatting[name])
            else:
                errorbar(fmid, mb, sb, label=name, **formatting[name])
    if show_mean_error and show_bias: subplot(211)
    if show_mean_error:
        xlabel('Frequency (%s)'%frequency_unit)
        ylabel(r'Mean error ($\mu$s)')
        if show_legend: legend(loc='upper right', ncol=2)
    if show_mean_error and show_bias: subplot(212)
    if show_bias:
        xlabel('Frequency (%s)'%frequency_unit)
        ylabel('Bias to centre (%)')
        axhline(0, ls='--', c='k')
        if show_legend and not show_mean_error:
            legend(loc='upper right', ncol=2)
    return plotlines

if __name__=='__main__':
    from models.mcalpine_guinea_pig import *
    #from models.joris_cat import *
    from analyse_model import *
    use_ideal_responses = False # Removes all response noise from the results
    num_shuffles = 5
    training_size = 400
    testing_size = 200
    num_frequency_pools = 3
    frequency_pool_width = (80, 40)
    #limit_frequency = (0*Hz, 900*Hz)
    num_frequency_points = 40
    training_filters = (
        'type=="whitenoise"',
        'i<training_size',
        )
    testing_filters = (
        'type=="whitenoise"',
        'i<testing_size',
        )
    acousticnoisemodel = NoAcousticNoise()
    extraname['acousticnoisemodel'] = acousticnoisemodel
    analysis = get_analysis_from_namespace()
    estimator_types = (
        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed Jeffress'),
        (MakeEstimator(PatternMatch), 'Pattern match'),
        (MakeEstimator(TwoChannel, PolyClosest(6), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
        #(MakeEstimator(TwoChannelCrossFrequency, PolyClosest(3), itdmax_extend=itdmax_extend), 'Two channel (xf)'),
        #(MakeEstimator(Vector, PolyClosest(4)), 'Vector'),
        )
    fig_performance_with_frequency(analysis, estimator_types,
                                   num_frequency_pools, frequency_pool_width,
                                   num_frequency_points,
                                   show_freq_bars=False,
                                   bandpass=False,
                                   bias_fraction=0.6,
                                   )
    show()
    