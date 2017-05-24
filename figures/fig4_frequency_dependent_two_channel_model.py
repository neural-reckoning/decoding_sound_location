from base import *
from frequency import *
from cells import *
from stimulialpha import *
from tonefreq import *
from binauraldisplay import *
from analyse_model import *
from colours import *
from matplotlib import cm

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.1
rcParams['figure.subplot.right'] = .98
rcParams['figure.subplot.bottom'] = .1
rcParams['figure.subplot.top'] = .91
rcParams['figure.subplot.wspace'] = 0.4
rcParams['figure.subplot.hspace'] = 0.3
# figsize is width x height in inches
figure(figsize=(7, 6), dpi=80)

show_tuning_curves = True
show_two_channel_examples = True
show_cross_freq_examples = True
show_frequency_plot = True
show_bandpass_plot = True
#bandpass_type = 'Tones'
bandpass_type = 'Band pass'
num_frequency_bands = 5
scale_bands = False # 0.5, True, False
#use_low_frequencies = False
colour_gradient = lambda c: (c, 0.3*(1-c), 1-c*0.5)
use_smooth = True
fill_between_alpha = 0.5

from models.mcalpine_guinea_pig import *

use_ideal_responses = False # Removes all response noise from the results
num_shuffles = 25
training_size = 400
testing_size = 800
num_frequency_pools = 3
frequency_pool_width = 40
num_frequency_points = 10

the_limit_frequency = (0*Hz, 1200*Hz)

all_estimator_types = {
    'Two channel': MakeEstimator(TwoChannel, PolyClosest(6), difference=True,
                                 itdmax_extend=itdmax_extend),
    'Two channel cross-frequency': MakeEstimator(TwoChannelCrossFrequency,
                                                 PolyClosest(6),
                                                 itdmax_extend=itdmax_extend),
    'Two channel cross-frequency alt': MakeEstimator(TwoChannelCrossFrequencyAlt,
                                                     PolyClosest(6),
                                                     itdmax_extend=itdmax_extend),
#    'Two channel banded': MakeEstimator(TwoChannelBanded),
#    'Two channel banded alt 10': MakeEstimator(TwoChannelBandedAlt, PolyClosest(6), bandsize=10, itdmax_extend=itdmax_extend),
    'Pattern match': MakeEstimator(PatternMatch),
    'Pattern match banded':MakeEstimator(PatternMatch, normalise_banded=40),
    }

formatting = {
    'Two channel':{'color':estimator_colours['twochannel'], 'ls':'--', 'lw':2},
#    'Two channel cross-frequency':{'color':estimator_colours['twochannel_xf'], 'ls':'--', 'lw':2},
    'Two channel cross-frequency alt':{'color':estimator_colours['twochannel'], 'ls':'-', 'lw':2},
#    'Two channel banded':{'color':'m'},
#    'Two channel banded alt 10':{'color':'c'},
    'Pattern match':{'color':estimator_colours['patternmatch'], 'ls':'--', 'lw':2},
    'Pattern match banded':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':2},
    }

the_filters = {
    'whitenoise': (('type=="whitenoise"', 'i<training_size'),
                   ('type=="whitenoise"', 'i<testing_size')),
    'Tones':       (('type=="whitenoise"', 'i<training_size'),
                   ('type=="tone"', 'i<testing_size')),
    'Band pass':   (('type=="whitenoise"', 'i<training_size'),
                   ('type=="bandpassnoise"', 'i<testing_size')),
   }

bandpass_keys = {
    'Tones':'tonefreq',
    'Band pass':'bpcfreq',
    }
bandpass_key = bandpass_keys[bandpass_type]

acousticnoisemodel = NoAcousticNoise()
extraname['acousticnoisemodel'] = acousticnoisemodel

@figcache
def get_sample_tuning_curves(num_tuning_curves):
    fm = analysis.fm_orig
    cf = fm.centre_frequencies[:len(fm.centre_frequencies)/2]
    I = argsort(cf)
    cfI = cf[I]
    training_responses, testing_responses = analysis.traintest_responses(None)
    R = []
    locations = []
    for response, location, meta, _ in training_responses:
        R.append(response/meta['duration'])
        locations.append(location)
    R = array(R)
    locations = array(locations)
    subplot(2, 3, 1)
    J = arange(0, len(I), len(I)/num_tuning_curves)
    return locations, R[:, I[J]], cf[I[J]]

@figcache
def get_ratio_plot_data(use_cross_frequency_two_channel,
                        num_frequency_bands,
                        scale_bands=False,
                        ):
    if use_cross_frequency_two_channel:
#        estimator_type = MakeEstimator(TwoChannelCrossFrequency,
#                                       PolyClosest(6),
#                                       itdmax_extend=itdmax_extend)
        estimator_type = MakeEstimator(TwoChannelCrossFrequencyAlt,
                                       PolyClosest(6),
                                       itdmax_extend=itdmax_extend)
#        estimator_type = MakeEstimator(TwoChannelBandedAlt,
#                                       PolyClosest(6),
#                                       bandsize=10,
#                                       itdmax_extend=itdmax_extend)
    else:
        estimator_type = MakeEstimator(TwoChannel, PolyClosest(6),
                                       difference=True,
                                       itdmax_extend=itdmax_extend)
    fm = analysis.fm_orig
    cf = fm.centre_frequencies[:len(fm.centre_frequencies)/2]
    I = argsort(cf)
    II = argsort(I)
    cfI = cf[I]
    estimator = estimator_type(fm)
    training_responses, testing_responses = analysis.traintest_responses(None)
    estimator.train(responses=testing_responses)
#    coeffs = estimator.coefficients
    locations = []
    ratios = []
    freqs = []
    for response, location, meta in estimator.data:
        response = meta['noiseless_response']
        response = response[I]
#        coefficients = estimator.coefficients[I]
#        rc = response*coefficients
        J = hstack((arange(0, len(cf), len(cf)/num_frequency_bands), len(cf)))
        for jstart, jend in zip(J[:-1], J[1:]):
            f = mean(cfI[jstart:jend])*Hz
#            denominator = sum(response[jstart:jend])
#            partialratio = sum(rc[jstart:jend])/denominator
            R = zeros(len(response))
            R[jstart:jend] = response[jstart:jend]
            R = R[II]
            partialratio = estimator.func(R)
            locations.append(location)
            ratios.append(partialratio)
            freqs.append(f)
    locations = array(locations)
    ratios = array(ratios)
    freqs = array(freqs)
    if scale_bands and use_cross_frequency_two_channel:
        ufreqs = unique(freqs)
        for uf in ufreqs:
            I = (freqs==uf)
            R = ratios[I]
            L = locations[I]
            a, b = polyfit(L, R, 1)
            if scale_bands==0.5:
                R = -(R-b)
            else:
                R = ((R-b)/a)/amax(abs(locations))
            ratios[I] = -R
    banded_ratios = []
    banded_locations = []
    banded_freqs = sorted(unique(freqs))
    for uf in banded_freqs:
        I = (freqs==uf)
        banded_ratios.append(ratios[I])
        banded_locations.append(locations[I])
    return locations, ratios, freqs, banded_locations, banded_ratios, banded_freqs

def smoothplot(X, Y, F, w, xu=1):
    for x, y, f in zip(X, Y, F):
        c = (f-amin(F))/(amax(F)-amin(F))
        c = colour_gradient(c)
        if use_smooth:
            x, y = smoothed_scatter(x/xu, y, w)
            plot(x, y, '-', color=c, lw=2)
        else:
            #plot(x/xu, y, '.', color=c, ms=4, mew=0, mec=(1,1,1,0))
            scatter(x/xu, y, color=c, s=2)

if show_tuning_curves:
    training_filters, testing_filters = the_filters['whitenoise']
    use_ideal_responses = True
    analysis = get_analysis_from_namespace()
    locations, R, cf = get_sample_tuning_curves(num_frequency_bands)
    subplot(2, 3, 1)
    smoothplot([locations]*len(cf), R.T, cf, 0.05, xu=msecond)
    xlabel('ITD (ms)')
    ylabel('Firing rate (Hz)')
    axis('tight')
    xlim(amin(locations/ms), amax(locations/ms))

if show_two_channel_examples or show_cross_freq_examples:
    training_filters, testing_filters = the_filters['whitenoise']
    use_ideal_responses = False
    analysis = get_analysis_from_namespace()
    
    if show_two_channel_examples:
        locations, ratios, freqs, banded_locations, banded_ratios, banded_freqs = get_ratio_plot_data(
                                                       False,
                                                       num_frequency_bands,
                                                       scale_bands)
        subplot(2, 3, 2)
        smoothplot(banded_locations, banded_ratios, banded_freqs, 0.05, xu=msecond)
        xlabel('ITD (ms)')
        ylabel('Hemispheric difference')
        axis('tight')
        xlim(amin(locations/ms), amax(locations/ms))
    
    if show_cross_freq_examples:
        locations, ratios, freqs, banded_locations, banded_ratios, banded_freqs = get_ratio_plot_data(
                                                       True,
                                                       num_frequency_bands,
                                                       scale_bands)
        subplot(2, 3, 3)
        smoothplot(banded_locations, [1000*r for r in banded_ratios],
                   banded_freqs, 0.05, xu=msecond)
        xlabel('ITD (ms)')
        ylabel(r'Frequency-weighted difference')
        axis('tight')
        xlim(amin(locations/ms), amax(locations/ms))
        ylim(amin(1000*ratios), amax(1000*ratios))

from models.mcalpine_guinea_pig import *

#if use_low_frequencies:
#    limit_frequency = (0*Hz, 900*Hz)

if scale_bands:
    extra_twochan = 'Two channel banded'
    #extra_twochan = 'Two channel banded alt 10'
else:
#    extra_twochan = 'Two channel cross-frequency'
    extra_twochan = 'Two channel cross-frequency alt'

estimator_types = tuple((all_estimator_types[n], n) for n in [
                        'Two channel', 'Pattern match banded',
                        extra_twochan, 'Pattern match',
                        ])
fill_between_pairs = [(0, 2), (1, 3)]

if show_frequency_plot:
    training_filters, testing_filters = the_filters['whitenoise']
    use_ideal_responses = False
    analysis = get_analysis_from_namespace()
    subplot(2, 3, 4)
    pwithfreq = fig_performance_with_frequency(analysis, estimator_types,
                                   num_frequency_pools, frequency_pool_width,
                                   num_frequency_points,
                                   do_figure=False,
                                   do_suptitle=False,
                                   show_mean_error=True,
                                   show_bias=False,
                                   show_legend=False,
                                   show_freq_bars=False,
                                   formatting=formatting,
                                   frequency_unit=kHz)
    ylim(ymax=100)
    
if show_bandpass_plot:
    ymax = 0
    limit_frequency = the_limit_frequency
    for i, namekey in enumerate([['Two channel', extra_twochan],
                                 ['Pattern match', 'Pattern match banded']]):
        estimator_types = tuple((all_estimator_types[n], n) for n in namekey)
        fill_between_pairs = [(0, 1)]
        acousticnoisemodel = IndependentWhiteAcousticNoise((level, level))
        extraname['acousticnoisemodel'] = acousticnoisemodel
        training_filters, testing_filters = the_filters[bandpass_type]
        use_ideal_responses = False
        analysis = get_analysis_from_namespace()
        subplot(2, 2, 3+i)
        ymine, ymaxe, yminb, ymaxb = fig_performance_with_stimuli_tonefreq(
                analysis, estimator_types,
                tonefreqkey=bandpass_key,
                do_figure=False, do_suptitle=False,
                show_mean_error=True,
                show_bias=False,
                show_legend=False,
                show_ylabel=(i==0),
                formatting=formatting,
                frequency_unit=kHz,
                fill_between_alpha=fill_between_alpha,
                fill_between_pairs=fill_between_pairs,
                )
        if ymaxe>ymax:
            ymax = ymaxe
        for name, (fmid, me, se, fwid) in pwithfreq.items():
            form = formatting[name].copy()
            form['ls'] = ':'
            plot(fmid, me/usecond, **form)
        axis('tight')
        ylim(0)
        xlim(xmax=limit_frequency[1]/kHz)
        ylim(ymax=100)

#    ymax = 0
#    for i, use_low_frequencies in enumerate([True, False]):
#        if use_low_frequencies:
#            limit_frequency = the_limit_frequency
#        acousticnoisemodel = IndependentWhiteAcousticNoise((level, level))
#        extraname['acousticnoisemodel'] = acousticnoisemodel
#        training_filters, testing_filters = the_filters[bandpass_type]
#        use_ideal_responses = False
#        analysis = get_analysis_from_namespace()
#        subplot(2, 2, 3+i)
#        ymine, ymaxe, yminb, ymaxb = fig_performance_with_stimuli_tonefreq(
#                analysis, estimator_types,
#                tonefreqkey=bandpass_key,
#                do_figure=False, do_suptitle=False,
#                show_mean_error=True,
#                show_bias=False,
#                show_legend=False,
#                show_ylabel=(i==0),
#                formatting=formatting,
#                frequency_unit=kHz,
#                fill_between_alpha=fill_between_alpha,
#                fill_between_pairs=fill_between_pairs,
#                )
#        if ymaxe>ymax:
#            ymax = ymaxe
#        for name, (fmid, me, se, fwid) in pwithfreq.items():
#            form = formatting[name].copy()
#            form['ls'] = ':'
#            plot(fmid, me/usecond, **form)
#        axis('tight')
#        ylim(0)
#        if use_low_frequencies:
#            xlim(xmax=limit_frequency[1]/kHz)
#        ylim(ymax=100)
        
#if use_low_frequencies:
#    limit_frequency = (0*Hz, 900*Hz)        
        
#    subplot(2, 3, 6)
#    fig_performance_with_stimuli_tonefreq(analysis, estimator_types,
#            tonefreqkey=bandpass_key,
#            do_figure=False, do_suptitle=False,
#            show_mean_error=False,
#            show_bias=True,
#            show_legend=False,
#            formatting=formatting,
#            frequency_unit=kHz,
#            )
#    #axis('tight')
#    if use_low_frequencies:
#        xlim(xmax=limit_frequency[1]/kHz)
#        ylim(ymax=75)

for n, figlet in [((2, 3, 1), 'A'),
                  ((2, 2, 3), 'B'),
                  #((2, 3, 5), 'C'),
                  ]:
    subplot(*n)
    text(-0.3, 1.1, figlet, size=16, transform=gca().transAxes)

subplot(2, 2, 3)
xlim(xmax=limit_frequency[1]/kHz)
subplot(2, 2, 4)
xlim(xmax=limit_frequency[1]/kHz)
    
show()
