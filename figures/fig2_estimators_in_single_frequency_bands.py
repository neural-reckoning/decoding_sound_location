from base import *
from frequency import *
from binauraldisplay import *
from analyse_model import *
from colours import *

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.08
rcParams['figure.subplot.right'] = .98
rcParams['figure.subplot.bottom'] = .08
rcParams['figure.subplot.top'] = .95
rcParams['figure.subplot.wspace'] = 0.4
rcParams['figure.subplot.hspace'] = 0.5
# figsize is width x height in inches
figure(figsize=(12, 9), dpi=70)

show_example_twochan_peak = True
show_twochan_ratios = True
show_guinea_pig = True
show_cat = True
show_example_pattern_match = True

use_bestphase = False
show_two_chan_ratio_fit = True
show_two_chan_ratio_points = False
markersize = 4
tightaxis = True
two_chan_ratio_point_colour = (.5, .5, .5)
pm_max_marker_size = 200
pm_axis = [0.4, 0.55, -0.6, 0.6]
example_itd = 200*usecond
example_centre = .7*kHz
example_max_marker_size = 200
example_pointcol = 0.4
example_min_ratio_col = 0.7
cpos = (1.0, 0.9, 0.0)
cneg = (1.0, 0.6, 0.0)

formatting = {
    'Peak coding':{'color':estimator_colours['peak'], 'ls':'--', 'lw':2},
    'Smoothed peak':{'color':estimator_colours['peak'], 'ls':'-', 'lw':2},
    'Two channel':{'color':estimator_colours['twochannel'], 'ls':'-', 'lw':2},
    'Pattern match':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':2},
    'Pattern match banded':{'color':estimator_colours['patternmatch'], 'ls':'--', 'lw':2},
    }

def set_all_estimator_types():
    global all_estimator_types
    all_estimator_types = {
        'Peak coding':MakeEstimator(Jeffress, NoisyMax()),
        'Smoothed peak':MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)),
        'Two channel':MakeEstimator(TwoChannel, PolyClosest(6), difference=True,
                                    itdmax_extend=itdmax_extend),
        'Pattern match':MakeEstimator(PatternMatch),
        'Pattern match banded':MakeEstimator(PatternMatch, normalise_banded=40),
        }

use_ideal_responses = False # Removes all response noise from the results
num_shuffles = 25
training_size = 400
testing_size = 800
num_frequency_pools = 3
frequency_pool_width = 40
num_frequency_points = 10
training_filters = (
    'type=="whitenoise"',
    'i<training_size',
    )
testing_filters = (
    'type=="whitenoise"',
    'i<testing_size',
    )

from models.mcalpine_guinea_pig import *
set_all_estimator_types()
acousticnoisemodel = NoAcousticNoise()
extraname['acousticnoisemodel'] = acousticnoisemodel
analysis = get_analysis_from_namespace()

#    text(-0.3, 1.1, figlet, size=14, transform=gca().transAxes)

# LAYOUT
ax_width = ax_height = 0.22
inset_width = inset_height = 0.1
inset_offset = -0.02
vspacing = 0.08
ratio_inset = 0.01
ratio_width = ax_width/2
ratio_height = ax_height/8

L1 = 0.25-0.75*ax_width
L2 = 0.5-0.5*ax_width
L3 = 0.75-0.25*ax_width
B1 = 1.0-vspacing-ax_height

#rh = 0.5-2*vspacing-0.5*ax_height-inset_height
rh = 0.5-2*vspacing-0.5*ax_height
L1i = L1+ax_width-inset_width-inset_offset
B2 = vspacing
#B3 = 2*vspacing+rh+inset_height
B3 = 2*vspacing+rh

ax_example = axes([L1, B1, ax_width, ax_height])
ax_ex_est = axes([L2, B1, ax_width, ax_height])
ax_twochan_ratio = axes([L2+ax_width-ratio_width-ratio_inset, B1+ratio_inset, ratio_width, ratio_height])
xticks([]); yticks([]);
ax_ratio = axes([L3, B1, ax_width, ax_height])

ax_pm_data = axes([L2, 0.5-ax_height-vspacing, ax_width, ax_height])
ax_pm_pat1 = axes([L3, 0.5-ax_height/2-vspacing/2, ax_width, ax_height])
ax_pm_pat2 = axes([L3, 0.5-1.5*ax_height-1.5*vspacing, ax_width, ax_height])

ax_res_gpig = axes([L1, B3, ax_width, rh])
ax_res_gpig_inset = axes([L1i, B3+rh-inset_height-inset_offset, inset_width, inset_height])

ax_res_cat = axes([L1, B2, ax_width, rh])
ax_res_cat_inset = axes([L1i, B2+rh-inset_height-inset_offset, inset_width, inset_height])

for ax, figlet in [(ax_example, 'A'),
                   (ax_res_gpig, 'B'),
                   (ax_pm_data, 'C')]:
    axes(ax)
    text(-0.3, 1.1, figlet, size=16, transform=gca().transAxes)

if True:

    if show_example_twochan_peak:
        @figcache
        def get_example():
            fm = analysis.fm_base
            cf = fm.centre_frequencies[:fm.cfN]
            closest_to_centre = argmin(abs(cf-example_centre))
            I = argsort(abs(cf-example_centre))
            curI = I[:frequency_pool_width]
            sound = whitenoise(100*ms)
            response, _, _, _ = fm.apply(sound, example_itd)
            return cf, bd, curI, response
        cf, bd, curI, response = get_example()
        # twochan and peak
        #subplot(3, 4, 1)
        axes(ax_example)
        cfmin = amin(cf[curI]/kHz)
        cfmax = amax(cf[curI]/kHz)
        bdmin = amin(bd[curI]/msecond)
        bdmax = amax(bd[curI]/msecond)
#        rtot = rpos+rneg
#        rpos = rpos/rtot
#        rneg = rneg/rtot
#        cpos = (example_min_ratio_col,)*3
#        cneg = (1-0.25*(1-example_min_ratio_col),)*3
#        poly1 = fill([cfmin, cfmax, cfmax, cfmin], [0, 0, bdmax, bdmax], color=cpos)
#        poly2 = fill([cfmin, cfmax, cfmax, cfmin], [0, 0, bdmin, bdmin], color=cneg)
#        for p in poly1+poly2:
#            p.set_zorder(0.9)
        display_response(cf, bd, response, itd=example_itd,
                         max_marker_size=example_max_marker_size,
                         mincol=example_pointcol, maxcol=example_pointcol,
                         coloured=False, colour=(0.2, 0.4, 0.9),
                         bginterp=None,
                         )
        axis([cfmin, cfmax, bdmin, bdmax])
        i = curI[argmax(response[curI])]
        plot([cf[i]/kHz], [bd[i]/msecond], 'o', color=(0.4, 0.6, 1.0),
             ms=sqrt(example_max_marker_size))
        plot([cf[i]/kHz], [bd[i]/msecond], 'x', color='k',
             ms=sqrt(example_max_marker_size))
        # smoothed peak    
        #subplot(3, 4, 2)
        rpos = mean(1.0*response[curI][bd[curI]>0])
        rneg = mean(1.0*response[curI][bd[curI]<0])
        axes(ax_ex_est)
        x = bd[curI]
        y = response[curI]
        poly1 = fill([0, rpos, rpos, 0], [0, 0, bdmax, bdmax], color=cpos)
        poly2 = fill([0, rneg, rneg, 0], [0, 0, bdmin, bdmin], color=cneg)
        for p in poly1+poly2:
            p.set_zorder(0.9)
        plot(y, x/msecond, '.', #color=(example_pointcol,)*3,
             ms=markersize*3,
             color=(0.2, 0.4, 0.9),
             mew=0)
        xs = linspace(amin(x), amax(x), 100)
        ys = gaussian_smooth(x, y, 0.15*space.itd_max, xs)
        plot(ys, xs/msecond, '-k', lw=2)
        i = argmax(y)
        plot([y[i]], [x[i]/msecond], 'ok', ms=3*markersize, mfc='none', mew=2)
        i = argmax(ys)
        plot([ys[i]], [xs[i]/msecond], 'ok', ms=3*markersize, mfc='none', mew=2)
        axhline(example_itd/msecond, ls='--', color='k', lw=2)
        ylabel('Best delay (ms)')
        xlabel('Spike count')
        ylim(ymin=amin(x)/msecond, ymax=amax(x)/msecond)
        axes(ax_twochan_ratio)
        r = rpos/(rpos+rneg)
        p1=fill([0, r, r, 0], [0, 0, 1, 1], color=cpos)
        p2=fill([r, 1, 1, r], [0, 0, 1, 1], color=cneg)
        for p in p1+p2:
            p.set_zorder(0.9)
        axvline(r, ls='-', color='k', lw=1)
        axis([0, 1, 0, 1])
        title('Ratio')
    
    def show_animal(subplot_offset):
        # performance for peak, 2chan, smoothed peak
        #subplot(3, 4, 3+subplot_offset)
        if subplot_offset==0:
            axes(ax_res_gpig)
        else:
            axes(ax_res_cat)
        title(animal_name)
        estimator_types = tuple((all_estimator_types[n], n) for n in [
                                'Peak coding', 'Smoothed peak', 'Two channel', 'Pattern match',
                                ])
        fig_performance_with_frequency(analysis, estimator_types,
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
        axis('tight')
        if subplot_offset==0:
            ylim(0, 200)
        else:
            ylim(0, 350)
        # BD distribution
        #subplot(3, 4, 7+subplot_offset)
        if subplot_offset==0:
            axes(ax_res_gpig_inset)
        else:
            axes(ax_res_cat_inset)
        best_delay_distribution(cf, bd, usephase=use_bestphase, tightaxis=tightaxis,
                                markersize=1.5)
    
    if show_guinea_pig:
        show_animal(0)
    
    if show_twochan_ratios:
        estimator_type = all_estimator_types['Two channel']
        def get_trained_estimator(centre, width):
            fm_orig = analysis.fm_orig
            cf = fm_orig.centre_frequencies[:fm_orig.cfN]
            closest_to_centre = argmin(abs(cf-centre))
            I = argsort(abs(cf-centre))
            curI = I[:width]
            fm = SubsetFixedModel(fm_orig, curI)
            analysis.fm = fm
            estimator = estimator_type(fm)
            training_responses, testing_responses = analysis.traintest_responses(None)
            estimator.train(responses=training_responses)
            return fm, estimator
        @figcache
        def get_two_chan_ratio_data(centre, width):
            fm, estimator = get_trained_estimator(centre, width)
            f = estimator.fitter
            return estimator.locations, estimator.values, f.x, f.y
        def do_two_chan_ratio_plot(centre, width, dotitle=False, colour='k'):
            l, v, x, y = get_two_chan_ratio_data(centre, width)
            l = l/usecond
            if show_two_chan_ratio_points:
                plot(l, v, '.', color=colour, ms=4, mew=0)
            if show_two_chan_ratio_fit:
                I = argsort(l)
                l = l[I]
                v = v[I]
                plot(l, gaussian_smooth(l, v, 10), color=colour, lw=2)
                #I = abs(x/usecond)<300
                #plot(x[I]/usecond, y[I], '-k')
            if dotitle: title('CF near %s'%centre)
            xlabel('ITD ($\mu$s)')
            ylabel('Ratio')
        # low frequency
        #subplot(3, 4, 5)
        axes(ax_ratio)
        do_two_chan_ratio_plot(example_centre, frequency_pool_width,
                               colour=(0.2, 0.4, 0.9))
        # high frequency (ambiguous)
        #subplot(3, 4, 6)
        do_two_chan_ratio_plot(1.3*kHz, frequency_pool_width,
                               colour=(0.7, 0, 0.5),
                               )#, dotitle=True)
        axhline(-0.85, ls='--', color='k', lw=2)
    
    if show_cat or show_example_pattern_match:
        from models.joris_cat import *
        set_all_estimator_types()
        acousticnoisemodel = NoAcousticNoise()
        extraname['acousticnoisemodel'] = acousticnoisemodel
        analysis = get_analysis_from_namespace()
        
        if show_cat:
            show_animal(1)
        
        if show_example_pattern_match:
            fm = analysis.fm_base
            trueitd = 400*usecond
            falseitd = -400*usecond
            @figcache
            def get_example_pattern_match(trueitd, falseitd):
                response, _, _, _ = fm.apply(whitenoise(100*ms), trueitd)
                _, _, meta, _ = fm.apply(whitenoise(100*ms), trueitd)
                truepattern = meta['noiseless_response']
                _, _, meta, _ = fm.apply(whitenoise(100*ms), falseitd)
                falsepattern = meta['noiseless_response']
                return response, truepattern, falsepattern
            response, truepattern, falsepattern = get_example_pattern_match(trueitd, falseitd)
            axes(ax_pm_data)
            col = (0.5, 0.5, 0.5)
            display_response(cf, bd, (1.0*response)/amax(response),
                             max_marker_size=pm_max_marker_size,
                             coloured=False, colour=col,
                             bginterp=None,
                             )
            axhline(trueitd/msecond, ls='--', color='k', lw=2)
            axis(pm_axis)
            title('Response')
            axes(ax_pm_pat1)
            display_response(cf, bd, falsepattern/amax(falsepattern),
                             max_marker_size=pm_max_marker_size,
                             coloured=False, colour=col,
                             bginterp=None,
                             )
            axhline(falseitd/msecond, ls='--', color='k', lw=2)
            axis(pm_axis)
            title('Pattern A')
            axes(ax_pm_pat2)
            display_response(cf, bd, truepattern/amax(truepattern),
                             max_marker_size=pm_max_marker_size,
                             coloured=False, colour=col,
                             bginterp=None,
                             )
            axis(pm_axis)
            axhline(trueitd/msecond, ls='--', color='k', lw=2)
            title('Pattern B')

for ax in [ax_res_cat_inset, ax_res_gpig_inset,
           ax_pm_data, ax_pm_pat1, ax_pm_pat2]:
    axes(ax)
    xticks([])
    yticks([])
    xlabel('')
    ylabel('')

show()
