import os
os.environ['SAMPLERATE'] = '192000' # guinea pig
#os.environ['SAMPLERATE'] = '97656' # cat
from base import *
from analyse_model import *
from cells import *
from binauraldisplay import *
from colours import *

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.07
rcParams['figure.subplot.right'] = .97
rcParams['figure.subplot.bottom'] = .1
rcParams['figure.subplot.top'] = .92
rcParams['figure.subplot.wspace'] = 0.4
rcParams['figure.subplot.hspace'] = 0.6
# figsize is width x height in inches
figure(figsize=(10, 5), dpi=80)

show_guinea_pig_photo = True
show_itd_freq_plots = True
show_binaural_displays = True
show_results = True
show_guinea_pig = True
show_cat = True
limit_frequency = (0*Hz, 1200*Hz)
show_itd_freq_range = (250*Hz, 1.5*kHz)
azims_to_show = (90, 60, 30)
#azimcols = [(0, 0, 1), (0.5, 0, 1), (1, 0, 1)]
azimcols = [(.4, .4, 1),
            (0, 0, 1),
            (0, 0, .5),
            ]
Nticks_skip = 2
min_marker_size, max_marker_size = 2.0, 20.0
#min_marker_size, max_marker_size = 0.0, 0.0
#mincol, maxcol = 0.2, 0.8
mincol, maxcol = 0.0, 1.0
bginterp = 'cells'

test_sound_duration = 1*second
test_sound_ideal = True
test_sound_no_limit_frequency = True

the_models = []
if show_guinea_pig:
    the_models.append(('ircam_mcalpine_guinea_pig', 0))
if show_cat:
    the_models.append(('joris_tollin_cat', 1))

def blankplot():
    xticks([])
    yticks([])
    gca().set_frame_on(False)

def set_all_estimator_types():
    global estimator_types
    azim = space.hrtfset.coordinates['azim']
    itd = space.itd_max*sin(azim*pi/180)
    estimator_types = (
#        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)),
#                       'Smoothed peak'),
#        (MakeEstimator(TwoChannel, PolyClosest(6, grid=itd), difference=True,
#                                   itdmax_extend=itdmax_extend), 'Two channel'),
#        (MakeEstimator(TwoChannelBanded), 'Two channel banded'),
#        (MakeEstimator(TwoChannelCrossFrequency, PolyClosest(6),
#                       itdmax_extend=itdmax_extend),
#                'Two channel cross-frequency'),
        (MakeEstimator(TwoChannelCrossFrequencyAlt,
                       #PolyClosest(6, grid=itd),
                       Closest(),
                       itdmax_extend=itdmax_extend),
                'Two channel cross-frequency alt'),
#        (MakeEstimator(PatternMatch), 'Pattern match'),
        (MakeEstimator(PatternMatch, normalise_banded=40),
                'Pattern match banded'),
        )

use_ideal_responses = False # Removes all response noise from the results
num_shuffles = 25
training_size = 400
testing_size = 800
#num_shuffles = 5
#training_size = 400
#testing_size = 200

training_filters = (
    'type=="whitenoise"',
    'i<training_size',
    )
testing_filters = (
    'type=="whitenoise"',
    'i<testing_size',
    )

formatting = {
    'Smoothed peak':{'color':estimator_colours['peak'], 'ls':'-', 'lw':2},
    'Two channel':{'color':estimator_colours['twochannel'], 'lw':2},
#    'Two channel cross-frequency':{'color':'g'},
    'Two channel cross-frequency alt':{'color':estimator_colours['twochannel'],
                                       'lw':2, 'ls':'-'},
#    'Two channel banded':{'color':'m'},
    'Pattern match':{'color':estimator_colours['patternmatch'],
                     'lw':2, 'ls':'-'},
    'Pattern match banded':{'color':estimator_colours['patternmatch'],
                     'lw':2, 'ls':'-'},
    }

the_animal_name = {
    'joris_tollin_cat':'Cat',
    'ircam_mcalpine_guinea_pig':'Guinea pig'
    }

@figcache
def get_plot_data(modelname, show_itd_freq_plots, show_binaural_displays,
                  show_results, #limit_frequency,
                  azims_to_get):
    global acousticnoisemodel, limit_frequency, binauralmodel
    orig_limit_frequency = limit_frequency
    if test_sound_no_limit_frequency:
        limit_frequency = False
    results = {}
    #exec 'from models.%s import *' % modelname in globals()
    execfile('../models/'+modelname+'.py', globals())
    set_all_estimator_types()
    acousticnoisemodel = NoAcousticNoise()
    extraname['acousticnoisemodel'] = acousticnoisemodel
    analysis = get_analysis_from_namespace()
    fmbase = analysis.fm_base
    if show_itd_freq_plots:
        sspace = space.subspace(lambda azim: azim in azims_to_get)
        all_bitds, all_cfind = sspace.get_bitds(fmbase,
                                                itdmax=1.0*ms,
                                                continuous=True,
                                                interpolation=True,
                                                )
        azims = sspace.hrtfset.coordinates['azim']
        plotpairs = []
        for azim, bitd, cfind in zip(azims, all_bitds, all_cfind):
            #if azim==0 or azim==45 or azim==90:
            if azim in azims_to_get:
                plotpairs.append((cf, bitd, azim))
            if int(azim)==azims_to_show[0]:
                show_bitd = bitd
        results['itd'] = plotpairs
    if show_binaural_displays:
        sound = whitenoise(test_sound_duration)
        itd = -sin(azims_to_show[0]*pi/180)*space.itd_max
        if modelname=='joris_tollin_cat':
            mulfac = -1
        else:
            mulfac = 1
        itd *= mulfac
        response, _, meta, _ = fmbase.apply(sound, itd)
        if test_sound_ideal:
            response = meta['noiseless_response']
        results['bd'] = (cf, bd, response, itd, -mulfac*show_bitd)
    if show_results:
#        try:
        limit_frequency = orig_limit_frequency
        #exec 'from models.%s import *' % modelname in globals()
        execfile('../models/'+modelname+'.py', globals())
        set_all_estimator_types()
        acousticnoisemodel = NoAcousticNoise()
        extraname['acousticnoisemodel'] = acousticnoisemodel
        analysis = get_analysis_from_namespace()
        plotcommands = fit_performance_with_cells(analysis, estimator_types,
            number_of_mso_neurons=number_of_mso_neurons[animal_name],
            do_figure=False, do_suptitle=False, do_legend=False,
            show_guesslevel=False,
            show_infinity=False,
            show_mso=False,
            show_fit=False,
            extra_label_spacing=3,
            #Nticks=Nticks,
            #itd_unit=ms,
            angles=True,
            exec_plotcommands=False,
            formatting=formatting,
            Nticks_skip=Nticks_skip,
            )
#        except ValueError:
#            plotcommands = ()
        results['performance'] = plotcommands
    return results

if show_guinea_pig_photo:
    subplot(2, 4, 1)
    I = imread('guineapig.jpg')
    imshow(I, origin='lower left')
    blankplot()
    space = IRCAMHRTFSpatial(subject=3012, azim_offset=-25,
                             adaptive_windowing=(512, 16, 256)
                             )    
    for i, azim in enumerate(azims_to_show):
        subplot(6, 4, 2+4*i)
        #azim = 45*(2-i)
        ir = space.hrtfset(azim=azim, elev=0).impulse_response[:1.5*ms]
        M = amax(abs(ir))
        plot(ir.times/ms, ir.left, color=(1, .5, 0))
        plot(ir.times/ms, ir.right, color=(.4, 0, .4))
        axis('tight')
        ylim(-M, M)
        yticks([])
        xt = [0, 0.5, 1.0, 1.5]
        #xt, _ = xticks()
        if i<2:
            xticks(xt, ['' for _ in xt])
        else:
            xticks(xt)
            xlabel('Time (ms)')
        ylabel('%d$^\circ$'%azim)

if show_itd_freq_plots or show_binaural_displays or show_results:
    for modelname, offset in the_models:
        
        thedata = get_plot_data(modelname, show_itd_freq_plots,
                                show_binaural_displays, show_results,
                                #limit_frequency,
                                azims_to_show)
    
        if show_itd_freq_plots:
            subplot(2, 4, 3+offset)
            plotpairs = thedata['itd']
            for (cf, bitd, azim), c in zip(plotpairs, azimcols):
                if the_animal_name[modelname]=='Cat':
                    plot(cf/kHz, bitd/usecond, label=str(azim), color=c, lw=2)
                else:
                    plot(cf/kHz, -bitd/usecond, label=str(azim), color=c, lw=2)
            axhline(0, ls='--', color='k')
            xlabel('CF (kHz)')
            ylabel('ITD ($\mu$s)')
            #axis('tight')
            rlow, rhigh = show_itd_freq_range
            if amin(cf)*Hz>rlow:
                rlow = amin(cf)*Hz
            if amax(cf)*Hz<rhigh:
                rhigh = amax(cf)*Hz
            xlim(rlow/kHz, rhigh/kHz)
            title(the_animal_name[modelname])
            
        if show_binaural_displays:
            subplot(2, 4, 5+offset)
            cf, bd, response, itd, show_itd = thedata['bd']
            display_response(cf, bd, response, #itd=-itd,
                             use_bestphase=False,#use_bestphase,
                             tightaxis=True,#tightaxis,
                             min_marker_size=min_marker_size,
                             max_marker_size=max_marker_size,
                             mincol=mincol, maxcol=maxcol,
                             show_xlabel=True,
                             show_ylabel=True,
                             bginterp=bginterp,
                             horiz_stretch=1.5,
                             )
            plot(cf/kHz, show_itd/msecond, '--k', lw=2)
            if limit_frequency is not False and not test_sound_no_limit_frequency:
                xlim(xmax=limit_frequency[1]/kHz)
            ylim(-1, 1)
            title(the_animal_name[modelname])
            
        if show_results:
            plotcommands = thedata['performance']
            subplot(2, 4, 7+offset)
            for f, args, kwds in plotcommands:
                f(*args, **kwds)
            title(the_animal_name[modelname])

for n, figlet in [(1, 'A'),
                  (3, 'B'),
                  (5, 'C'),
                  (7, 'D')]:
    subplot(2, 4, n)
    text(-0.3, 1.1, figlet, size=16, transform=gca().transAxes)

show()
