from base import *
from binauraldisplay import *
from models.mcalpine_guinea_pig import *

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.13
rcParams['figure.subplot.right'] = .98
rcParams['figure.subplot.bottom'] = .12
rcParams['figure.subplot.top'] = .95
rcParams['figure.subplot.wspace'] = 0.3
rcParams['figure.subplot.hspace'] = 0.5
# figsize is width x height in inches
figure(figsize=(7, 12), dpi=70)

show_best_delay_distribution = True
show_model = True
show_example_response_curve = True
show_example_model_responses = True

markersize = 4
#mincol, maxcol = 0.2, 0.8
mincol, maxcol = 0.0, 1.0
use_bestphase = False
tightaxis = True
formatting = None

bginterp = 'cells'
#min_marker_size, max_marker_size = 2.0, 20.0
#min_marker_size, max_marker_size = 30.0, 31.0
min_marker_size, max_marker_size = 1, 2

example_response_curve_freq = 1*kHz
example_response_curve_bd = 100*usecond
example_response_curve_itd_range = (-1*ms, 1*ms, 100)

############## BEST DELAY DISTRIBUTION

if show_best_delay_distribution:
    subplot(4, 2, 1)
    best_delay_distribution(cf, bd, usephase=use_bestphase, tightaxis=tightaxis,
                            markersize=markersize)

if show_model:
    start = .05*second
    duration = 4*ms
    subplot(4, 2, 3)
    sound = whitenoise(.1*second)
    s = sound[start:][:duration]
    plot(s.times, s, '-k')
    xticks([])
    yticks([])
    gca().set_frame_on(False)
    subplot(4, 2, 4)
    sound = Sound(Gammatone(sound, [1*kHz]).process())[start:][:duration]
    plot(sound.times, sound, '-k')
    plot(sound.times+300*usecond, sound, '--k')
    xticks([])
    yticks([])
    gca().set_frame_on(False)

    
if show_example_model_responses:
    i = 5
    stimuli = {
        'White noise':whitenoise_stimuli(space).next()[0],
        'Natural sound':natural_sounds_stimuli(space).next()[0]#, duration=5*second)
        }
    acousticnoisemodel = NoAcousticNoise()
    fm = FixedModel(space, acousticnoisemodel, binauralmodel,
                    responsenoisemodel, binauraldistribution,
                    generate_extended_datamanager_name(basename, extraname),
                    gammatone_params=gammatone_params,
                    compression=compression,
                    rms_normalisation=rms_normalisation,
                    )
    @figcache
    def example_response(itd, soundtype):
        sound = stimuli[soundtype]
        response, _, _, _ = fm.apply(sound, itd)
        return response
    # TODO switch the order round?
    for u, itd in enumerate([-250*usecond, 100*usecond]):
        for v, soundtype in enumerate(['White noise', 'Natural sound']):
            subplot(4, 2, i)
            response = example_response(itd, soundtype)
            if v==0:
                show_ylabel = 'ITD = %g $\mathrm{\mu}$s'%int(itd/usecond)
            else:
                show_ylabel = False  
            if u==1:
                show_xlabel = soundtype
            else:
                show_xlabel = False
            display_response(cf, bd, response, itd=itd,
                             use_bestphase=use_bestphase,
                             tightaxis=tightaxis,
                             min_marker_size=min_marker_size,
                             max_marker_size=max_marker_size,
                             mincol=mincol, maxcol=maxcol,
                             show_xlabel=show_xlabel,
                             show_ylabel=show_ylabel,
                             bginterp=bginterp,
                             formatting=formatting,
                             )
            i += 1

if show_example_response_curve:
    itdmin, itdmax, itdN = example_response_curve_itd_range
    subplot(4, 2, 2)
    cf = example_response_curve_freq*ones(itdN)
    bd = linspace(*example_response_curve_itd_range)
    @figcache
    def get_example_response_curve():
        responsenoisemodel = NoResponseNoise()
        extraname['responsenoisemodel'] = responsenoisemodel
        dl, dr = delays_from_bd(bd)
        binauraldistribution = BinauralDistribution(zip(cf, dl, dr), cfrepeat=cfrepeat)
        acousticnoisemodel = NoAcousticNoise()
        gammatone_params = {'ear_Q':beta*(hstack((cf, cf))/kHz)**alpha,
                            'min_bw':0.0}
        fm = FixedModel(space, acousticnoisemodel, binauralmodel,
                        responsenoisemodel, binauraldistribution,
                        generate_extended_datamanager_name(basename, extraname),
                        gammatone_params=gammatone_params,
                        compression=compression,
                        rms_normalisation=rms_normalisation,
                        )
        sound = whitenoise(.1*second).atlevel(70*dB)[:.3*second]
        R, _, _, _ = fm.apply(sound, example_response_curve_bd)
        return R, sound.duration
    R, duration = get_example_response_curve()
    if use_bestphase:
        plot(2*pi*example_response_curve_freq*bd, R/duration, '-k')
        axvline(x=2*pi*example_response_curve_freq*example_response_curve_bd, ls='--', color='k')
        xlabel(r'IPD (rad)')
    else:
        plot(bd/ms, R/duration, '-k')
        axvline(x=example_response_curve_bd/msecond, ls='--', color='k')
        xlabel(r'ITD (ms)')
    ylabel('Firing rate (Hz)')
    if tightaxis: axis('tight')

for fignum in xrange(3):
    subplot(4, 2, fignum*2+1)
    text(-0.3, 1.1, chr(ord('A')+fignum), size=14, transform=gca().transAxes) 

show()
