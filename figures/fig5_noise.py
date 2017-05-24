from base import *
from acousticnoise import *
from analyse_model import *
from colours import *

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.15
rcParams['figure.subplot.right'] = .97
rcParams['figure.subplot.bottom'] = .07
rcParams['figure.subplot.top'] = .95
rcParams['figure.subplot.wspace'] = 0.35
rcParams['figure.subplot.hspace'] = 0.25
# figsize is width x height in inches
figure(figsize=(6, 6), dpi=100)

show_protocol = True
show_performance = True
show_guinea_pig = True
show_cat = True
#limit_frequency = (0*Hz, 900*Hz)
#limit_frequency = (0*Hz, 600*Hz)
limit_frequency = (0*Hz, 1200*Hz)
log_scale_error = False
log_scale_bias = False

SNR = 10*dB
duration = 2*ms
delay = 300*usecond
freq = 1000*Hz

formatting = {
    'Smoothed peak':{'color':estimator_colours['peak'], 'ls':'-', 'lw':2},
#    'Two channel':{'color':estimator_colours['twochannel'], 'ls':'--', 'lw':2},
    'Two channel cross-frequency':{'color':estimator_colours['twochannel'], 'ls':'-', 'lw':2},
#    'Two channel banded':{'color':'m'},
    'Pattern match':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':2},
    'Pattern match banded':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':2},
    }

def set_all_estimator_types():
    global estimator_types
    estimator_types = (
        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)),
                       'Smoothed peak'),
#        (MakeEstimator(TwoChannel, PolyClosest(6), difference=True,
#                                   itdmax_extend=itdmax_extend), 'Two channel'),
#        (MakeEstimator(TwoChannelBanded), 'Two channel banded'),
        (MakeEstimator(TwoChannelCrossFrequencyAlt, PolyClosest(6),
                       itdmax_extend=itdmax_extend),
                'Two channel cross-frequency'),
#        (MakeEstimator(PatternMatch), 'Pattern match'),
        (MakeEstimator(PatternMatch, normalise_banded=40),
                'Pattern match banded'),
        )

use_ideal_responses = False # Removes all response noise from the results
num_shuffles = 25
training_size = 400
testing_size = 800

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

the_models = []
if show_guinea_pig:
    the_models.append(('mcalpine_guinea_pig', 2))
if show_cat:
    the_models.append(('joris_cat', 4))

def blankplot():
    xticks([])
    yticks([])
    gca().set_frame_on(False)
    #ylim(-1.2, 1.2)
    
if show_protocol:
    sound = tone(freq, duration+delay)
    sound_left = sound[:duration]
    sound_right = sound[delay:]
    sound = sound[:duration]
    noise_left = whitenoise(duration)*gain(-SNR)
    noise_right = whitenoise(duration)*gain(-SNR)
    combined_left = sound_left+noise_left
    combined_right = sound_right+noise_right
    ymax = max([amax(abs(sound)),
                amax(abs(noise_left)),
                amax(abs(noise_right)),
                amax(abs(combined_left)),
                amax(abs(combined_right))])
    
    subplot(6,4,1) # original signal
    plot(sound.times, sound, '-k')
    title('Source')
    ylim(-ymax, ymax)
    blankplot()
    subplot(12,4,13)
    title('TODO:\nAlign and arrows')
    blankplot()
    
    subplot(6,4,2) # ITD left
    plot(sound_left.times, sound_left, '-k')
    ylim(-ymax, ymax)
    title('With ITD')
    blankplot()
    subplot(6,4,6) # ITD right
    plot(sound_right.times, sound_right, '-k')
    ylim(-ymax, ymax)
    blankplot()
    
    subplot(6,4,3) # noise left
    plot(noise_left.times, noise_left, '-k')
    ylim(-ymax, ymax)
    title('Noise')
    blankplot()
    subplot(6,4,7) # noise right
    plot(noise_right.times, noise_right, '-k')
    ylim(-ymax, ymax)
    blankplot()
    
    subplot(6,4,4) # combined left
    plot(combined_left.times, combined_left, '-k')
    title('Signal at ears')
    blankplot()
    subplot(6,4,8) # combined right
    plot(combined_right.times, combined_right, '-k')
    blankplot()

if show_performance:
    for modelname, offset in the_models:
        exec 'from models.%s import *' % modelname
        set_all_estimator_types()
        acousticnoisemodel = IndependentWhiteAcousticNoise(tuple(level-sn for sn in sn_level_range[::-1]))
        extraname['acousticnoisemodel'] = acousticnoisemodel
        analysis = get_analysis_from_namespace()
        subplot(3, 2, 1+offset) # mean error
        fig_performance_with_acoustic_noise(analysis, estimator_types,
                                            do_figure=False,
                                            do_suptitle=False,
                                            show_mean_error=True,
                                            show_bias=False,
                                            show_legend=False,
                                            show_xlabel=(offset!=2),
                                            show_ylabel=True,
                                            formatting=formatting,
                                            log_y_axis=log_scale_error,
                                            )
#        axis('tight')
        ylabel(animal_name+'\n'+gca().get_ylabel())
        if offset==2:
            title('Mean error')
        ylim(ymin=0)
        subplot(3, 2, 2+offset) # bias
        fig_performance_with_acoustic_noise(analysis, estimator_types,
                                            do_figure=False,
                                            do_suptitle=False,
                                            show_mean_error=False,
                                            show_bias=True,
                                            show_legend=False,
                                            show_xlabel=(offset!=2),
                                            show_ylabel=True,
                                            formatting=formatting,
                                            log_y_axis=log_scale_bias,
                                            )
        if offset==2:
            title('Bias')
            ylim(ymin=-10)
        else:
            ylim(ymin=-5)

subplot(6, 4, 1)
text(-0.3, 1.1, 'A', size=16, transform=gca().transAxes)
subplot(3, 2, 3)
text(-0.3, 1.1, 'B', size=16, transform=gca().transAxes)

show()
