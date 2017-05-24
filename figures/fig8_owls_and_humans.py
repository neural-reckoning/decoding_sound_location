import os
#os.environ['SAMPLERATE'] = '96000' # owls

from base import *
from frequency import *
from cells import *
from stimulialpha import *
from tonefreq import *
from binauraldisplay import *
from analyse_model import *
from colours import *

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.15
rcParams['figure.subplot.right'] = .98
rcParams['figure.subplot.bottom'] = .15
rcParams['figure.subplot.top'] = .9
rcParams['figure.subplot.wspace'] = 0.4
rcParams['figure.subplot.hspace'] = 0.4
# figsize is width x height in inches
figure(figsize=(6, 5), dpi=80)

show_models = [
    'wagner_owl',
    'ircam_human_uniform_bipd',
    'ircam_human_mcalpine_bd',
    'ircam_human_tollin_bd',
    ]

lw = 2
formatting = {
    'Smoothed peak':{'color':estimator_colours['peak'], 'ls':'-', 'lw':lw},
    'Two channel':{'color':estimator_colours['twochannel'], 'ls':'-', 'lw':lw},
    'Two channel cross-frequency alt':{'color':estimator_colours['twochannel'], 'ls':'-', 'lw':lw},
    'Pattern match':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':lw},
    'Pattern match banded':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':lw},
    }

def set_all_estimator_types():
    global estimator_types
    azim = space.hrtfset.coordinates['azim']
    itd = space.itd_max*sin(azim*pi/180)
    estimator_types = (
#        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed peak'),
#        (MakeEstimator(TwoChannel, PolyClosest(6, grid=itd), difference=True,
#                                    itdmax_extend=itdmax_extend), 'Two channel'),
        (MakeEstimator(TwoChannelCrossFrequencyAlt,
                       PolyClosest(6, grid=itd),
                       itdmax_extend=itdmax_extend),
            'Two channel cross-frequency alt'),
#        #(MakeEstimator(TwoChannelBanded), 'Two channel banded'),
#        (MakeEstimator(PatternMatch), 'Pattern match'),
        (MakeEstimator(PatternMatch, normalise_banded=40), 'Pattern match banded'),
        )

use_ideal_responses = False # Removes all response noise from the results
num_shuffles = 25
training_size = 400
testing_size = 800
#num_shuffles = 5
#training_size = 400
#testing_size = 200

theaxes = {
    'wagner_owl': subplot(2, 2, 1),
    'ircam_human_uniform_bipd': subplot(2, 2, 2),
    'ircam_human_mcalpine_bd': subplot(2, 2, 3),
    'ircam_human_tollin_bd': subplot(2, 2, 4),
    }

training_filters = (
    'type=="whitenoise"',
    'i<training_size',
    )
testing_filters = (
    'type=="whitenoise"',
    'i<testing_size',
    )

@figcache
def get_plotcommands(modelname):
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
        #exec_plotcommands=False,
        formatting=formatting,
        Nticks_skip=2,
        )
    return plotcommands

for modelname in show_models:
    sca(theaxes[modelname])
    plotcommands = get_plotcommands(modelname)
    for f, args, kwds in plotcommands:
        f(*args, **kwds)
    ylim(ymin=-5)
    if 'human' in modelname:
        ylim(ymax=40)

for n, figlet in [(1, 'A'),
                  (2, 'B'),
                  (3, 'C'),
                  (4, 'D')]:
    subplot(2, 2, n)
    text(-0.3, 1.1, figlet, size=16, transform=gca().transAxes)
    
show()
