from base import *
from location import *
from analyse_model import *
from colours import *
from matplotlib import cm

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.15
rcParams['figure.subplot.right'] = .97
rcParams['figure.subplot.bottom'] = .107
rcParams['figure.subplot.top'] = 0.75
rcParams['figure.subplot.wspace'] = 0.25
rcParams['figure.subplot.hspace'] = 0.29
# figsize is width x height in inches
figure(figsize=(4, 7))

from models.mcalpine_guinea_pig import *

use_ideal_responses = False # Removes all response noise from the results
num_shuffles = 25
training_size = 400
testing_size = 2000
#num_shuffles = 5
#training_size = 400
#testing_size = 800

def lesionfunc(bd):
    return bd>0

def doplot(lesionfunc, estimator_name):
    all_estimator_types = {
        'Smoothed peak': MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max),
                                       lesionfunc=lesionfunc),
        'Pattern match': MakeEstimator(PatternMatch, lesionfunc=lesionfunc),
        'Pattern match banded': MakeEstimator(PatternMatch, normalise_banded=40, lesionfunc=lesionfunc),
        }
    
    if lesionfunc is None:
        ls = '--'
    else:
        ls = '-'
    
    formatting = {
        'Smoothed peak':{'color':estimator_colours['peak'], 'ls':ls, 'lw':2},
        'Pattern match':{'color':estimator_colours['patternmatch'], 'ls':ls, 'lw':2},
        'Pattern match banded':{'color':estimator_colours['patternmatch'], 'ls':ls, 'lw':2},
        }
    
    acousticnoisemodel = NoAcousticNoise()
    extraname['acousticnoisemodel'] = acousticnoisemodel
    
    training_filters = ('type=="whitenoise"', 'i<training_size')
    testing_filters = ('type=="whitenoise"', 'i<testing_size')
    analysis = get_analysis_from_namespace()
    
    estimator_types = tuple((all_estimator_types[n], n) for n in [
                            estimator_name,
                            ])
    
    show_performance_with_location(analysis, estimator_types, dofigure=False,
                                   which='error.smooth',
                                   formatting=formatting)

subplot(211)
#doplot(None, 'Pattern match')
#doplot(lesionfunc, 'Pattern match')
doplot(None, 'Pattern match banded')
doplot(lesionfunc, 'Pattern match banded')
xlabel('')
ylim(ymin=0, ymax=25)
subplot(212)
doplot(None, 'Smoothed peak')
doplot(lesionfunc, 'Smoothed peak')

show()
