import os
from base import *
from analyse_model import *
from cells import *
from binauraldisplay import *
from colours import *

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.13
rcParams['figure.subplot.right'] = .97
rcParams['figure.subplot.bottom'] = .107
rcParams['figure.subplot.top'] = 1.0
rcParams['figure.subplot.wspace'] = 0.25
rcParams['figure.subplot.hspace'] = 0.29
# figsize is width x height in inches
figure(figsize=(5, 7))

show_bd = True
show_results = True
std_factors = [1.5, 1, 0.75, 0.5, 0.25, 0.15]#, 0.05]
show_std_factors = [0.25, 1.5]
#limit_frequency = (0*Hz, 900*Hz)
#limit_frequency = (0*Hz, 600*Hz)
limit_frequency = (0*Hz, 1200*Hz)

#figure(figsize=(2*len(std_factors), 4), dpi=80)

max_error = 90

def set_all_estimator_types():
    global estimator_types
    estimator_types = (
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

formatting = {
    'Two channel':{'color':estimator_colours['twochannel'], 'ls':'--', 'lw':2},
    'Two channel cross-frequency':{'color':estimator_colours['twochannel'], 'ls':'-', 'lw':2},
#    'Two channel banded':{'color':'m'},
    'Pattern match':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':2},
    'Pattern match banded':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':2},
    }

training_filters = (
    'type=="whitenoise"',
    '"acoustic_noise_level" not in locals()',
    'i<training_size',
    )
testing_filters = (
    'type=="whitenoise"',
    '"acoustic_noise_level" in locals()',
    '-5<float(level-acoustic_noise_level)<5',
    'i<testing_size',
    )

use_ideal_responses = False # Removes all response noise from the results
num_shuffles = 25
training_size = 400
testing_size = 800

from models.mcalpine_guinea_pig import *

#acousticnoisemodel = NoAcousticNoise()
acousticnoisemodel = IndependentWhiteAcousticNoise(tuple(level-sn for sn in sn_level_range[::-1]))
extraname['acousticnoisemodel'] = acousticnoisemodel

if show_bd:
    for i, std_factor in enumerate(show_std_factors):
        seed(32409822)
        curbd = generate_random_mcalpine_et_al_2001_bds(cf, std_factor=std_factor)
        seed()
        #subplot(2, len(std_factors), i+1)
        subplot(4, 2, i+3)
        best_delay_distribution(cf, curbd, markersize=2)
        title('Spread factor '+str(std_factor))
        ylim(-1.5, 1.5)
        if i>0:
            ylabel('')

if show_results:
    best_error = defaultdict(list)
    best_std = defaultdict(list)
    for i, std_factor in enumerate(std_factors):
        if std_factor!=1:
            extraname['std_factor'] = std_factor
            seed(32409822)
            bd = generate_random_mcalpine_et_al_2001_bds(cf, std_factor=std_factor)
            seed()
            dl, dr = delays_from_bd(bd)
            binauraldistribution = BinauralDistribution(zip(cf, dl, dr), cfrepeat=cfrepeat)
        set_all_estimator_types()
        analysis = get_analysis_from_namespace()
        
        for estimator_type, name in estimator_types:
            Nbinuserange, error, error_std = analyse_performance_with_cells(analysis, estimator_type)
            best_error[name].append(error[-1])
            best_std[name].append(error_std[-1])
    subplot(2,1,2)
    for estimator_type, name in estimator_types:
        errorbar(std_factors,
                 array(best_error[name])/usecond,
                 array(best_std[name])/usecond,
                 label=name,
                 **formatting[name])
        i = (array(std_factors)==1).nonzero()[0][0]
        axhline(float(best_error[name][i]/usecond), ls=':', color=formatting[name]['color'])
    axvline(1, ls=':', color='k')
    xlabel('Spread factor')
    ylabel(r'Mean error ($\mu$s)')
    ylim(ymin=0)


#        subplot(2, len(std_factors), i+1+len(std_factors))
#        fit_performance_with_cells(analysis, estimator_types,
#            number_of_mso_neurons=number_of_mso_neurons[animal_name],
#            do_figure=False, do_suptitle=False, do_legend=False,
#            show_guesslevel=False,
#            show_infinity=False,
#            extra_label_spacing=3,
#        #        Nticks=Nticks,
#            #itd_unit=ms,
#            )
#        ylim(ymax=max_error)

show()
