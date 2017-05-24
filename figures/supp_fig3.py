from base import *
from frequency import *
from cells import *
from cutoff_frequency import *
from stimulialpha import *
from tonefreq import *
from binauraldisplay import *
from analyse_model import *
from colours import *
from scikits.learn import *

rcParams['font.size'] = 10
rcParams['figure.subplot.left']  = 0.08
rcParams['figure.subplot.right'] = .98
rcParams['figure.subplot.bottom'] = .1
rcParams['figure.subplot.top'] = .9
rcParams['figure.subplot.wspace'] = 0.4
rcParams['figure.subplot.hspace'] = 0.3
# figsize is width x height in inches
figure(figsize=(10, 4), dpi=80)

show_guinea_pig = True
show_cat = True
show_all_frequencies = True
show_low_frequencies = True
show_broadband = True
show_bandpass = True
show_ratio = False
show_alpha = True
#bandpass_type = 'Tones'
bandpass_type = 'Band pass'
use_bandpass_noise_ratio = False
use_denominator = True
partial_denominator = False
num_frequency_bands = 10
use_cross_frequency_two_channel = False
num_ratio_plot_points = 200

lw = 2
formatting = {
    'Smoothed peak':{'color':estimator_colours['peak'], 'ls':'-', 'lw':lw},
    'Two channel':{'color':estimator_colours['twochannel'], 'ls':'-', 'lw':lw},
    'Two channel cross-frequency alt':{'color':estimator_colours['twochannel'], 'ls':':', 'lw':lw},
    'Pattern match':{'color':estimator_colours['patternmatch'], 'ls':'-', 'lw':lw},
    'Pattern match banded':{'color':estimator_colours['patternmatch'], 'ls':':', 'lw':lw},
    'Linear regression':{'color':'m', 'ls':'--', 'lw':lw},
    'Ridge regression':{'color':'m', 'ls':'-', 'lw':lw},
    'NN':{'color':'k', 'ls':'-', 'lw':lw},
    'Vector':{'color':'y', 'ls':'-', 'lw':lw},
    'MLE':{'color':'c', 'ls':'-', 'lw':lw},
    }

def set_all_estimator_types():
    global estimator_types
    estimator_types = (
        (MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed peak'),
        (MakeEstimator(TwoChannel, PolyClosest(6), difference=True,
                                    itdmax_extend=itdmax_extend), 'Two channel'),
        (MakeEstimator(Vector, PolyClosest(6),
                                    itdmax_extend=itdmax_extend), 'Vector'),
        (MakeEstimator(FitDistributionMLE), 'MLE'),
        (MakeEstimator(PatternMatch), 'Pattern match'),
        #(MakeEstimator(ScikitsLearnEstimator, linear_model.LinearRegression()), 'Linear regression'),
        (MakeEstimator(ScikitsLearnEstimator, linear_model.Ridge()), 'Ridge regression'),
        (MakeEstimator(ScikitsLearnEstimator, neighbors.NeighborsRegressor(5)), 'NN'),
        )

use_ideal_responses = False # Removes all response noise from the results
#num_shuffles = 5
#training_size = 200
#testing_size = 400
num_shuffles = 25
training_size = 400
testing_size = 800

the_filters = {
    'whitenoise': (('type=="whitenoise"', 'i<training_size'),
                   ('type=="whitenoise"', 'i<testing_size')),
    'coloured':   (('type=="whitenoise"', 'i<training_size'),
                   ('type=="powerlawnoise"', 'i<testing_size')),
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

the_models = []
if show_guinea_pig:
    the_models.append(('mcalpine_guinea_pig', 0))
if show_cat:
    the_models.append(('joris_cat', 12))
    
the_frequencies = []
if show_all_frequencies:
    the_frequencies.append((False, 0, ''))
if show_low_frequencies:
    #the_frequencies.append(((0*Hz, 900*Hz), 6, ' (low freq)'))
    the_frequencies.append(((0*Hz, 1200*Hz), 6, ' (low freq)'))
    #the_frequencies.append(((0*Hz, 600*Hz), 6, ' (low freq)'))

@figcache
def get_ratio_plot_data(use_cross_frequency_two_channel,
                        use_bandpass_noise_ratio,
                        use_denominator,
                        partial_denominator,
                        num_frequency_bands,
                        ):
    if use_cross_frequency_two_channel:
        estimator_type = MakeEstimator(TwoChannelCrossFrequency,
                                       PolyClosest(6),
                                       itdmax_extend=itdmax_extend)
    else:
        estimator_type = MakeEstimator(TwoChannel, PolyClosest(6),
                                       difference=True,
                                       itdmax_extend=itdmax_extend)
    fm = analysis.fm_orig
    cf = fm.centre_frequencies[:len(fm.centre_frequencies)/2]
    I = argsort(cf)
    estimator = estimator_type(fm)
    training_responses, testing_responses = analysis.traintest_responses(None)
    estimator.train(responses=testing_responses)
    coeffs = estimator.coefficients
    locations = []
    ratios = []
    freqs = []
    for response, location, meta in estimator.data:
        response = meta['noiseless_response']
        if use_bandpass_noise_ratio:
            f = meta[bandpass_key]
            ratio = estimator.func(response)
            locations.append(location)
            ratios.append(ratio)
            freqs.append(f)
        else:
            response = response[I]
            coefficients = estimator.coefficients[I]
            rc = response*coefficients
            if not partial_denominator:
                denominator = sum(response)
            J = hstack((arange(0, len(cf), len(cf)/num_frequency_bands), len(cf)))
            for jstart, jend in zip(J[:-1], J[1:]):
                f = mean(cf[I][jstart:jend])*Hz
                if partial_denominator:
                    denominator = sum(response[jstart:jend])
                if not use_denominator:
                    denominator = 1.0
                partialratio = sum(rc[jstart:jend])/denominator
                locations.append(location)
                ratios.append(partialratio)
                freqs.append(f)
    locations = array(locations)
    ratios = array(ratios)
    freqs = array(freqs)
    return locations, ratios, freqs
        
for modelname, offset in the_models:
    exec 'from models.%s import *' % modelname
    set_all_estimator_types()
    acousticnoisemodel = NoAcousticNoise()
    extraname['acousticnoisemodel'] = acousticnoisemodel
    
    subplot(2, 4, 3+4*(modelname=='joris_cat'))
    limit_frequency = False
    training_filters, testing_filters = the_filters['whitenoise']
    analysis = get_analysis_from_namespace()
    fig_performance_with_cutoff_frequency(analysis, estimator_types,
              formatting=formatting)
    if modelname=='joris_cat':
        ylim(ymin=0, ymax=100)
    else:
        ylim(ymin=0, ymax=50)
        xlabel('')
    xlim(xmin=0.1, xmax=1.5)

    for limit_frequency, offset2, extralabel in the_frequencies:

        if show_broadband:
            training_filters, testing_filters = the_filters['whitenoise']
            analysis = get_analysis_from_namespace()
            Ncells = analysis.fm.cfN
            if Ncells==480:
                Nticks = [40, 160, 320, 480]
            else:
                Nticks = [40, 200, Ncells]
            # error with num cells
            subplotnum = 1
            if limit_frequency is not False:
                subplotnum += 2
            if modelname=='joris_cat':
                subplotnum += 4
                ymax = 70
            else:
                ymax = 30
            subplot(2, 4, subplotnum)
            #subplot(4, 6, 1+offset+offset2)
            if limit_frequency is False:
                fit_performance_with_cells(analysis, estimator_types,
                    number_of_mso_neurons=number_of_mso_neurons[animal_name],
                    do_figure=False, do_suptitle=False, do_legend=False,
                    show_guesslevel=False,
                    show_infinity=False,
                    show_mso=False,
                    extra_label_spacing=3,
                    Nticks=Nticks,
                    formatting=formatting,
                    show_fit=False,
                    do_xlabel=(modelname=='joris_cat'),
                    #itd_unit=ms,
                    )
                ylim(ymax=ymax)
                ylabel(animal_name+extralabel+'\n'+gca().get_ylabel())
#            if offset+offset2==0:
#                title('Broadband')
                    
        if show_ratio:
            training_filters, testing_filters = the_filters['whitenoise']
            
        if show_alpha and limit_frequency is not False and modelname!='joris_cat':
            training_filters, testing_filters = the_filters['coloured']
            analysis = get_analysis_from_namespace()
            #subplot(4, 6, 5+offset+offset2)
            subplot(2, 4, 4)
            fig_performance_with_stimuli_alpha(analysis, estimator_types,
                do_figure=False, do_suptitle=False, do_legend=False,
                show_mean_error=True, show_bias=False,
                show_xlabel=False, show_ylabel=True,
                formatting=formatting,
                )
#            if offset+offset2==0:
#                title('Coloured\nError')
            #subplot(4, 6, 6+offset+offset2)
            subplot(2, 4, 8)
            fig_performance_with_stimuli_alpha(analysis, estimator_types,
                do_figure=False, do_suptitle=False, do_legend=False,
                show_mean_error=False, show_bias=True,
                show_xlabel=True, show_ylabel=True,
                formatting=formatting,
                )
#            if offset+offset2==0:
#                title('Coloured\nBias')
                
    if show_ratio:
        if use_bandpass_noise_ratio:
            acousticnoisemodel = IndependentWhiteAcousticNoise((level, level))
            extraname['acousticnoisemodel'] = acousticnoisemodel
            training_filters, testing_filters = the_filters[bandpass_type]
        else:
            training_filters, testing_filters = the_filters['whitenoise']            
        for limit_frequency, offset2, extralabel in the_frequencies:
            _testing_size, testing_size = testing_size, num_ratio_plot_points
            analysis = get_analysis_from_namespace()
            testing_size = _testing_size
            locations, ratios, freqs = get_ratio_plot_data(use_cross_frequency_two_channel,
                                                           use_bandpass_noise_ratio,
                                                           use_denominator,
                                                           partial_denominator,
                                                           num_frequency_bands)
            subplot(4, 6, 4+offset+offset2)
            nfreqs = (freqs-amin(freqs))/(amax(freqs)-amin(freqs))
            cols = [(c, 1-c, 0, 1) for c in nfreqs]
            scatter(locations/msecond, ratios, facecolors=cols, edgecolors='none',
                    s=3)
            xlabel('ITD (ms)')
            ylabel('Ratio')
            axis('tight')
            if offset+offset2==0:
                if use_cross_frequency_two_channel:
                    title('Cross frequency\nTwo channel ratio')
                else:
                    title('Two channel ratio')

    if show_bandpass:
        acousticnoisemodel = IndependentWhiteAcousticNoise((level, level))
        extraname['acousticnoisemodel'] = acousticnoisemodel
        for limit_frequency, offset2, extralabel in the_frequencies:
            if limit_frequency is not False:
                continue
            training_filters, testing_filters = the_filters[bandpass_type]
            analysis = get_analysis_from_namespace()
            # error
            subplot(2, 4, 2+4*(modelname=='joris_cat'))
            #subplot(4, 6, 2+offset+offset2)
            fig_performance_with_stimuli_tonefreq(analysis, estimator_types,
                    tonefreqkey=bandpass_key,
                    do_figure=False, do_suptitle=False,
                    show_mean_error=True,
                    show_bias=False,
                    show_legend=False,
                    frequency_unit=kHz,
                    formatting=formatting,
                    show_xlabel=(modelname=='joris_cat'),
                    )
            xlim(0.15, 1.45)
            #axis('tight')
            if limit_frequency:
                #xlim(*array(limit_frequency)/kHz)
                xlim(xmax=limit_frequency[1]/kHz)
            ylim(0)
#            if offset+offset2==0:
#                title('%s\nError'%bandpass_type)
#            # bias
#            subplot(4, 6, 3+offset+offset2)
#            fig_performance_with_stimuli_tonefreq(analysis, estimator_types,
#                    tonefreqkey=bandpass_key,
#                    do_figure=False, do_suptitle=False,
#                    show_mean_error=False,
#                    show_bias=True,
#                    show_legend=False,
#                    frequency_unit=kHz,
#                    formatting=formatting,
#                    )
#            #axis('tight')
#            if limit_frequency:
#                #xlim(*array(limit_frequency)/kHz)
#                xlim(xmax=limit_frequency[1]/kHz)
#            if offset+offset2==0:
#                title('%s\nBias'%bandpass_type)
    

for ax, figlet in enumerate('ABCD'):
    subplot(2, 4, ax+1)
    text(-0.3, 1.1, figlet, size=16, transform=gca().transAxes)
      
show()
