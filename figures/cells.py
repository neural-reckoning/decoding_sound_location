#import os
#os.environ['SAMPLERATE'] = '97656' # cat
from base import *
from scipy import stats
from external_data.number_of_mso_neurons import *

def analyse_performance_with_cells(analysis, estimator_type, angles=False):
    Nbinuserange = analysis.default_Nbinuserange
    meanerror = analysis(analysis.cells_results, estimator_type, angles=angles)
    error = []
    error_std = []
    for Nbinuse in Nbinuserange:
        me = meanerror[Nbinuse]
        me = [e for e in me if not isnan(e)]
        if len(me):
            mestd = std(me)
            me = mean(me)
        else:
            mestd = nan
            me = nan
        error.append(me)
        error_std.append(mestd)
    return Nbinuserange, array(error), array(error_std)

def fig_performance_with_cells(analysis, estimator_types, angles=False):
    if angles:
        errorunit = 1
    else:
        errorunit = usecond
    basename = analysis.settings['basename']
    figure()
    suptitle(basename)
    if angles:
        axhline(60, ls='--', color='k')
    else:
        axhline((2./3)*float(analysis.settings['space'].itd_max/usecond),
                ls='--', color='k')
    for estimator_type, name in estimator_types:
        Nbinuserange, error, error_std = analyse_performance_with_cells(analysis,
                                                estimator_type, angles=angles)
        errorbar(Nbinuserange, error/errorunit, error_std/errorunit, label=name)
    legend(loc='upper right', ncol=2)
    xlabel('Number of cells')
    if angles:
        ylabel('Mean error (deg)')
    else:
        ylabel('Mean error (us)')
    ylim(ymin=0)

def fit_performance_with_cells(analysis, estimator_types, number_of_mso_neurons=5000,
                               do_figure=True, do_suptitle=True, do_legend=False,
                               do_xlabel=True, do_ylabel=True,
                               show_guesslevel=True,
                               colours=None,
                               itd_unit=usecond,
                               show_mso=True, show_infinity=True,
                               show_fit=True,
                               extra_label_spacing=1,
                               Nticks=None,
                               Nticks_skip=1,
                               angles=False,
                               exec_plotcommands=True,
                               formatting=None,
                               tight=True):
    if formatting is None: formatting = dict((name, {}) for _, name in estimator_types)
    plotcommands = []
    if itd_unit is usecond:
        itd_unit_name = '$\mu$s'
    else:
        itd_unit_name = str(itd_unit)
    if angles:
        itd_unit = 1
        itd_unit_name = 'deg'
    basename = analysis.settings['basename']
    if do_figure: plotcommands.append((figure, (), {}))#figure()
    if do_suptitle: plotcommands.append((suptitle, (basename,), {}))#suptitle(basename)
    if colours is None:
        colours = dict((name, c) for (_, name), c in zip(estimator_types, 'bgrcmyk'))
    #cols = 'bgrcmy'
    if show_guesslevel:
        if angles:
            #axhline(60, ls='--', color='k')
            plotcommands.append((axhline, (60,), {'ls':'--', 'color':'k'}))
        else:
            #axhline((2./3)*float(analysis.settings['space'].itd_max/itd_unit),
            #        ls='--', color='k')
            plotcommands.append((axhline,
                                 ((2./3)*float(analysis.settings['space'].itd_max/itd_unit),),
                                 {'ls':'--', 'color':'k'}))
    for estimator_type, name in estimator_types:
        col = colours[name]
        Nbinuserange, error, error_std = analyse_performance_with_cells(analysis,
                                                estimator_type, angles=angles)
        I = -isnan(error)
        x = Nbinuserange[I]
        y = error[I]/itd_unit
        f = lambda (a, b, c), x: a+b*x**c
        w = 1.0/error_std[I]#**2
        #w = x
        g = lambda p: (f(p, x)-y)*w
        ssmin = Inf
        p = None
        for a0 in linspace(0, amin(y), 1000):
            slope, intercept, _, _, _ = stats.linregress(log(x), log(y-a0))
            b0 = exp(intercept)
            c0 = slope
            ss = sum(g((a0, b0, c0))**2)
            if ss<ssmin or p is None:
                ssmin = ss
                p = (a0, b0, c0)
        #print name, p
        #errorbar(x, y, yerr=error_std[I]/itd_unit, color=col, ls='--', marker='.')
        fmt = {'yerr':error_std[I]/itd_unit,
                             #'color':col,
                             'ls':'--',
                             #'marker':'.'
                             }
        fmt.update(formatting[name])
        #print errorbar, (x, y), fmt
        plotcommands.append((errorbar, (x, y), fmt))
        d = diff(Nbinuserange[I])[0]
        xt = Nbinuserange[I]
        Nv = Nbinuserange[I]
        if Nticks is None:
            tx = Nbinuserange[I][::Nticks_skip]
            xt0 = Nbinuserange[I][::Nticks_skip]
        else:
            tx = Nticks
            xt0 = Nticks
        extra_labels = []
        if show_mso:
            Nv = hstack((Nv, number_of_mso_neurons))
            xt = hstack((xt, xt[-1]+extra_label_spacing*d))
            tx = hstack((tx, xt[-1]))
            extra_labels.append('MSO')
        if show_infinity:
            Nv = hstack((Nv, 1000000000000))
            xt = hstack((xt, xt[-1]+extra_label_spacing*d))
            tx = hstack((tx, xt[-1]))
            extra_labels.append('$\\infty$')
        #plot(xt, f(p, Nv), label=name, color=col, ls='-', marker='.')
        #xticks(tx, map(str, xt0)+extra_labels)
        fmt = {'label':'name',# 'color':col,
                            'ls':'-',# 'marker':'.'
                            }
        fmt.update(formatting[name])
        if show_fit:
            plotcommands.append((plot, (xt, f(p, Nv)), fmt))
        plotcommands.append((xticks, (tx, map(str, xt0)+extra_labels), {}))
    if do_legend:
        plotcommands.append((legend, (), {'loc':'upper right', 'ncol':1}))
        #legend(loc='upper right', ncol=1)
    if do_xlabel:
        #xlabel('Number of cells')
        plotcommands.append((xlabel, ('Number of cells',), {}))
    if do_ylabel:
        if angles:
            #ylabel('Mean error (deg)')
            plotcommands.append((ylabel, ('Mean error (deg)',), {}))
        elif itd_unit is usecond:
            plotcommands.append((ylabel, (r'Mean error ($\mu$s)',), {}))
        else:
            #ylabel('Mean error (%s)' % itd_unit)
            plotcommands.append((ylabel, ('Mean error (%s)'%itd_unit,), {}))
    #ylim(ymin=0)
    plotcommands.append((ylim, (), {'ymin':0}))
    if exec_plotcommands:
        for f, args, kwds in plotcommands:
            f(*args, **kwds)
    return plotcommands

if __name__=='__main__':
    from analyse_model import * # load default values of some parameters
    # change this to change the model
    #from models.joris_cat import *
    #from models.joris_tollin_cat import *
    #from models.ircam_human_uniform_bitd import *
    #from models.ircam_human_uniform_bipd import *
    from models.mcalpine_guinea_pig import *
    use_ideal_responses = False # Removes all response noise from the results
    num_shuffles = 5
    training_size = 400
    testing_size = 200
    limit_frequency = (0*Hz, 900*Hz)
    acousticnoisemodel = NoAcousticNoise()
    extraname['acousticnoisemodel'] = acousticnoisemodel
    training_filters = (
        'type=="whitenoise"',
        'i<training_size',
        )
    testing_filters = (
        'type=="whitenoise"',
        'i<testing_size',
        )
    estimator_types = (
    ##### JEFFRESS LIKE #################
        #(MakeEstimator(Jeffress, SmoothedMax(0.02*space.itd_max), phasemultiply=True), 'Smoothed Jeffress (phasemult)'),
        #(MakeEstimator(Jeffress, SmoothedMax(0.15*space.itd_max)), 'Smoothed Jeffress'),
            # Use next for IRCAM IPD
        #(MakeEstimator(TrainedJeffress, DiscreteMax(window_len=13), bdmaximiser=DiscreteMax(13)), 'Jeffress (trained, smoothed, discrete)'),
        #(MakeEstimator(TrainedJeffress, SmoothedMax(0.05*space.itd_max), bdmaximiser=FitSmoothedMax(0.05*space.itd_max, neighbourhood=0.075)), 'Jeffress (trained, smoothed)'),
            # Use next TrainedJeffress for Tollin cat
        #(MakeEstimator(TrainedJeffress, SmoothedMax(0.05*space.itd_max), bdmaximiser=FitSmoothedMax(0.05*space.itd_max, neighbourhood=0.15)), 'Jeffress (trained'),
    ##### PATTERN MATCHING ###############
        (MakeEstimator(PatternMatch), 'Pattern match'),
    ##### TWO CHANNEL ####################
        (MakeEstimator(TwoChannel, PolyClosest(6), itdmax_extend=itdmax_extend), 'Two channel'), # BEST
        (MakeEstimator(TwoChannelCrossFrequencyAlt, PolyClosest(6), itdmax_extend=itdmax_extend), 'Two channel cross-frequency alt'),
    ##### SCIKITS.LEARN REGRESSION #######
      ### Linear Regression ###
        #(MakeEstimator(ScikitsLearnEstimator, LinearRegression()), 'Linear regression'), # OK
        #(MakeEstimator(ScikitsLearnEstimator, RidgeCV()), 'Ridge regression CV'), # GOOD
      ### Nearest neighbours ###
        #(MakeEstimator(ScikitsLearnEstimator, NeighborsRegressor()), 'Nearest neighb'), # Excellent
        )

    analysis = get_analysis_from_namespace()
    
#    fig_performance_with_cells(analysis, estimator_types, angles=True)
    fit_performance_with_cells(analysis, estimator_types,
                               number_of_mso_neurons=number_of_mso_neurons[animal_name],
                               angles=True,
                               )
    try:
        title('std_factor=%s'%std_factor)
    except:
        pass
    show()
    