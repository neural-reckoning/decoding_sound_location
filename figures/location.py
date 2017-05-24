from base import *
from confusionmatrices import *

def analyse_performance_with_location(analysis, estimator):
    bins = analysis.moresettings['itd_bins']
    binmids = 0.5*(bins[1:]+bins[:-1])
    dxij = reshape(binmids, (len(binmids), 1))-reshape(binmids, (1, len(binmids)))
    absdxij = abs(dxij)
    I = compute_confusion_matrix(analysis, estimator)
    bx = sum(dxij*I, axis=0)
    mx = sum(absdxij*I, axis=0)
    return mx, bx

def show_performance_with_location(analysis, estimator_types,
                                   dofigure=True,
                                   which=None,
                                   remove_ends=2,
                                   formatting=None,
                                   ):
    if formatting is None: formatting = dict((name, {}) for _, name in estimator_types)
    basename = analysis.settings['basename']
    bins = analysis.moresettings['itd_bins']
    binmids = 0.5*(bins[1:]+bins[:-1])
    if remove_ends:
        binmids = binmids[remove_ends:-remove_ends]
    if dofigure:
        figure()
        suptitle(basename)
    for f, name in estimator_types:
        estimator = f(analysis.fm_orig)
        mx, bx = analyse_performance_with_location(analysis, estimator)
        if remove_ends:
            mx = mx[remove_ends:-remove_ends]
            bx = bx[remove_ends:-remove_ends]
        if which is None: subplot(221)
        if which is None or which=='bias':
            axhline(0, ls='--', c='k')
            plot(binmids/usecond, bx/usecond, label=name, **formatting[name])
        if which is None: subplot(223)
        if which is None or which=='bias.smooth':
            axhline(0, ls='--', c='k')
            plot(binmids/usecond, smooth(bx/usecond, 5), label=name, **formatting[name])
        if which is None: subplot(222)
        if which is None or which=='error':
            plot(binmids/usecond, mx/usecond, label=name, **formatting[name])
        if which is None: subplot(224)
        if which is None or which=='error.smooth':
            plot(binmids/usecond, smooth(mx/usecond, 5), label=name, **formatting[name])
    if which is None: subplot(221)
    if which is None: legend(loc='upper right', ncol=2)    
    if which is None or which=='bias':
        ylabel(r'Bias ($\mu$s)')
    if which is None: subplot(223)
    if which is None or which=='bias.smooth':
        ylabel(r'Smoothed bias ($\mu$s)')
        xlabel(r'ITD ($\mu$s)')
    if which is None: subplot(222)
    if which is None or which=='error':
        ylabel(r'Mean error ($\mu$s)')
    if which is None: subplot(224)
    if which is None or which=='error.smooth':
        ylabel(r'Smoothed mean error ($\mu$s)')
    xlabel(r'ITD ($\mu$s)')
