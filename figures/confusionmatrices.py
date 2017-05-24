from base import *

def compute_confusion_matrix(analysis, estimator):
    num_shuffles = analysis.settings['num_shuffles']
    bins = analysis.moresettings['itd_bins']
    binmids = 0.5*(bins[1:]+bins[:-1])
    confmat = zeros((len(bins)-1, len(bins)-1))
    shuffled_results = analysis(analysis.shuffled_results, estimator, num_shuffles)
    for shufflenum in xrange(num_shuffles):
        results, trueval, guessedval, testing_responses = shuffled_results[shufflenum]
        trueval = digitize(trueval, bins)-1
        guessedval = digitize(guessedval, bins)-1
        trueval[trueval==len(bins)-1] = len(bins)-2
        guessedval[guessedval==len(bins)-1] = len(bins)-2
        for correct, guessed in zip(trueval, guessedval):
            confmat[guessed, correct] += 1
    sumconfmat = reshape(sum(confmat, axis=0), (1, len(binmids)))
    sumconfmat[sumconfmat==0] = 1
    confmat /= sumconfmat
    return confmat

def show_confusion_matrices(analysis, estimator_types):
    basename = analysis.settings['basename']
    bins = analysis.moresettings['itd_bins']
    binmids = 0.5*(bins[1:]+bins[:-1])
    figure()
    suptitle(basename)
    for estnum, (f, name) in enumerate(estimator_types):
        estimator = f(analysis.fm_orig)
        subplot(*subplot_size[len(estimator_types)]+(estnum+1,))
        confmat = compute_confusion_matrix(analysis, estimator)
        imshow(confmat, origin='lower left', interpolation='nearest',
               aspect='auto', extent=(-amin(binmids)*second/usecond,
                                      -amax(binmids)*second/usecond,
                                      -amin(binmids)*second/usecond,
                                      -amax(binmids)*second/usecond))
        ylabel('Estimated location ($\mu$s)')
        xlabel('Location ($\mu$s)')
        title(name)
