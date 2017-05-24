from base import *

__all__ = ['BDDistribution']

class BDDistribution(object):
    '''
    ``symmetry``
        Can be 'none' (BITDs will be biased towards one side), 'swap' (each
        BITD will be on the left or the right randomly) or 'duplicate' (each
        BITD is copied to both the left and the right).
    '''
    def __init__(self, f, bitd, frange=None, symmetry='swap'):
        if frange is not None:
            flow, fhigh = frange
        else:
            flow, fhigh = 0, Inf
        self.flow, self.fhigh = flow, fhigh
        I, = logical_and(f>flow, f<fhigh).nonzero()
        if symmetry=='swap':
            side = randint(2, size=len(I))*2-1
        elif symmetry=='duplicate':
            side = hstack((ones(len(I)), -ones(len(I))))
            I = hstack((I, I))
        elif symmetry=='none':
            side = ones(len(I))
        self.f = f[I]
        self.bitd = side*bitd[I]
        # filter NAN values
        I = logical_and(-isnan(self.bitd), -isnan(self.f))
        self.f = self.f[I]
        self.bitd = self.bitd[I]
        # Compute KDE estimates
        self.kde_f = gaussian_kde(self.f)
        self.kde = gaussian_kde(vstack((self.f, self.bitd)))
        # precompute for sample_bitd
        blow, bhigh = amin(self.bitd), amax(self.bitd)
        bmean = (blow+bhigh)/2
        bwidth = (bhigh-blow)/2
        blow, bhigh = bmean-2*bwidth, bmean+2*bwidth # double the range
        brange = linspace(blow, bhigh, 100)
        self.brange = brange
    
    def sample_f(self, n, frange=None):
        if frange is None:
            frange = self.flow, self.fhigh
        flow, fhigh = frange
        f = array(())
        ns = n
        while ns>0:
            newf = self.kde_f.resample(ns)
            I = logical_and(newf>flow, newf<fhigh)
            newf = newf[I]
            f = hstack((f, newf))
            ns = n-len(f)
        return f
    
    def sample_f_bitd(self, n, frange=None):
        if frange is None:
            frange = self.flow, self.fhigh
        flow, fhigh = frange
        f, bitd = array(()), array(())
        ns = n
        while ns>0:
            newf, newbitd = self.kde.resample(ns)
            I = logical_and(newf>flow, newf<fhigh)
            newf = newf[I]
            newbitd = newbitd[I]
            f = hstack((f, newf))
            bitd = hstack((bitd, newbitd))
            ns = n-len(f)
        return f, bitd
    
    def sample_bitd(self, f):
        allf = repeat(f, len(self.brange))
        allb = tile(self.brange, len(f))
        data = vstack((allf, allb))
        probs = self.kde.evaluate(data)
        probs.shape = (len(f), len(self.brange))
        cumprobs = cumsum(probs, axis=1)
        cumprobs /= reshape(amax(cumprobs, axis=1), (len(f), 1))
        p = rand(len(f))
        allb = []
        for i, (freq, prob) in enumerate(zip(f, p)):
            cp = hstack((0, cumprobs[i, :], 1))
            j = amax((prob>cp).nonzero()[0]) # max index where p>cdf
            plow, phigh = cp[j], cp[j+1]
            q = (prob-plow)/(phigh-plow) # should be in [0,1]
            blow, bhigh = self.brange[[j, j+1]]
            b = blow+q*(bhigh-blow)
            allb.append(b)
        allb = array(allb)
        return allb
