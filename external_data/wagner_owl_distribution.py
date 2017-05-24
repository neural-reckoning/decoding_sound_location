from base import *
from bd_distribution import *

__all__ = ['get_wagner_bd_data',
           'WagnerOwlDistribution',
           ]

def get_wagner_bd_data():
    base, _ = os.path.split(__file__)
    fname = open(os.path.join(base, 'best_delays_owl.txt'))
    data = loadtxt(fname, delimiter=',')
    cf, bd = data.T
    return cf, bd*usecond


class WagnerOwlDistribution(BDDistribution):
    def __init__(self, frange=None, symmetry='swap'):
        f, bitd = get_wagner_bd_data()
        BDDistribution.__init__(self, f, bitd, frange=frange, symmetry=symmetry)


if __name__=='__main__':
    if 0:
        dist = WagnerOwlDistribution(symmetry='duplicate',
                                     frange=(100*Hz, 1.5*kHz))
        cf = erbspace(100*Hz, 1.5*kHz, 480)
        bd = dist.sample_bitd(cf)
#        plot(cf, bd, '.')

        fmin, fmax = 100*Hz, 1500*Hz
        tmin, tmax = -1.5, 1.5
        frange = linspace(fmin, fmax, 100)
        trange = linspace(tmin*ms, tmax*ms, 100)

        pf, pt = meshgrid(frange, trange)
        s = pf.shape
        pf.shape = pf.size
        pt.shape = pt.size
        pp = 2*pi*pf*pt

        f, bitd = dist.f, dist.bitd
        print len(f), len(f)/2
        g = dist.kde
        img = g.evaluate(vstack((pf, pt)))
        img.shape = s
        img = img/sum(img, axis=0)
        subplot(121)
        imshow(img, origin='lower left',
               extent=(fmin, fmax, tmin, tmax),
               aspect='auto')
        plot(f, bitd/ms, '.', color=(1,1,1))
        title('CF/BITD conditional on CF')
        subplot(122)
        plot(cf, bd, '.')
        show()
    if 1:
        cf, bd = get_wagner_bd_data()
        plot(cf, bd, '.')
        show()
