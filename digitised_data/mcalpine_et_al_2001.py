from brian import *
from brian.hears import *
from scipy.stats import truncnorm

__all__ = ['generate_random_mcalpine_et_al_2001_bds']

def fixup(s):
    s = s.replace(',', '.')
    s = s.split('\n')
    s = [map(float, w.split('    ')) for w in s if w]
    f, bitd = zip(*s)
    f = array(f)*kHz
    bitd = array(bitd)*usecond
    return f, bitd

# IGNORE THE F VALUES HERE, THEY ARE WRONG
fig_2a_means = '''
0,09455    707,10712
0,16542    520,84442
0,23318    361,37778
0,29635    277,76535
0,35333    232,09654
0,41458    182,66420
0,46000    163,59335
0,51884    205,06943
0,57556    148,14299
0,61844    113,97392
0,68096    147,91190
0,75553    117,48437
0,80553    121,18188
0,99987    109,52809
'''

fig_2a_means_plus_stds = '''
0,09879    1125,42432
0,19757    819,93372
0,30073    604,84766
0,39557    412,23495
0,49462    412,60233
0,59540    333,41052
0,68949    242,79839
0,78939    307,37531
0,89622    250,80063
0,97863    201,73302
1,09955    209,49567
1,23526    228,61478
1,34885    179,54718
1,75320    191,33490
'''

_, mean_bitd = fixup(fig_2a_means)
f, bitd_mean_plus_std = fixup(fig_2a_means_plus_stds)
std_bitd = bitd_mean_plus_std-mean_bitd

def generate_random_mcalpine_et_al_2001_bds(cf, std_factor=1.0):
    fmid = 0.5*(f[1:]+f[:-1])
    I = digitize(cf, fmid)
    mu = mean_bitd[I]
    sigma = std_bitd[I]*std_factor
    side = 2*(rand(len(mu))<.5)-1
    return (randn(len(mu))*sigma+mu)*side

if __name__=='__main__':
    if 1:
        cf = erbspace(50*Hz, 1.5*kHz, 480)
        bd = generate_random_mcalpine_et_al_2001_bds(cf)
        bd2 = (pi/4*rand(480)+pi/8)/(2*pi*cf)*(2*randint(2, size=480)-1)
        bd3 = generate_random_mcalpine_et_al_2001_bds(cf, 0.25)
        fill([amin(cf), amax(cf), amax(cf), amin(cf)],
             [-350, -350, 350, 350], color=(0.8, 0.8, 0.8))
        plot(cf, bd/usecond, '.')
        plot(cf, bd2/usecond, '.')
        plot(cf, bd3/usecond, '.')
        plot(cf, (1/(2*cf))/usecond, '-k')
        plot(cf, (1/(8*cf))/usecond, '-k')
        plot(cf, -(1/(2*cf))/usecond, '-k')
        plot(cf, -(1/(8*cf))/usecond, '-k')
        ylim(amin(bd/usecond), amax(bd/usecond))
        show()
    if 0:
        subplot(221)
        semilogx(f/kHz, mean_bitd/usecond, '.')
        subplot(223)
        semilogx(f/kHz, std_bitd/usecond, '.')
        subplot(222)
        semilogx(f/kHz, mean_bitd*f, '.')
        ylim(0, 1)
        subplot(224)
        semilogx(f/kHz, std_bitd*f, '.')
        ylim(-0.2, 0.4)
        show()
        