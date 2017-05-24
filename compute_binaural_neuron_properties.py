# Make sure to reduce the size of Nbin etc. in the model you load
from models.mcalpine_guinea_pig import *

duration = 100*ms
level = 70*dB

responsenoisemodel = NoResponseNoise()

acousticnoisemodel = NoAcousticNoise()

fm = FixedModel(space, acousticnoisemodel, binauralmodel,
                responsenoisemodel, binauraldistribution,
                generate_extended_datamanager_name('binaural_'+basename,
                                                   extraname),
                gammatone_params=gammatone_params,
                compression=compression,
                rms_normalisation=rms_normalisation,
                )

X = []

sound = whitenoise(duration).atlevel(level)
itds = linspace(-space.itd_max*8, space.itd_max*8, 50)
for i, itd in enumerate(itds):
    print i,
    R, _, _, _ = fm.apply(sound, itd*second, store=False)
    print 'done'
    R /= duration
    X.append(R)
X = array(X) # shape (numlocs, numcells)
d = space.itd_max/usecond
fill([-d, d, d, -d], [0, 0, amax(X), amax(X)], color=(0.8, 0.8, 0.8))
for i, x in enumerate(X.T):
    c = float(i)/(len(X.T)-1)
    plot(itds/usecond, x, color=(c, 1-c, 0))
    plot(itds[argmax(x)]*second/usecond, x[argmax(x)], 'o', color=(c, 1-c, 0))
    plot(itds/usecond, amax(x)*(0.5+0.5*cos(2*pi*cf[i]*(itds-bd[i])))**4, '--', color=(c, 1-c, 0))
show()
