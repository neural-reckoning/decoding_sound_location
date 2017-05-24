from base import *

basename = 'joris_cat'
animal_name = 'Cat'
extraname = {}
space = Spatial(itd_max=400*usecond)
#space = DiscreteSpatial(itd_max=400*usecond, N=19)
#space = SphericalHeadSpatial(radius=67*mmetre, distance=2*metre, azim_ear=110, elev_ear=40)
extraname['space'] = space
Nbin = 480
extraname['Nbin'] = Nbin
cfrepeat = 1 # should divide into Nbin
extraname['cfrepeat'] = cfrepeat
duplicate = True
seed(32409822) # ensure that bd and cf are always the same
if duplicate:
    cat_dist = JorisCatDistribution((100*Hz, 1.5*kHz), symmetry='duplicate')
    extraname['cat_dist_symmetry_duplicate'] = True
else:
    cat_dist = JorisCatDistribution((100*Hz, 1.5*kHz))
use_erbspace = True
if use_erbspace:
    cf = erbspace(100*Hz, 1.5*kHz, Nbin/cfrepeat)
    extraname['use_erbspace'] = True
else:
    cf = cat_dist.sample_f(Nbin/cfrepeat)
cf = repeat(cf, cfrepeat)
bd = cat_dist.sample_bitd(cf)
seed() # proper random numbers from here on

# These ERB parameters are from Shera et al. 2001 (citing Liberman 1990, Tsuji
# and Liberman 1997). Note we hstack((cf, cf)) because it is a binaural model
alpha = 0.37
beta = 5.0
gammatone_params = {'ear_Q':beta*(hstack((cf, cf))/kHz)**alpha,
                    'min_bw':0.0}
extraname['SheraERB'] = (alpha, beta)

# Compression and normalisation parameters set to give a firing rate between
# 0 and 200 Hz for the binaural neurons
compression = None
binauralmodel = PowerBinauralModel(1, 0, 0, power=4, rectify=False)
max_firing_rate = 200*Hz
rms_normalisation = 0.5*float(max_firing_rate)**(1./binauralmodel.power)
extraname['compression'] = compression
extraname['binauralmodel'] = binauralmodel
extraname['rms_normalisation'] = rms_normalisation

responsenoisemodel = PoissonResponseNoise()
extraname['responsenoisemodel'] = responsenoisemodel

print 'CF range (Hz):', amin(cf), amax(cf)
print 'ITD range (usecond):', amin(bd/usecond), amax(bd/usecond)

dl, dr = delays_from_bd(bd)
binauraldistribution = BinauralDistribution(zip(cf, dl, dr), cfrepeat=cfrepeat)

#if __name__=='__main__':
#    sound = whitenoise(.1*second).atlevel(70*dB)
#    acousticnoisemodel = NoAcousticNoise()
#    fm = FixedModel(space, acousticnoisemodel, binauralmodel,
#                    responsenoisemodel, binauraldistribution,
#                    generate_extended_datamanager_name(basename, extraname),
#                    gammatone_params=gammatone_params,
#                    compression=compression,
#                    rms_normalisation=rms_normalisation,
#                    )
#    R, _, _, _ = fm.apply(sound, space.itd_max/2)
#    plot(R, '.')
#    show()
#    