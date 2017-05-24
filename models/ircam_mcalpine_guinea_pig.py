from base import *
from digitised_data.mcalpine_et_al_2001 import *

if int(samplerate)!=192000:
    raise ValueError('Samplerate for IRCAM guinea pig should be 192 kHz')

basename = 'ircam_mcalpine_guinea_pig'
animal_name = 'Guinea pig'
extraname = {}
space = IRCAMHRTFSpatial(subject=3012, azim_offset=-25,
                         adaptive_windowing=(512, 16, 256)
                         )
extraname['space'] = space
Nbin = 480
extraname['Nbin'] = Nbin
cfrepeat = 1 # should divide into Nbin
extraname['cfrepeat'] = cfrepeat
seed(32409822) # ensure that bd and cf are always the same
cf = erbspace(100*Hz, 1.5*kHz, Nbin/cfrepeat)
cf = repeat(cf, cfrepeat)
#bd = (pi/4*rand(Nbin)+pi/8)/(2*pi*cf)*(2*randint(2, size=Nbin)-1)
bd = generate_random_mcalpine_et_al_2001_bds(cf)
extraname['measured_bd'] = True
seed() # proper random numbers from here on

# These ERB parameters are from Shera et al. 2001 (citing Liberman 1990, Tsuji
# and Liberman 1997). Note we hstack((cf, cf)) because it is a binaural model
alpha = 0.35
beta = 4.0
gammatone_params = {'ear_Q':beta*(hstack((cf, cf))/kHz)**alpha,
                    'min_bw':0.0}
extraname['SheraERB'] = (alpha, beta)

# Compression and normalisation parameters set to give a firing rate between
# 0 and 200 Hz for the binaural neurons
compression = None
binauralmodel = PowerBinauralModel(1, 0, 0, power=8, rectify=False)
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
#    sound = whitenoise(.1*second).atlevel(70*dB)[:.3*second]
#    acousticnoisemodel = NoAcousticNoise()
#    fm = FixedModel(space, acousticnoisemodel, binauralmodel,
#                    responsenoisemodel, binauraldistribution,
#                    generate_extended_datamanager_name(basename, extraname),
#                    gammatone_params=gammatone_params,
#                    compression=compression,
#                    rms_normalisation=rms_normalisation,
#                    )
#    R, _, _, _ = fm.apply(sound, space.itd_max/2)
#    plot(bd, R, '.')
#    show()
#    