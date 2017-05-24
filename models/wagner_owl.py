#import os
#os.environ['SAMPLERATE'] = '96000'
from base import *

if int(samplerate)!=96000:
    raise ValueError('Samplerate for Wagner owl data should be 96000 Hz')

basename = 'wagner_owl'
animal_name = 'Owl'
extraname = {}
space = WagnerOwlHRTFSpatial(adaptive_windowing=(512, 16, 256))
extraname['space'] = space
Nbin = 480
extraname['Nbin'] = Nbin
cfrepeat = 1 # should divide into Nbin
extraname['cfrepeat'] = cfrepeat
duplicate = True
seed(32409822) # ensure that bd and cf are always the same

# We use a lower limit of 2 kHz because owls are assumed not to use ITD cues
# below this frequency, and an upper limit of 8 kHz because we do not have 
# BD distribution above this frequency
flow, fhigh = 2*kHz, 8*kHz
if duplicate:
    owl_dist = WagnerOwlDistribution((flow, fhigh), symmetry='duplicate')
    extraname['owl_dist_symmetry_duplicate'] = True
else:
    owl_dist = WagnerOwlDistribution((flow, fhigh))
use_erbspace = True
if use_erbspace:
    cf = erbspace(flow, fhigh, Nbin/cfrepeat)
    extraname['use_erbspace'] = True
else:
    cf = owl_dist.sample_f(Nbin/cfrepeat)
cf = repeat(cf, cfrepeat)
bd = owl_dist.sample_bitd(cf)
seed() # proper random numbers from here on

# Parameters from Koppl 1997 "Frequency tuning and spontaneous activity in the
# auditory nerve and cochlear nucleus magnocellularis of the barn owl Tyto alba"
# Original formula, Q10 = 0.074*CF**0.504, but for gammatone we have
# BW10 = 1.76b = 1.76*1.019*ERB = 1.76*1.019*((f/ear_Q)**erb_order+min_bw),
# and Q10 = f/BW10
# and solving we get the parameters below
gammatone_params = {'ear_Q':0.132*hstack((cf, cf))**0.504,
                    'min_bw':0.0, 'erb_order':1}

# Compression and normalisation parameters set to give a firing rate between
# 0 and 200 Hz for the binaural neurons
compression = None
# Use model from Fischer, Christianson, Pena 2008
binauralmodel = PowerBinauralModel(1, 0, 0, power=2, rectify=False)
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
#    from figures.binauraldisplay import *
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
#    subplot(211)
#    best_delay_distribution(cf, bd)
#    subplot(212)
#    display_response(cf, bd, R, itd=None, min_marker_size=1, max_marker_size=2,
#                     mincol=0.0, maxcol=1.0)
#    show()
#    