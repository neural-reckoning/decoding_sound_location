'''
This utility should be run after generating the data from runsimschedule.sh,
it copies the data from the noise-free runs to the data directory of the noisy
runs so that noiseless data can be used for training.
'''
from common import *
import shutil

basepath = os.path.join(datapath, 'optimal_localisation_decoding')

_, dirnames, _ = os.walk(basepath).next()

files = {}
nonoise = {}

for dirname in dirnames:
    _, _, filenames = os.walk(os.path.join(basepath, dirname)).next()
    filenames.sort()
    files[dirname] = filenames
    if 'acousticnoisemodel_NoAcousticNoise' in filenames:
        nonoise[dirname] = [f for f in filenames if not f.startswith('acousticnoisemodel')]

copies_to_make = []

for dirname, basefnames in nonoise.items():
    for dirname2, fnames in files.items():
        fnames = [f for f in fnames if not f.startswith('acousticnoisemodel')]
        if len(fnames)<len(files[dirname2]): # i.e. there was an acousticnoisemodel
            if str(fnames)==str(basefnames) and dirname2!=dirname:
                copies_to_make.append((dirname, dirname2))

for frompath, topath in copies_to_make:
    frompath = os.path.join(basepath, frompath, 'data.data')
    topath = os.path.join(basepath, topath, 'data.data')
    for name in os.listdir(frompath):
        src = os.path.join(frompath, name)
        dest = os.path.join(topath, name+'.nonoise')
        try:
            os.symlink(src, dest)
            print 'Created symlink', dest, '->', src
        except:
            shutil.copy(src, dest)
            print 'Copied', src, '->', dest
