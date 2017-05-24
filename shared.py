'''
Shared imports, etc. modified from newlocalisation.
'''

from brian import *
from brian.hears import *
from brian.tools import datamanager
from brian.utils.progressreporting import *
import os, sys, time, multiprocessing, Tkinter, pickle, gc, glob, datetime
import random
import numpy, scipy
import playdoh
from scipy.io import *
from scipy.stats import gaussian_kde
from collections import defaultdict

# shared samplerate, we use this for everything resampling if necessary
if 'SAMPLERATE' in os.environ:
    samplerate = int(os.environ['SAMPLERATE'])*Hz
else:
    samplerate = 44.1*kHz
    #samplerate = 192*kHz # Victor's IRCAM data
    #samplerate = 97656.25*Hz # Tollin cat
set_default_samplerate(samplerate)

# this value is common to all simulations and says that there will never be
# a delay longer than this, which we use to normalise the application of noise
# to left/right delayed sounds so that the total amount of noise compared to
# signal is always the same. The value should be as small as possible. This
# value covers the delay for the two channel model down to 50Hz assuming best
# phases go up to 3pi/8. 
global_max_delay = 5*ms
# In addition, gammatone filtered sounds continue to ring after the end of the
# sound, in order to equalise between two delayed versions, we have to allow the
# ringing to completely decay, this time allows a 50Hz signal to almost entirely
# decay to 0.
max_decay_time = 100*ms

# base path for data, and derived DataManager class which uses it
datapath, _ = os.path.split(__file__)
datapath = os.path.normpath(os.path.join(datapath, 'data'))
class DataManager(datamanager.DataManager):
    def __init__(self, name):
        datamanager.DataManager.__init__(self, name, datapath)

def search_for_path(paths, errmsg=None):
    for path in paths:
        if os.path.exists(path):
            return path
    else:
        s = 'Cannot find any of the specified paths'
        if errmsg is not None:
            s += ' ('+errmsg+')'
        s += ', add file locations to shared.py'
        raise IOError(s)

# convenience function to get the IRCAM database, if everyone just adds the
# location of their IRCAM database to the ircam_locations list, then it will
# automatically be loaded on all of our computers without having to keep a
# configuration file.
ircam_locations = [
    r'/home/bertrand/Data/Measurements/HRTF/IRCAM',               
    r'D:\HRTF\IRCAM',
    r'C:\Documents and Settings\dan\My Documents\Programming\IRCAM',
    r'C:\HRTF\IRCAM',
    r'F:\HRTF\IRCAM',
    r'/home/dan/programming/hrtf/IRCAM',
    r'C:\Users\Dan\programming\HRTF\IRCAM',
    ]
def get_ircam():
    path = search_for_path(ircam_locations, 'IRCAM')
    ircam = IRCAM_LISTEN(path)
    return ircam

# Similar to the above for Dan Tollin's CAT HRTF database
dan_tollin_cat_hrtf_locations = [
    r'D:\HRTF\DanTollinCat',
    r'F:\HRTF\DanTollinCat',
    r'/home/dan/programming/hrtf/DanTollinCat',
    r'C:\Users\Dan\programming\HRTF\DanTollinCat',
    ]

# Wagner owl locations
wagner_owl_locations = [
    r'/home/dan/programming/hrtf/wagner_owl',
    r'C:\Users\Dan\programming\HRTF'
    ]

# From http://www.scipy.org/Cookbook/SignalSmooth
def smooth(x,window_len=11,window='hanning'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        window_len = x.size
        #raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def gaussian_smooth(x, y, w, xout=None):
    if xout is None:
        xout = x
    xout = reshape(xout, (1, len(xout)))
    x = reshape(x, (len(x), 1))
    weight = exp(-(x-xout)**2/(2*w**2))
    sumweight = sum(weight, axis=0)
    sumweight = reshape(sumweight, (1, len(sumweight)))
    weight /= sumweight
    xout = xout.flatten()
    y = reshape(y, (len(y), 1))
    return sum(y*weight, axis=0)

def generate_extended_datamanager_name(basename, extraname, manbasepath='optimal_localisation_decoding/'):
    if int(samplerate)!=44100:
        extraname['samplerate'] = samplerate
    for k in extraname.keys():
        extraname[k] = str(extraname[k])
    extraname_hash = hash(tuple(sorted(extraname.items())))
    ext_name = basename+'_'+str(extraname_hash)
    temp_dataman = DataManager(manbasepath+ext_name)
    basepath = temp_dataman.basepath
    del temp_dataman
    gc.collect()
    for k, v in extraname.iteritems():
        p = os.path.join(basepath, k+'_'+v)
        f = open(p, 'w')
    return ext_name+'.data'+'/data'

subplot_size = {1:(1, 1), 2:(1, 2), 3:(1, 3), 4:(2, 2),
                5:(2, 3), 6:(2, 3), 7:(2, 4), 8:(2, 4),
                9:(3, 3), 10:(3, 4), 11:(3, 4), 12:(3, 4),
                13:(4, 4), 14:(4, 4), 15:(4, 4), 16:(4, 4)}

if __name__=='__main__':
    x = randn(100)
    x.sort()
    x += .1*randn(100)
    y = smooth(x, window_len=50)
    plot(x, '.')
    plot(y)
    show()