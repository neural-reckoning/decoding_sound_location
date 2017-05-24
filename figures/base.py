from brian import *
from brian.hears import *
from matplotlib import cm, rcParams
import sys, os, glob
sys.path = [os.path.normpath(os.path.join(os.path.split(__file__)[0], '../'))]+sys.path 
from common import *

def figcache(f):
    basename = os.path.split(inspect.getsourcefile(f))[1].split('.')[0]
    basename += '.'+f.__name__
    basedir = os.path.normpath(os.path.join(os.path.split(__file__)[0], '../../../data/optimal_localisation_decoding/figcache/'))
    def func(*args, **kwds):
        name = ','.join(map(str, args)+[k+'='+str(v) for k, v in kwds.iteritems()])
        fname = os.path.join(basedir, basename+'_'+name)
        if os.path.exists(fname):
            print 'Using cached value of %s(%s)'%(basename, name)
            return pickle.load(open(fname, 'rb'))
        print 'Generating data for %s(%s)'%(basename, name)
        val = f(*args, **kwds)
        pickle.dump(val, open(fname, 'wb'), -1) 
        return val
    return func

def smoothed_scatter(x, y, w):
    I = argsort(x)
    x = x[I]
    y = y[I]
    y = gaussian_smooth(x, y, w)
    return x, y

if __name__=='__main__':
    @figcache
    def test(x, y, z=5):
        return [1,2,x,y,z]
    print test(3, 4, z=6)
