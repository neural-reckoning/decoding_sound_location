from shared import *
import itertools
import mdp
import pickle
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy import optimize
from scipy.stats import poisson
from scipy.special import gammaln
from scipy.stats.kde import gaussian_kde
from copy import copy
import hashlib
import inspect
from scipy import weave

_have_sample_figs = True
def no_more_sample_figs():
    global _have_sample_figs
    _have_sample_figs = False

class MakeEstimator(object):
    def __init__(self, estclass, *args, **kwds):
        self.estclass = estclass
        self.args = args
        self.kwds = kwds
    def __call__(self, fixedmodel):
        return self.estclass(*((fixedmodel,)+self.args), **self.kwds)
    def __str__(self):
        kwds = map(str, self.kwds.items())
        kwds.sort()
        strhash = str(hash(str((self.args, kwds))))
        return self.estclass.__name__+strhash

class Estimator(object):
    '''
    Base class for machine-learning estimators.
    
    Estimator class should implement the following methods:
    
    ``__init__``
        Should set the value of the attribute ``fixedmodel`` if the
        default methods below are to be used. It can do this by calling
        ``Estimator.__init__(fixedmodel)``.
        
    ``start_training()``
        Optional: called at the beginning of the training session.
        
    ``end_training()``
        Optional: called at the end of the training session.
        
    ``train_response(response, location)``
        Incremental training, the given ``response`` corresponds to the given
        ``location``.
        
    ``test_response(response)``
        Should return the estimated location based on ``response``.
        
    ``save(filename[, datapath])``
        Save the training data. The datapath can be changed, by default it
        points to ../../data.
        
    ``load(filename[, datapath])``
        Load the training data. The datapath can be changed, by default it
        points to ../../data.
    
    The following methods are implemented by default:
    
    ``train(responses=None, start_training=True, end_training=True)``
        Train the estimator with a set of stimuli and responses. By default,
        uses the saved responses of the fixedmodel, but you can provide a 
        different iterable container of elements of the form
        ``(response, location, meta, sound)``.
        Calls the method ``train_response`` which should be
        implemented by the Estimator, and, if start_training and end_training
        are True (the default) will call the start_training and end_training
        methods at the beginning and end. If you want to train with multiple
        stimuli sets, you should call with only start_training=True,
        end_training=False for the first stimuli set, both False for all
        but the last stimuli set, and values (False, True) for the last one.
        
    ``test(responses=None)``
        Test the estimator on a set of stimuli and responses, as for ``train()``
        it uses the saved responses of the fixedmodel by default.
        Should return a list of
        elements ``(location, estimate, meta, sound)`` where ``location`` is
        the true location of the sound, ``estimate`` is the estimated location,
        ``meta`` is the stimulus meta-data, and ``sound`` is the sound file
        itself (or None if store_sounds=False). Calls the method
        ``test_response`` which should be implemented by the Estimator.    
    '''
    def __init__(self, fixedmodel,
                 maximiser=None, samplegrid=400,
                 samplefigs=0):
        self.fixedmodel = fixedmodel
        self.bd = bd = self.fixedmodel.binauraldistribution
        numcell = len(bd.best_frequencies)
        self.numcell = numcell
        self.numfigs = 0
        self.samplefigs = samplefigs
        if maximiser is None:
            maximiser = NoisyMax()
        self.maximiser = maximiser
        self.samplegrid = samplegrid

    def start_training(self):
        self.data = []
    
    def end_training(self):
        pass

    def train_response(self, response, location, meta):
        self.data.append((response, location, meta))
    
    def test_response(self, response):
        pass
    
    def save(self, filename, datapath=datapath):
        pass
    
    def load(self, filename, datapath=datapath):
        pass
    
    def train(self, responses=None, start_training=True, end_training=True):
        if start_training:
            self.start_training()
        if responses is None:
            responses = self.fixedmodel.isaved_responses()
        for response, location, meta, sound in responses:
            self.train_response(response, location, meta)
        if end_training:
            self.end_training()
        
    def test(self, responses=None, store_sounds=False):
        if responses is None:
            responses = self.fixedmodel.isaved_responses()
        results = []
        for response, location, meta, sound in responses:
            self.numfigs += 1
            if _have_sample_figs and self.numfigs<=self.samplefigs:
                figure()
            self._current_test_location = float(location) # for diagnostics of estimators
            estimate = self.test_response(response)
            if not store_sounds:
                sound = None
            results.append((location, estimate, meta, sound))
        return results

    def plotsamplefig(self, bd, response, *args, **kwds):
        if _have_sample_figs and self.numfigs<=self.samplefigs:
            plot(bd, response, *args, **kwds)
        if self.numfigs>0 and self.numfigs==self.samplefigs:
            show()
            exit()
            
    def __str__(self):
        return '%s(%s,%d)'%(self.__class__.__name__, str(self.maximiser),
                            self.samplegrid)
    
    def __repr__(self):
        return str(self)

class ScikitsLearnEstimator(Estimator):
    '''
    Use scikits.learn regression.
    '''
    def __init__(self, fixedmodel, regressor, normalisation=True):
        Estimator.__init__(self, fixedmodel)
        self.normalisation = normalisation
        self.regressor = regressor
    
    def start_training(self):
        self.data = []

    def train_response(self, response, location, meta):
        if self.normalisation:
            self.data.append((response/amax(response), location))
        else:
            self.data.append((response, location))

    def end_training(self):
        responses, locations = zip(*self.data)
        self.regressor.fit(responses, locations)
    
    def test_response(self, response):
        if self.normalisation:
            return self.regressor.predict([response/amax(response)])[0]
        else:
            return self.regressor.predict([response])[0]
        
    def __str__(self):
        strhash = str(hash(str((self.regressor, self.normalisation))))
        return self.regressor.__class__.__name__+strhash

class NoisyMax(object):
    def setparams(self, estimator, bd):
        self.estimator = estimator
        self.bd = bd
        self.itdmax = itdmax = float(estimator.fixedmodel.space.itd_max)
        if hasattr(self, 'itdmax_extend'):
            itdmax *= self.itdmax_extend
        self.x = linspace(-itdmax, itdmax, estimator.samplegrid)
    def plotsamplefig(self, *args, **kwds):
        plotcmd = kwds.pop('plot', plot)
        if _have_sample_figs and self.estimator.numfigs<=self.estimator.samplefigs:
            plotcmd(*args, **kwds)
            axvline(self.estimator._current_test_location, ls='--', color='m')
    def __call__(self, response):
        bd = self.bd
        i = argmax(response)
        if hasattr(self, 'useblue'):
            self.plotsamplefig(bd, response, '.b')
        else:
            self.plotsamplefig(bd, response, '.k')
        self.plotsamplefig([bd[i]], [response[i]], 'or')
        return bd[i]
    def __repr__(self):
        return str(self)
    def __str__(self):
        return self.__class__.__name__
    
class DiscreteMax(NoisyMax):
    def __init__(self, window_len=None):
        self.window_len = window_len
    def setparams(self, estimator, bd):
        NoisyMax.setparams(self, estimator, bd)
        ubd = unique(bd)
        ubd.sort()
        self.ubd = ubd = hstack((ubd, Inf))
        self.histcount, _ = histogram(bd, ubd)
    def __call__(self, response):
        tot, _ = histogram(self.bd, self.ubd, weights=response)
        dm = tot/self.histcount
        self.plotsamplefig(self.bd, response, '.k')
        self.plotsamplefig(self.ubd[:-1], dm, '-b')
        if self.window_len is not None:
            dm = smooth(dm, window_len=self.window_len)
            self.plotsamplefig(self.ubd[:-1], dm, '-c')
        i = argmax(dm)
        m = self.ubd[i]
        self.plotsamplefig([m], [dm[i]], 'or')
        return m
    def __str__(self):
        return '%s(%s)'%(self.__class__.__name__, self.window_len)

class SmoothedMax(NoisyMax):
    def __init__(self, width, itdmax_extend=None):
        self.width = width
        if itdmax_extend is not None:
            self.itdmax_extend = itdmax_extend
    def setparams(self, estimator, bd):
        NoisyMax.setparams(self, estimator, bd)
        itd = self.x
        itd = reshape(itd, (1, len(itd)))
        bd = reshape(bd, (len(bd), 1))
        weight = exp(-(bd-itd)**2/(2*self.width**2))
        sumweight = sum(weight, axis=0)
        sumweight = reshape(sumweight, (1, len(sumweight)))
        weight /= sumweight
        self.weight = weight
    def compute_smoothed(self, response):
        bd = self.bd
        self.plotsamplefig(bd, response, '.k')
        response = reshape(response, (len(response), 1))
        response = sum(response*self.weight, axis=0)
        self.plotsamplefig(self.x, response, '.')
        return self.x, response
    def __call__(self, response):
        x, response = self.compute_smoothed(response)
        i = argmax(response)
        self.plotsamplefig([x[i]], [response[i]], 'or')
        return x[i]
    def __str__(self):
        s = self.__class__.__name__+'('+str(self.width)
        if hasattr(self, 'itdmax_extend'):
            s += ', '+str(self.itdmax_extend)
        s += ')'
        return s

# Local linear regression smoothed max
class LLRSmoothedMax(SmoothedMax):
    def __init__(self, width, discrete=False, itdmax_extend=None):
        SmoothedMax.__init__(self, width, itdmax_extend)
        self.discrete = discrete
    def setparams(self, estimator, bd):
        NoisyMax.setparams(self, estimator, bd)
        if self.discrete:
            self.x = unique(self.bd)
            self.x.sort()
        U = self.x
        width = self.width
        self.xr = xr = reshape(bd, (len(bd), 1))
        self.u = u = reshape(U, (1, len(U)))
        w = (1-abs((xr-u)/width)**3)*(abs(xr-u)<width)
        w /= reshape(sum(w, axis=0), (1, len(U)))
        self.xbar = sum(xr*w, axis=0)
        self.xxbar = sum(xr**2*w, axis=0)
        self.invwcovx = 1/(self.xxbar-self.xbar**2)
        self.weight = w
    def compute_smoothed(self, response):
        bd = self.bd
        self.plotsamplefig(bd, response, '.k')
        y = response
        # vectorised local linear regression
        xr = self.xr
        yr = reshape(y, (len(y), 1))
        u = self.u
        w = self.weight
        xbar = self.xbar
        ybar = sum(yr*w, axis=0)
        xybar = sum((xr*yr)*w, axis=0)
        xxbar = self.xxbar
        #a = (xybar-xbar*ybar)/(xxbar-xbar**2)
        a = (xybar-xbar*ybar)*self.invwcovx
        b = ybar-a*xbar
        #a = (sum(xr*yr*w, axis=0)-sum(yr*w, axis=0)*sum(xr*w, axis=0))/(sum(xr**2*w, axis=0)-sum(xr*w, axis=0)**2)
        #b = sum(yr*w, axis=0)-a*sum(xr*w, axis=0)
        response = a*self.x+b
        self.plotsamplefig(self.x, response, '-')
        return self.x, response
    def __str__(self):
        s = self.__class__.__name__+'('+str(self.width)+', '+str(self.discrete)
        if hasattr(self, 'itdmax_extend'):
            s += ', '+str(self.itdmax_extend)
        s += ')'
        return s

# Local quadratic regression smoothed max
class LQRSmoothedMax(SmoothedMax):
    def __init__(self, width, discrete=False, itdmax_extend=None):
        SmoothedMax.__init__(self, width, itdmax_extend)
        self.discrete = discrete
    def setparams(self, estimator, bd):
        NoisyMax.setparams(self, estimator, bd)
        if self.discrete:
            self.x = unique(self.bd)
            self.x.sort()
        U = self.x
        width = self.width
        self.xr = xr = reshape(bd, (len(bd), 1))
        self.u = u = reshape(U, (1, len(U)))
        w = (1-abs((xr-u)/width)**3)*(abs(xr-u)<width)
        w /= reshape(sum(w, axis=0), (1, len(U)))
        self.X = X = [sum(xr**k*w, axis=0) for k in xrange(5)]
        self.Dinv = 1.0/(-X[0]*X[2]*X[4]+X[0]*X[3]**2+X[1]**2*X[4]-2*X[1]*X[2]*X[3]+X[2]**3)
        self.weight = w
    def compute_smoothed(self, response):
        bd = self.bd
        self.plotsamplefig(bd, response, '.k')
        y = response
        # vectorised local quadratic regression
        xr = self.xr
        yr = reshape(y, (len(y), 1))
        u = self.u
        w = self.weight
        Y = [sum(xr**k*yr*w, axis=0) for k in xrange(3)]
        X = self.X
        Dinv = self.Dinv

        a = Dinv*(Y[2]*X[1]**2-X[3]*Y[0]*X[1]-X[2]*Y[1]*X[1]+X[2]**2*Y[0]+X[0]*X[3]*Y[1]-X[0]*X[2]*Y[2])
        b = Dinv*(Y[1]*X[2]**2-X[3]*Y[0]*X[2]-X[1]*Y[2]*X[2]+X[1]*X[4]*Y[0]-X[0]*X[4]*Y[1]+X[0]*X[3]*Y[2])
        c = Dinv*(Y[2]*X[2]**2-X[4]*Y[0]*X[2]-X[3]*Y[1]*X[2]+X[3]**2*Y[0]+X[1]*X[4]*Y[1]-X[1]*X[3]*Y[2])

        response = a*self.x**2+b*self.x+c
        self.plotsamplefig(self.x, response, '-')
        return self.x, response
    def __str__(self):
        s = self.__class__.__name__+'('+str(self.width)+', '+str(self.discrete)
        if hasattr(self, 'itdmax_extend'):
            s += ', '+str(self.itdmax_extend)
        s += ')'
        return s
    
class FitMax(NoisyMax):
    def __init__(self, type='poly', degmax=2, algorithm='leastsq',
                 neighbourhood=0.25):
        self.type = type
        self.degmax = degmax
        self.algorithm = algorithm
        self.neighbourhood = neighbourhood
    def fitmax(self, x, y, xout=None):
        self.plotsamplefig(x, y, '.')
        if xout is None:
            xout = self.x
        if self.type=='poly':
            degmax = self.degmax
            polys = [polyfit(x, y, deg) for deg in range(1, degmax+1)]
            err = [sum((polyval(p, x)-y)**2) for p in polys]
            p = polys[argmin(err)]
            yout = polyval(p, xout)
            self.plotsamplefig(xout, yout, '.')
            i = argmax(yout)
            self.plotsamplefig([xout[i]], [yout[i]], 'og')
            return xout[i]
        elif self.type=='gaussian':
            sumy = sum(y)
            y = y/sumy
            f = lambda a, mu, sigma: a*exp(-(x-mu)**2/(2*sigma**2))
            g = lambda p: f(*p)-y
            ymax = amax(y)
            absxmax = amax(abs(self.x))
            a0 = mean(y)
            mu0 = sum(x*y)
            std0 = sqrt(sum(y*(x-mu0)**2))
            if self.algorithm=='leastsq':
                p, _ = leastsq(g, (a0, mu0, std0))
            elif self.algorithm=='fmin_cobyla':
                h = lambda p: sum(g(p)**2)
                cons = [lambda (a, mu, sigma): absxmax-abs(mu),
                        lambda (a, mu, sigma): a,
                        lambda (a, mu, sigma): ymax-a,
                        lambda (a, mu, sigma): sigma,
                        lambda (a, mu, sigma): 4*absxmax-sigma]
                p = optimize.fmin_cobyla(h, (a0, mu0, std0), cons)
            else:
                raise ValueError('Algorithm %s unknown.'%self.algorithm)
            a, mu, sigma = p
            self.plotsamplefig(x, sumy*f(a0, mu, sigma), '.')
            self.plotsamplefig(x, sumy*f(a, mu, sigma), '.')
            mu = clip(mu, -absxmax, absxmax)
            return mu
        else:
            raise ValueError('type %s unknown'%self.type)
    def __call__(self, response):
        bd = self.bd
        i = argmax(response)
        xwidth = amax(bd)-amin(bd)
        I = abs(bd-bd[i])<xwidth*self.neighbourhood
        if sum(I)>5:
            bd = bd[I]
            response = response[I]
        return self.fitmax(bd, response)
    def __str__(self):
        strhash = str(hash(str((self.type, self.degmax, self.algorithm,
                                self.neighbourhood))))
        return self.__class__.__name__+strhash

class FitSmoothedMax(SmoothedMax, FitMax):
    def __init__(self, width, type='poly', degmax=2, algorithm='leastsq',
                 neighbourhood=0.25):
        SmoothedMax.__init__(self, width)
        FitMax.__init__(self, type=type, degmax=degmax, algorithm=algorithm,
                        neighbourhood=neighbourhood)
    def __call__(self, response):
        bd = self.bd
        x, y = self.compute_smoothed(response)
        i = argmax(y)
        xwidth = amax(x)-amin(x)
        I = abs(bd-x[i])<xwidth*self.neighbourhood
        if sum(I)>5:
            bd = bd[I]
            response = response[I]
        return self.fitmax(bd, response,
                           linspace(amin(bd), amax(bd), len(self.x)))
    def __str__(self):
        strhash = str(hash(str((self.type, self.degmax, self.algorithm,
                                self.neighbourhood, self.width))))
        return self.__class__.__name__+strhash

class KDEModeMax(NoisyMax):
    def __init__(self, fitter=None, grid=100):
        if fitter is None:
            fitter = NoisyMax()
        self.fitter = fitter
        self.ygrid = grid
    def setparams(self, estimator, bd):
        NoisyMax.setparams(self, estimator, bd)
        self.fitter.setparams(estimator, self.x)
    def __call__(self, response):
        bd = self.bd
        rmin, rmax = amin(response), amax(response)
        y = linspace(rmin, rmax, self.estimator.samplegrid)
        X, Y = meshgrid(self.x, y)
        shape = X.shape
        X.shape = X.size
        Y.shape = Y.size
        g = gaussian_kde(vstack((self.bd, response)))
        I = g.evaluate(vstack((X, Y)))
        I.shape = shape
        newresponse = y[argmax(I, axis=0)]
        self.plotsamplefig(I, origin='lower left',
                           extent=(-self.itdmax, self.itdmax, rmin, rmax),
                           aspect='auto', plot=imshow)
        self.plotsamplefig(bd, response, '.w')
        self.plotsamplefig(self.x, newresponse, 'ow')
        return self.fitter(newresponse)
    

def phase_multiply(cfs, bds, itdmax):
#    plot(cfs, bds/usecond, '.')
    itdmax = float(itdmax)
    cfout = []
    bdout = []
    I = []
    for i, (cf, bd) in enumerate(zip(cfs, bds)):
        bp = (2*pi*cf*bd)%(2*pi) 
        phimax = 2*pi*cf*itdmax
        n = max(ceil((phimax-bp)/(2*pi)), -floor((-phimax-bp)/(2*pi)))+1
        bp = bp+2*pi*arange(-n, n+1)
        bp = bp[abs(bp)<phimax]
        if len(bp):
            bd = bp/(2*pi*cf)
        else:
            bd = array([bd])
        cfout.extend([cf]*len(bd))
        bdout.extend(bd)
        I.extend([i]*len(bd))
    I = array(I, dtype=int)
    cfout = array(cfout)
    bdout = array(bdout)
#    plot(cfout+1, bdout/usecond, '.')
#    show()
#    exit()
    return cfout, bdout, I


class Jeffress(Estimator):
    def __init__(self, fixedmodel,
                 maximiser=None, samplegrid=400,
                 phasemultiply=False,
                 lesionfunc=None,
                 samplefigs=0,):
        Estimator.__init__(self, fixedmodel,
                           maximiser=maximiser, samplegrid=samplegrid,
                           samplefigs=samplefigs)
        self.phasemultiply = phasemultiply
        bd = self.bd.best_delay
        if phasemultiply:
            cf = self.bd.best_frequencies
            itdmax = self.fixedmodel.space.itd_max
            cf, bd, I = phase_multiply(cf, bd, itdmax)
            self.pm_indices = I
        self.maximiser.setparams(self, bd)
        self.lesionfunc = lesionfunc
        if lesionfunc is not None:
            self.lesionarr = lesionfunc(bd)
    def test_response(self, response):
        if self.lesionfunc is not None:
            response = response*self.lesionarr
        if self.phasemultiply:
            response = response[self.pm_indices]
        return self.maximiser(response)
    def __str__(self):
        s = '%s(%s,%d,%s)'%(self.__class__.__name__, str(self.maximiser),
                            self.samplegrid, str(self.phasemultiply))
        if self.lesionfunc is not None:
            s = s+'_'+str(hash(inspect.getsource(self.lesionfunc)))
        return s

class JeffressAssemblies(Estimator):
    def __init__(self, fixedmodel, bins=None, acc=1e-10):
        Estimator.__init__(self, fixedmodel)
        itdmax = float(self.fixedmodel.space.itd_max)
        bd = self.bd.best_delay
        if isinstance(bins, int):
            bins = linspace(-itdmax-1e-20, itdmax+1e-20, bins+1)
        elif bins is None:
            # In this case, we use one bin per best delay up to fixed accuracy
            bdint = array(bd/acc, dtype=int)
            I = unique(bdint)
            bins = []
            for i in I:
                thisbd = bd[bdint==i]
                bins.extend([amin(thisbd)-1e-20, amax(thisbd)+1e-20])
            bins.sort()
            bins = array(bins)
        self.bins = bins
        self.numinbin, _ = histogram(bd, bins)
        self.numinbin[self.numinbin==0] = 1
        self.binmids = 0.5*(bins[1:]+bins[:-1])
    def test_response(self, response):
        bd = self.bd.best_delay
        R, _ = histogram(bd, self.bins, weights=response)
        R = R/self.numinbin
        return self.binmids[argmax(R)]


class UsesIdealResponses(object):
    def compute_ideal_responses(self):
        if not hasattr(self, 'idealresponses'):
            self.idealresponses = []
            self.locations = []
            for _, location, meta in self.data:
                self.idealresponses.append(meta['noiseless_response'])
                self.locations.append(location)
            self.locations = array(self.locations)
            self.idealresponses = array(self.idealresponses) # shape (numdata, numcells)


class TrainedJeffress(Estimator, UsesIdealResponses):
    def __init__(self, fixedmodel,
                 maximiser=None, samplegrid=400,
                 bdmaximiser=None,
                 phasemultiply=False,
                 samplefigs=0,):
        Estimator.__init__(self, fixedmodel,
                           maximiser=maximiser, samplegrid=samplegrid,
                           samplefigs=samplefigs)
        self.phasemultiply = phasemultiply
        self.bdmaximiser = bdmaximiser
    def end_training(self):
        self.compute_ideal_responses()
        if self.bdmaximiser is None:
            bd = self.locations[argmax(self.idealresponses, axis=0)]
        else:
            bdm = self.bdmaximiser
            if self.samplefigs==0:
                bdm.plotsamplefig = lambda *args, **kwds: None
            self._current_test_location = 0
            bdm.setparams(self, self.locations)
            bd = []
            for i, R in enumerate(self.idealresponses.T):
                bd.append(bdm(R))
                if i==0:
                    bdm.plotsamplefig = lambda *args, **kwds: None
            bd = array(bd)
        if self.phasemultiply:
            cf = self.bd.best_frequencies
            itdmax = self.fixedmodel.space.itd_max
            cf, bd, I = phase_multiply(cf, bd, itdmax)
            self.pm_indices = I
        self.maximiser.setparams(self, bd)
    def test_response(self, response):
        if self.phasemultiply:
            response = response[self.pm_indices]
        return self.maximiser(response)
    def __str__(self):
        s = ''
        if self.phasemultiply is True:
            s += ',phasemultiply=True'
        return '%s(%s,%d,%s%s)'%(self.__class__.__name__, str(self.maximiser),
                               self.samplegrid, str(self.bdmaximiser), s)


class FitDistributionMLE(Estimator, UsesIdealResponses):
    def end_training(self):
        self.compute_ideal_responses()
        kmax = max(amax(self.idealresponses)*2, 10000)
        self.logkfac = gammaln(arange(kmax)+1)
        self.logidealresponses = log(self.idealresponses)
        self.maximiser.setparams(self, self.locations)
    def test_response(self, response):
        response = array(response, dtype=int)
        logkfac = self.logkfac[response]
        k = reshape(response, (1, len(response)))
        response = sum(k*self.logidealresponses-logkfac-self.idealresponses, axis=1)
        bd = self.locations
        return self.maximiser(response)


class InequalityPatternMatch(Estimator, UsesIdealResponses):
    def __init__(self, fixedmodel, num_inequalities=10000,
                 maximiser=None, samplegrid=400,
                 ):
        Estimator.__init__(self, fixedmodel,
                           maximiser=maximiser, samplegrid=samplegrid,
                           samplefigs=0)
        self.num_inequalities = num_inequalities
        
    def end_training(self):
        self.compute_ideal_responses()
        self.maximiser.setparams(self, self.locations)
        
    def ipm_response(self, response):
        # response has shape (numcells,)
        # idealresponses has shape (numdata, numcells)
        # output should have shape (numdata,)
        idealresponses = self.idealresponses
        numdata, numcells = idealresponses.shape
        output = zeros(numdata)
        ns = {'response': response, 'idealresponses':idealresponses, 'numdata':numdata, 'numcells':numcells,
              'output': output, 'num_inequalities': self.num_inequalities}
        weave.inline('''
        // i = data
        // j = cell LHS
        // k = cell RHS
        for(int i=0; i<numdata; i++)
        {
            // this will be the number of pairs with stored pattern equal to response pattern
            int num_equal = 0;
            srand(32943249);
            for(int s=0; s<num_inequalities; s++)
            {
                int j = rand() % numcells;
                int k = rand() % numcells;
            //for(int j=0; j<numcells; j++)
            //{
            //    for(int k=0; k<numcells; k++)
            //    {
                    // compute inequality for response
                    bool ineq_response = response[j]<response[k];
                    // compute inequality for stored patterns
                    bool stored_response = idealresponses[j+i*numcells]<idealresponses[k+i*numcells];
                    num_equal += (int)(ineq_response==stored_response);
            //    }
            //}
            }
            output[i] = (double)num_equal;
        }
        ''', ns.keys(), ns, compiler='gcc', extra_compile_args=['-O3', '-ffast-math', '-march=native'])
        return output
#        response = reshape(response, (1, len(response)))
#        response = sum(response*self.idealresponses, axis=1)
#        return response

    def test_response(self, response):
        return self.maximiser(self.ipm_response(response))
    
    def __str__(self):
        return self.__class__.__name__+str(self.maximiser)+'_'+str(self.num_inequalities)


def pm_default_normfunc(x):
    return sqrt(sum(x**2)) 
class PatternMatch(Estimator, UsesIdealResponses):
    def __init__(self, fixedmodel,
                 maximiser=None, samplegrid=400,
                 pools=None,
                 normalise=False,
                 normalise_response=False,
                 normalise_banded=False,
                 debugir=False,
                 use_analytic=False,
                 irfit=False,
                 normfunc=pm_default_normfunc,
                 raise_to_power=1,
                 norm_p=2,
                 log_pattern=False,
                 apply_func=None,
                 samplefigs=0,
                 lesionfunc=None,
                 ):
        Estimator.__init__(self, fixedmodel,
                           maximiser=maximiser, samplegrid=samplegrid,
                           samplefigs=samplefigs)
        self.pools = pools
        self.debugir = debugir
        self.normalise = normalise
        self.normalise_response = normalise_response
        self.normalise_banded = normalise_banded
        self.use_analytic = use_analytic
        self.irfit = irfit
        self.normfunc = normfunc
        self.raise_to_power = raise_to_power
        self.log_pattern = log_pattern
        self.lesionfunc = lesionfunc
        if lesionfunc is not None:
            bd = self.bd.best_delay
            self.lesionarr = lesionfunc(bd)
        if apply_func is None:
            apply_func = lambda x: x
        else:
            if '\n' in apply_func:
                ns = {}
                exec apply_func in globals(), ns
                apply_func = ns['apply_func']
            else:
                apply_func = eval('lambda x: '+apply_func, globals())
        self.apply_func = apply_func
        self.norm_p = norm_p
    def end_training(self):
        p = self.norm_p
        self.compute_ideal_responses()
        if self.lesionfunc is not None:            
            lesion = reshape(self.lesionarr, (1, -1))
            self.idealresponses *= lesion
        if self.raise_to_power!=1:
            self.idealresponses **= self.raise_to_power
        IRdiff = zeros_like(self.idealresponses)
        if self.normalise and self.bd.cfrepeat>1:
            cfrepeat = self.bd.cfrepeat
            bdp = self.bd.best_delay+repeat(arange(len(self.bd.best_delay)/cfrepeat), cfrepeat)
            J = argsort(bdp)
            I = argsort(self.locations)
            IR = self.idealresponses[I, :]
            IR = IR[:, J]
            IRnew = reshape(IR, (IR.shape[0], IR.shape[1]/cfrepeat, cfrepeat, 1))
            IRnew2 = swapaxes(IRnew, 0, 3)
            dp = mean(sum(IRnew*IRnew2, axis=2), axis=2)
            IRnew /= reshape(dp, (IR.shape[0], IR.shape[1]/cfrepeat, 1, 1))
            IR = reshape(IRnew, IR.shape)
            IR = IR[argsort(I), :]
            IR = IR[:, argsort(J)]
            self.idealresponses = IR
        if self.irfit:
            cf = self.bd.best_frequencies
            bd = self.bd.best_delay
            loc = self.locations
            IR = self.idealresponses
            for j in xrange(IR.shape[1]):
                IR[:, j] = polyval(polyfit(loc, IR[:, j], self.irfit), loc)
        if self.use_analytic:
            cf = self.bd.best_frequencies
            bd = self.bd.best_delay
            loc = self.locations
            numdata = len(loc)
            numcells = len(bd)
            bd = reshape(bd, (1, numcells))
            cf = reshape(cf, (1, numcells))
            loc = reshape(loc, (numdata, 1))
            IR = cos(pi*cf*(loc-bd))**2
            IRdiff = self.idealresponses
            self.idealresponses = IR
        if self.pools is not None:
            IR = self.idealresponses
            self.idealresponses = hstack(tuple(reshape(sum(IR[:, pool], axis=1), (IR.shape[0], 1)) for pool in self.pools))
        if self.debugir:
            I = argsort(self.locations)
            IR = self.idealresponses
            figure()
            IRphi = reshape(IR-IRdiff, IR.shape+(1,)) # phi, i, (theta)
            IRtheta = swapaxes(IRphi, 0, 2)        # (phi), i, theta
            ER = sum(IRtheta*IRphi, axis=1)/sqrt(sum(IRphi**2, axis=1))
            Rvar = sum(IRtheta*IRphi**2, axis=1)/sum(IRphi**2, axis=1)
            subplot(221)
            imshow(ER)
            colorbar()
            subplot(222)
            imshow(sqrt(Rvar))
            colorbar()
            subplot(223)
            plot(self.locations[I], self.idealresponses[I, :]-IRdiff[I,:])
        if self.normalise_banded:
            cf = self.bd.best_frequencies
            bandsize = self.normalise_banded
            x = self.idealresponses
            J = argsort(cf)
            x = x[:, J]
            for i in xrange(0, x.shape[1], bandsize):
                nf = ((sum(x[:, i:i+bandsize]**p, axis=1))**(1.0/p))[:, newaxis]
                nf[nf<1e-10] = 1
                x[:, i:i+bandsize] /= nf
            x /= ((sum(x**p, axis=1))**(1.0/p))[:, newaxis]
            x = x[:, argsort(J)]
            self.idealresponses = x
        else:
            self.idealresponses /= reshape(((sum(self.idealresponses**p, axis=1)**(1.0/p))),
                                           (len(self.data), 1))
        if self.log_pattern:
            self.idealresponses = log(self.idealresponses)
        if self.debugir:
            subplot(224)
            plot(self.locations[I], self.idealresponses[I, :])
            show()
            exit()
        self.maximiser.setparams(self, self.locations)
        
    def pm_response(self, response):
        response = self.apply_func(response)
        if self.lesionfunc is not None:
            response = response*self.lesionarr
        if self.raise_to_power!=1:
            response = response**self.raise_to_power
        if self.normalise_response and self.bd.cfrepeat>1:
            cfrepeat = self.bd.cfrepeat
            bdp = self.bd.best_delay+repeat(arange(len(self.bd.best_delay)/cfrepeat), cfrepeat)
            J = argsort(bdp)
            response = response[J]
            for j in xrange(len(J)/cfrepeat):
                response[j*cfrepeat:(j+1)*cfrepeat] /= self.normfunc(response[j*cfrepeat:(j+1)*cfrepeat])
            response = response[argsort(J)]
        if self.pools is not None:
            response = hstack(tuple(sum(response[pool]) for pool in self.pools))
        response = reshape(response, (1, len(response)))
        response = sum(response*self.idealresponses, axis=1)
        return response

    def test_response(self, response):
        return self.maximiser(self.pm_response(response))
    
    def __str__(self):
        props = (self.samplegrid, self.pools, self.normalise,
                 self.normalise_response, self.use_analytic)
        if self.lesionfunc is not None:
            props = props+('lesionfunc=\'\'\''+inspect.getsource(self.lesionfunc)+'\'\'\'',)
        if self.raise_to_power!=1:
            props = props+('raise_to_power=%s'%self.raise_to_power,)
        if self.log_pattern:
            props = props+('log_pattern',)
        if self.normalise_banded:
            props = props+('normalise_banded=%s'%self.normalise_banded, 'v1')
        if self.norm_p!=2:
            props = props+('norm_p=%d'%self.norm_p,)
        strhash = str(hash(str(props)))
        return '%s(%s)_%s'%(self.__class__.__name__, str(self.maximiser),
                            strhash)
        

class PatternMatchMatch(PatternMatch):
    def __init__(self, *args, **kwds):
        PatternMatch.__init__(self, *args, **kwds)
        self._args = args
        self._kwds = kwds
    def end_training(self):
        PatternMatch.end_training(self)
        self.matcher = PatternMatch(*self._args, **self._kwds)
        self.matcher.start_training()
        for response, location, meta in self.data:
            response = self.pm_response(response)
            meta = copy(meta)
            meta['noiseless_response'] = response
            self.matcher.data.append((response, location, meta))
        self.matcher.end_training()
        self.matcher.maximiser.useblue = True
    def test_response(self, response):
        self.matcher.numfigs = self.numfigs
        response = self.pm_response(response)
        return self.matcher.test_response(response)


class TwoChannelPatternMatch(PatternMatch):
    def __init__(self, fixedmodel,
                 maximiser=None, samplegrid=400,
                 normalise=True,
                 samplefigs=0, debugir=False):
        bd = fixedmodel.binauraldistribution
        self.left_channel_indices, = (bd.best_delay<0).nonzero()
        self.right_channel_indices, = (bd.best_delay>0).nonzero()
        pools = [self.left_channel_indices, self.right_channel_indices]
        PatternMatch.__init__(self, fixedmodel,
                              maximiser=maximiser, samplegrid=samplegrid,
                              pools=pools, normalise=normalise,
                              samplefigs=samplefigs, debugir=debugir)
   
def grid_to_str(grid):
    if isinstance(grid, ndarray):
        strgrid = hashlib.sha1(grid.view(uint8)).hexdigest()
    else:
        strgrid = str(grid)
    return strgrid

class NoisyInverseCurveFit(object):
    '''
    These classes are used to find the inverse of a noisy curve (x, y)
    
    Issues to note are that you can get bias at the edges of the range.
    '''
    def __init__(self, grid=False):
        self.grid = grid
    def setgrid(self, xmin, xmax, xn):
        if self.grid is True:
            self.x = linspace(xmin, xmax, xn)
        else:
            self.x = self.grid
    def train(self, x, y):
        if self.grid is False:
            self.x = x
        self.y = y
    def __call__(self, y):
        pass
    def __str__(self):
        strgrid = grid_to_str(self.grid)
        return self.__class__.__name__+'('+strgrid+')'
    def __repr__(self):
        return str(self)
    
class Closest(NoisyInverseCurveFit):
    def __call__(self, y):
        return self.x[argmin(abs(y-self.y))]
    
class SmoothedClosest(Closest):
    def __init__(self, width, grid=True):
        self.width = width
        self.grid = grid
    def train(self, x, y):
        self.y = gaussian_smooth(x, y, self.width, xout=self.x)
    def __str__(self):
        strgrid = grid_to_str(self.grid)
        return '%s(%s,%s)'%(self.__class__.__name__, str(self.width),
                            strgrid)

class PolyClosest(Closest):
    def __init__(self, degmax, grid=True):
        self.degmax = degmax
        self.grid = grid
    def train(self, x, y):
        degmax = self.degmax
        polys = [polyfit(x, y, deg) for deg in range(1, degmax+1)]
        err = [sum((polyval(p, x)-y)**2) for p in polys]
        p = polys[argmin(err)]
        self.y = polyval(p, self.x)
    def __str__(self):
        return '%s(%s,%s)'%(self.__class__.__name__, str(self.degmax),
                            grid_to_str(self.grid))


class OneDimensionalFunctionEstimator(Estimator, UsesIdealResponses):
    def __init__(self, fixedmodel, fitter=None, samplegrid=400,
                 itdmax_extend=1.0,
                 samplefigs=0, debugfit=False,
                 ):
        Estimator.__init__(self, fixedmodel, samplefigs=samplefigs,
                           samplegrid=samplegrid)
        self.debugfit = debugfit
        if fitter is None:
            fitter = Closest()
        self.fitter = fitter
        itdmax = float(self.fixedmodel.space.itd_max)*itdmax_extend
        self.itdmax_extend = itdmax_extend
        fitter.setgrid(-itdmax, itdmax, self.samplegrid)
        
    def func(self, response):
        raise NotImplementedError
        
    def end_training(self):
        self.compute_ideal_responses()
        self.values = array([self.func(R) for R in self.idealresponses])
        if _have_sample_figs and self.debugfit:
            figure()
            plot(self.locations, self.values, '.k')
        self.fitter.train(self.locations, self.values)
        if _have_sample_figs and self.debugfit:
            plot(self.fitter.x, self.fitter.y, '-or')
        
    def test_response(self, response):
        value = self.func(response)
        loc = self.fitter(value)
        self.plotsamplefig(self.locations, self.values, '.k')
        self.plotsamplefig([loc], [value], 'or')
        return loc

    def __str__(self):
        s = '%s(%s,%d,%s)'%(self.__class__.__name__, str(self.fitter),
                            self.samplegrid, str(self.itdmax_extend))
        return s


class LinearOneDimensionalFunctionEstimator(OneDimensionalFunctionEstimator):
    def __init__(self, fixedmodel, coefficients, fitter=None, samplegrid=400,
                 itdmax_extend=1.0,
                 samplefigs=0, debugfit=False):
        OneDimensionalFunctionEstimator.__init__(self, fixedmodel,
                                    fitter=fitter, samplegrid=samplegrid,
                                    itdmax_extend=itdmax_extend,
                                    samplefigs=samplefigs, debugfit=debugfit)
        self.coefficients = coefficients
    def func(self, response):
        sumresponse = float(sum(response))
        if sumresponse<1e-10:
            sumresponse = 1e-10
        return sum(response*self.coefficients)/sumresponse


class TwoChannel(LinearOneDimensionalFunctionEstimator):
    def __init__(self, fixedmodel, fitter=None, difference=False,
                 samplegrid=400,
                 itdmax_extend=1.0,
                 samplefigs=0, debugfit=False):
        self.difference = difference
        if difference:
            coefficients = (fixedmodel.binauraldistribution.best_delay<0)*2.0-1.0
        else:
            coefficients = (fixedmodel.binauraldistribution.best_delay<0)*1.0
        LinearOneDimensionalFunctionEstimator.__init__(self, fixedmodel,
                                    coefficients,
                                    fitter=fitter, samplegrid=samplegrid,
                                    itdmax_extend=itdmax_extend,
                                    samplefigs=samplefigs, debugfit=debugfit)


class TwoChannelCrossFrequency(LinearOneDimensionalFunctionEstimator):
    def __init__(self, fixedmodel, fitter=None,
                 samplegrid=400,
                 itdmax_extend=1.0,
                 samplefigs=0, debugfit=False):
        bd = fixedmodel.binauraldistribution.best_delay
        cf = array(fixedmodel.binauraldistribution.best_frequencies)
        coefficients = ((bd<0)*2.0-1.0)*(1.0/cf)
        LinearOneDimensionalFunctionEstimator.__init__(self, fixedmodel,
                                    coefficients,
                                    fitter=fitter, samplegrid=samplegrid,
                                    itdmax_extend=itdmax_extend,
                                    samplefigs=samplefigs, debugfit=debugfit)


class TwoChannelCrossFrequencyAlt(OneDimensionalFunctionEstimator):
    def __init__(self, fixedmodel, fitter=None, samplegrid=400,
                 itdmax_extend=1.0,
                 samplefigs=0, debugfit=False):
        OneDimensionalFunctionEstimator.__init__(self, fixedmodel,
                                    fitter=fitter, samplegrid=samplegrid,
                                    itdmax_extend=itdmax_extend,
                                    samplefigs=samplefigs, debugfit=debugfit)
        bd = fixedmodel.binauraldistribution.best_delay
        cf = array(fixedmodel.binauraldistribution.best_frequencies)
        self.numerator_coefficients = ((bd<0)*2.0-1.0)
        self.denominator_coefficients = cf
    def func(self, response):
        numer = sum(response*self.numerator_coefficients)
        denom = sum(response*self.denominator_coefficients)
        if denom<1e-10:
            denom = 1e-10
        return numer/denom
#        sumresponse = float(sum(response))
#        if sumresponse<1e-10:
#            sumresponse = 1e-10
#        return sum(response*self.coefficients)/sumresponse
    

class Vector(LinearOneDimensionalFunctionEstimator):
    def __init__(self, fixedmodel, fitter=None,
                 samplegrid=400,
                 itdmax_extend=1.0,
                 samplefigs=0, debugfit=False):
        coefficients = fixedmodel.binauraldistribution.best_delay
        LinearOneDimensionalFunctionEstimator.__init__(self, fixedmodel,
                                    coefficients,
                                    fitter=fitter, samplegrid=samplegrid,
                                    itdmax_extend=itdmax_extend,
                                    samplefigs=samplefigs, debugfit=debugfit)


class TwoChannelBanded(Estimator, UsesIdealResponses):
    def __init__(self, fixedmodel, bandsize=40,
                 samplefigs=0):
        Estimator.__init__(self, fixedmodel, samplefigs=samplefigs)
        self.bandsize = bandsize
        
    def end_training(self):
        self.compute_ideal_responses()
        cf = self.fixedmodel.centre_frequencies[:len(self.fixedmodel.centre_frequencies)/2]
        bestdelay = self.bd.best_delay
        I = argsort(cf)
        cfI = cf[I]
        bestdelayI = bestdelay[I]
        self.I = I
        self.cf = cf
        self.cfI = cfI
        self.bestdelay = bestdelay
        self.bestdelayI = bestdelayI
        IR = self.idealresponses # (numdata, numcells)
        IR = IR[:, I]
        if IR.shape[1]%self.bandsize!=0:
            extendresponses = self.bandsize*(1+IR.shape[1]//self.bandsize)-IR.shape[1]
            self.extendresponses = extendresponses
            IR = concatenate((IR,
                              zeros((IR.shape[0],
                                     extendresponses))),
                                    axis=1)
        else:
            self.extendresponses = 0
        IR.shape = (IR.shape[0], -1, self.bandsize)
        L = self.locations
        coeffs = (bestdelayI<0)*2-1
        coeffs = concatenate((coeffs,
                              zeros(self.extendresponses)),
                              axis=0)
        self.coefficients = coeffs
        C = coeffs.reshape(1, -1, self.bandsize)
        IRC = sum(IR*C, axis=2)/sum(IR, axis=2) # (numdata, numbands)
        A = []
        B = []
        for i in xrange(IRC.shape[1]):
            R = IRC[:, i]
            a, b = polyfit(L, R, 1)
            A.append(a)
            B.append(b)
        A = array(A)
        B = array(B)
        self.a = A
        self.b = B
        
    def test_response(self, response):
        cf, I, cfI = self.cf, self.I, self.cfI
        bestdelay = self.bestdelay
        bestdelayI = self.bestdelayI
        r = response[I]*1.0
        r = concatenate((r, zeros(self.extendresponses)), axis=0)
        rs = r*self.coefficients
        r.shape = (-1, self.bandsize)
        rs.shape = (-1, self.bandsize)
        sumr1 = sum(r, axis=1)
        loc = sum(sumr1*(sum(rs, axis=1)/sumr1-self.b)/self.a)/sum(r)
        return loc

    def __str__(self):
        s = '%s(%d,v6)'%(self.__class__.__name__, self.bandsize)
        return s


class TwoChannelBandedAlt(LinearOneDimensionalFunctionEstimator):
    def __init__(self, fixedmodel, fitter=None,
                 bandsize=40,
                 samplegrid=400,
                 itdmax_extend=1.0,
                 samplefigs=0, debugfit=False):
        self.bandsize = bandsize
        coefficients = (fixedmodel.binauraldistribution.best_delay<0)*2.0-1.0
        LinearOneDimensionalFunctionEstimator.__init__(self, fixedmodel,
                                    coefficients,
                                    fitter=fitter, samplegrid=samplegrid,
                                    itdmax_extend=itdmax_extend,
                                    samplefigs=samplefigs, debugfit=debugfit)

    def end_training(self):
        # set coefficients
        self.compute_ideal_responses()
        cf = self.fixedmodel.centre_frequencies[:len(self.fixedmodel.centre_frequencies)/2]
        bestdelay = self.bd.best_delay
        I = argsort(cf)
        cfI = cf[I]
        bestdelayI = bestdelay[I]
        self.I = I
        self.cf = cf
        self.cfI = cfI
        self.bestdelay = bestdelay
        self.bestdelayI = bestdelayI
        IR = self.idealresponses # (numdata, numcells)
        IR = IR[:, I]
        if IR.shape[1]%self.bandsize!=0:
            extendresponses = self.bandsize*(1+IR.shape[1]//self.bandsize)-IR.shape[1]
            self.extendresponses = extendresponses
            IR = concatenate((IR,
                              zeros((IR.shape[0],
                                     extendresponses))),
                                    axis=1)
        else:
            self.extendresponses = 0
        IR.shape = (IR.shape[0], -1, self.bandsize)
        L = self.locations
        coeffs = (bestdelayI<0)*2-1
        coeffs = concatenate((coeffs,
                              zeros(self.extendresponses)),
                              axis=0)
        C = coeffs.reshape(1, -1, self.bandsize)
        IRC = sum(IR*C, axis=2)/sum(IR, axis=2) # (numdata, numbands)
        A = []
        B = []
        for i in xrange(IRC.shape[1]):
            R = IRC[:, i]
            a, b = polyfit(L, R, 1)
            A.append(a)
            B.append(b)
        A = array(A)
        B = array(B)
        A = repeat(A, self.bandsize)
        B = repeat(B, self.bandsize)
        II = argsort(I)
        A = A[II]
        B = B[II]
        self.coefficients = -(self.coefficients-B)/A
        LinearOneDimensionalFunctionEstimator.end_training(self)

    def __str__(self):
        return LinearOneDimensionalFunctionEstimator.__str__(self)+'(bandsize=%d)'%self.bandsize
