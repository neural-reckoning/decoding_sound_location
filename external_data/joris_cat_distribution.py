from base import *
from bd_distribution import *

def load_joris_f2():
    '''
    Explanation of the matlab file JorisF2.mat giving the data from Figure 2 in:
    
    Joris, P., Van de Sande, B., Louage, D., van der Heijden, M. (2006).
    Binaural and cochlear disparities. Proceedings of the National Academy of
    Sciences of the United States of America, 103 (34), 12917-22.
    
    With the struct_as_record=False and squeeze_me=True options for loadmat,
    the object data is a struct with the following relevant attributes:
    
    ['CF_BDdif',         # length 144
     'CF_BDdif_2',
     'CF_BDdif_3',
     'CF_BDntd',         # length 47
     'CF_BDntd_5',
     'CF_BDntd_6',
     'CFman_BDdif',      # length 9
     'CFman_BDdif_8',
     'CFman_BDdif_9',
     'CFman_BDntd',      # length 9
     'CFman_BDntd_11',
     'CFman_BDntd_12',
     'DFdif_BDdif',      # length 9
     'DFdif_BDdif_14',
     'DFdif_BDdif_15',
     'DFntd_BDntd',      # length 1 (forced)
     'DFntd_BDntd_17',
     'DFntd_BDntd_18']
     
    The names without _number appended are the main peaks and the others are
    the lower and upper sidepeaks.
    
    Each data.attr is an array (we have to force this in the case of DFntd_BDntd
    because it is length 1 which means that it is squeezed out by
    squeeze_me=True).
    
    Each item of the array has one of the following sets of attributes:
    
    diff.bestitd       where difcor was done
    diff.secpeaks
    thr.cf             where CF was available
    nitdp.fitbestitd   where difcor not available
    nitdp.secpeaks
    thr.cfman          not sure?
    diff.fft.df        dominant frequency = peak of FFT of difcor
    
    Each value is either a float, int or occasionally a nan if the value wasn't
    available.
    
    We return a simplified data structure which is just:
    
    (f, bitd, lowpeak, highpeak)
    
    Each being an array of the same length. In this,
    CF and DF are treated as equivalent, and different sources are joined
    together. Note that nans are not filtered, and indicate missing values.
    '''
    fnamebase = datapath+'/optimal_localisation_decoding/JorisF2'
    # Check if loaded and decoded values have been cached, and if so load them
    if os.path.exists(fnamebase+'.pkl'):
        return pickle.load(open(fnamebase+'.pkl', 'rb'))
    fname = fnamebase+'.mat'
    matfile = loadmat(fname, struct_as_record=False, squeeze_me=True)
    data = matfile['T']
    # Force DFntd_BDntd* to be a list like the other attributes
    for attr in ['DFntd_BDntd'+ext for ext in ['', '_17', '_18']]:
        setattr(data, attr, [getattr(data, attr)])
    # each set of variables gives the best peak, low peak and high peaks
    attrsets = [('CF_BDdif', '_2', '_3'),
                ('CF_BDntd', '_5', '_6'),
                ('CFman_BDdif', '_8', '_9'),
                ('CFman_BDntd', '_11', '_12'),
                ('DFdif_BDdif', '_14', '_15'),
                ('DFntd_BDntd', '_18', '_18')
                ]
    # search for F/BITD values in this order
    f_attrs = ['thr.cf', 'thr.cfman', 'diff.fft.df', 'nitdp.fft.df']
    bitd_attrs = ['diff.bestitd', 'diff.secpeaks', 'nitdp.fitbestitd', 'nitdp.secpeaks']
    all_f = zeros(0)
    all_bitd = zeros(0)
    all_lowpeak = zeros(0)
    all_highpeak = zeros(0)
    for baseattr, ext1, ext2 in attrsets:
        curf = []
        curbitd = []
        # this loops through best/low/high peaks
        for attr in [baseattr+ext for ext in ['', ext1, ext2]]:
            var = getattr(data, attr)
            N = len(var)
            f = []
            bitd = []
            for i in xrange(N):
                # This searches the object o=var[i] for values such as
                # o.thr.cf, o.thr.cfman, etc. until it finds an existing,
                # non-nan value
                o = var[i]
                for fa in f_attrs:
                    try:
                        val = eval('o.'+fa)
                        if isnan(val):
                            continue
                        f.append(val)
                        break
                    except:
                        pass
                else:
                    f.append(nan)
                # similarly, but this time for the BITD values rather than the
                # CF/DF values
                for ba in bitd_attrs:
                    try:
                        val = eval('o.'+ba)
                        if isnan(val):
                            continue
                        bitd.append(val)
                        break
                    except:
                        pass
                else:
                    bitd.append(nan)
            curf.append(array(f, dtype=float))
            curbitd.append(array(bitd, dtype=float))
        curf = array(curf) # shape (3, numpoints)
        curbitd = array(curbitd)
        # if BITD is a nan, then the corresponding f will be a nan too, but we
        # want to use non-nan values for f if they are available, and only
        # signal invalid values with the BITDs, so we copy the f value for
        # the lowpeak/highpeak to the main peak column if they are not nans
        for i in [1, 2]:
            I = logical_and(isnan(curf[0, :]), -isnan(curf[i, :]))
            curf[0, I] = curf[i, I]
        all_f = hstack((all_f, curf[0, :]))
        all_bitd = hstack((all_bitd, curbitd[0, :]))
        all_lowpeak = hstack((all_lowpeak, curbitd[1, :]))
        all_highpeak = hstack((all_highpeak, curbitd[2, :]))
    # values are given in Hz
    all_f = array(all_f)*Hz
    # values are given in ms
    all_bitd = array(all_bitd)*ms
    all_lowpeak = array(all_lowpeak)*ms
    all_highpeak = array(all_highpeak)*ms
    return_vals = (all_f, all_bitd, all_lowpeak, all_highpeak)
    # Cache decoded values
    pickle.dump(return_vals, open(fnamebase+'.pkl', 'wb'), -1)
    return return_vals

class JorisCatDistribution(BDDistribution):
    def __init__(self, frange=None, symmetry='swap'):
        f, bitd, lowpeak, highpeak = load_joris_f2()
        BDDistribution.__init__(self, f, bitd, frange=frange, symmetry=symmetry)

if __name__=='__main__':
    if 0:
        joris = JorisCatDistribution(symmetry='swap')
        f, bitd = joris.sample_f_bitd(1000)
        plot(f, bitd, '.')
        f = repeat(joris.sample_f(20), 50)
        #f = repeat([500, 1000, 1500, 200, 2500, 3000], 100)
        bitd = joris.sample_bitd(f)
        plot(f, bitd, '.')
        show()
    if 1:
        # JORIS PNAS FIGURE 2 DATA
        f, bitd, lowpeak, highpeak = load_joris_f2()
#        #plot(f, 2*pi*f*bitd, '.')
#        #plot(f, (2*pi*f*bitd+pi)%(2*pi)-pi, '.')
#        phi = (2*pi*f*bitd+pi)%(2*pi)-pi
#        phi = 2*pi*f*bitd
#        #phi=hstack((phi, -phi))
#        #f = hstack((f,f))
#        phi1=phi[f<1000]
#        phi1 = phi1[-isnan(phi1)] 
#        phi2=phi[f>=1000]
#        phi2 = phi2[-isnan(phi2)]
#        subplot(211) 
#        hist(phi1,20)
#        subplot(212) 
#        hist(phi2,20)
#        show()
#        exit()
        figure()
        subplot(331)
        print 'Number of data points (should be 219=162+57):', len(f)
        plot(f, bitd/ms, 'ob')
        plot(f, lowpeak/ms, 'vg')
        plot(f, highpeak/ms, '^r')
        fmin, fmax = xlim()
        tmin, tmax = ylim()
        fill([fmin, fmax, fmax, fmin], [0.4, 0.4, tmax, tmax], color=(0.8, 0.8, 0.8))
        fill([fmin, fmax, fmax, fmin], [-0.4, -0.4, tmin, tmin], color=(0.8, 0.8, 0.8))
        frange = linspace(fmin+1, fmax, 100)
        plot(frange, 1.0/(2.0*frange)/ms, '-k')
        plot(frange, -1.0/(2.0*frange)/ms, '-k')
        axhline(y=0, color='k')
        ylim(tmin, tmax)
        ylabel('PEAK DELAY (ms)')
        xlabel('CHARACTERISTIC FREQUENCY (Hz)')
        title('Joris PNAS Figure 2')
        # Distribution of peak ITD
        I = logical_and(-isnan(bitd), -isnan(f))
        f = f[I]
        bitd = bitd[I]
        print 'Number of BITD points:', len(f)
        f = hstack((f, f))
        bitd = hstack((bitd, -bitd))
        frange = linspace(fmin, fmax, 100)
        trange = linspace(tmin*ms, tmax*ms, 100)
        subplot(332)
        gf = gaussian_kde(f)
        fprob = gf.evaluate(frange)
        plot(frange, fprob)
        xlabel('CF')
        ylabel('Prob')
        subplot(333)
        prange = linspace(-2*pi, 2*pi, 100)
        gp = gaussian_kde(2*pi*f*bitd)
        pprob = gp.evaluate(prange)
        plot(prange*180/pi, pprob)
        xlabel('BIPD')
        ylabel('Prob')
        
        pf, pt = meshgrid(frange, trange)
        s = pf.shape
        pf.shape = pf.size
        pt.shape = pt.size
        pp = 2*pi*pf*pt
        # Best phase based probability
        subplot(334)
        img = gp.evaluate(pp)
        img.shape = s
        imshow(img, origin='lower left',
               extent=(fmin, fmax, tmin, tmax),
               aspect='auto')
        plot(f, bitd/ms, '.', color=(1,1,1))
        title('BIPD')
        
        # Joint F/ITD based probability
        subplot(335)
        data = vstack((f, bitd))
        g = gaussian_kde(data)
        img = g.evaluate(vstack((pf, pt)))
        img.shape = s
        imshow(img, origin='lower left',
               extent=(fmin, fmax, tmin, tmax),
               aspect='auto')
        plot(f, bitd/ms, '.', color=(1,1,1))
        title('CF/BITD')
        subplot(338)
        img = img/amax(img, axis=0)
        imshow(img, origin='lower left',
               extent=(fmin, fmax, tmin, tmax),
               aspect='auto')
        plot(f, bitd/ms, '.', color=(1,1,1))
        title('CF/BITD conditional on CF')
        gfitd = g
        
        # Joint F/IPD based probability
        subplot(336)
        data = vstack((f, 2*pi*f*bitd))
        g = gaussian_kde(data)
        img = g.evaluate(vstack((pf, pp)))
        img.shape = s
        imshow(img, origin='lower left',
               extent=(fmin, fmax, tmin, tmax),
               aspect='auto')
        plot(f, bitd/ms, '.', color=(1,1,1))
        title('CF/BIPD')
        subplot(339)
        img = img/amax(img, axis=0)
        imshow(img, origin='lower left',
               extent=(fmin, fmax, tmin, tmax),
               aspect='auto')
        plot(f, bitd/ms, '.', color=(1,1,1))
        title('CF/BIPD conditional on CF')
        gfipd = g

        # Resampled data
        subplot(337)
        plot(f, bitd/ms, '.', label='Original')
        rf, rt = gfitd.resample(len(f))
        plot(rf, rt/ms, '.', label='CF/BITD')
        rf, rp = gfipd.resample(len(f))
        rt = rp/(2*pi*rf)
        plot(rf, rt/ms, '.', label='CF/BIPD')
        ylim(-3, 3)
        legend(loc='lower right')
        title('Resampled data')
        show()
