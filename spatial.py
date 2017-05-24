from shared import *
from sphericalhead import *
from external_data.tollin_cat import *
from external_data.wagner_owl_hrtfs import *
from estimation import phase_multiply
from scipy.signal import hann
import hashlib
sys.path.append(os.path.normpath(os.path.join(os.path.split(__file__)[0], '../spatial/')))


class Spatial(object):
    def __init__(self, itd_max, angular=False, fractional_itd=True,
                 fractional_delay_filterlength=2048
                 ):
        self.itd_max = itd_max
        self.itd_max_int = int(itd_max*samplerate)
        self.angular = angular
        self.fractional_itd = fractional_itd
        self.fractional_delay_filterlength = fractional_delay_filterlength
        if angular is True and fractional_itd is False:
            raise ValueError('If angular sampling is used, fractional_itd must be True.')

    def get_random_itd(self):
        if self.angular:
            theta = rand()*pi-pi/2
            itd = sin(azim)*self.itd_max
        elif self.fractional_itd:
            itd = rand()*2*self.itd_max-self.itd_max
        else:
            itd_int = randint(self.itd_max_int*2+1)-self.itd_max_int
            itd = itd_int/samplerate
        return itd
    
    def apply(self, sound, itd):
        right = sound
        if itd<0:
            itd = -itd
            swap = True
        else:
            swap = False
        left = sound.shifted(itd, fractional=self.fractional_itd,
                             filter_length=self.fractional_delay_filterlength)
        if swap:
            left, right = right, left
        nsamples = max(left.nsamples, right.nsamples)
        return Sound([left[:nsamples], right[:nsamples]])
    
    def __str__(self):
        return 'Spatial(%s, %s, %s, %d)'%(str(self.itd_max), str(self.angular),
                                          str(self.fractional_itd),
                                          self.fractional_delay_filterlength)
    __repr__ = __str__


class DiscreteSpatial(Spatial):
    def __init__(self, itd_max, N, angular=False,
                 fractional_delay_filterlength=2048
                 ):
        self.N = N
        Spatial.__init__(self, itd_max, angular=angular,
                         fractional_itd=True,
                         fractional_delay_filterlength=fractional_delay_filterlength,
                         )
        
    def get_random_itd(self):
        if self.angular:
            theta = randint(self.N)/float(self.N-1)*pi-pi/2
            itd = sin(azim)*self.itd_max
        else:
            itd = randint(self.N)/float(self.N-1)*2*self.itd_max-self.itd_max
        return itd
        
    def __str__(self):
        return 'DiscreteSpatial(%s, %d, %s, %s, %d)'%(str(self.itd_max),
                                                      self.N,
                                                      str(self.angular),
                                                      str(self.fractional_itd),
                                                      self.fractional_delay_filterlength)
    __repr__ = __str__

class HRTFSpatial(Spatial):
    def __init__(self, hrtfset, itd_max):
        self.hrtfset = hrtfset
        self.itd_max = itd_max
        self.bitd_cache_directory = os.path.normpath(os.path.join(datapath, 'optimal_localisation_decoding/cached_bitds'))
        if not os.path.exists(self.bitd_cache_directory):
            os.mkdir(self.bitd_cache_directory)
    def subspace(self, subsetfunc):
        newhrtfset = self.hrtfset.subset(subsetfunc)
        return HRTFSpatial(newhrtfset, self.itd_max)
    def get_random_itd(self):
        index = randint(len(self.hrtfset))
        azim = self.hrtfset.coordinates[index]['azim']
        itd = self.itd_max*sin(azim*pi/180)
        return itd
    def apply(self, sound, itd):
        # pad with zeros to get HRTF applied to the whole sound
        sound = sound[:sound.nsamples+self.hrtfset.num_samples]
        azim = arcsin(itd/self.itd_max)*180/pi
        if azim<0:
            azim += 360
        index = argmin(abs(self.hrtfset.coordinates['azim']-azim))
        hrtf = self.hrtfset[index]
        return hrtf.apply(sound)
    def get_bitds(self, fixedmodel, itdmax=None, phasemultiply=False,
                  continuous=False, interpolation=False, duration_min=50*ms):
        if continuous and phasemultiply:
            raise ValueError('Cannot use continuous and phasemultiply')
        if phasemultiply and itdmax is None:
            raise ValueError('Must specify itdmax if phasemultiply used')
        centre_frequencies = fixedmodel.centre_frequencies
        gammatone_params = fixedmodel.gammatone_params
        cf = centre_frequencies[:len(centre_frequencies)/2]
        # parameters are cf, hrtfset.data, so we hash these and cache to
        # a file
        hrtfset = self.hrtfset
        hash1 = hashlib.sha1(hrtfset.data.flatten().view(uint8)).hexdigest()
        hash2 = hashlib.sha1(hrtfset.coordinates.flatten().view(uint8)).hexdigest()
        hash3 = hashlib.sha1(cf.view(uint8)).hexdigest()
        hash4 = hash((hash1, hash2, hash3))
        fname = os.path.join(self.bitd_cache_directory,
                             str(hash4)+'_'+str((itdmax,
                                                 phasemultiply,
                                                 continuous,
                                                 interpolation,
                                                 duration_min))+'.pkl')
        if os.path.exists(fname):
            return pickle.load(open(fname, 'rb'))
        print 'Computing ITDs (make take a few minutes)'
        all_itds = []
        all_cfind = []
        for i in xrange(len(hrtfset)):
            ir = self.hrtfset[i].impulse_response
            gfb = Gammatone(Repeat(ir, len(cf)), centre_frequencies,
                            **gammatone_params)
            duration = max(ir.nsamples*2/gfb.samplerate, duration_min)
            #print ir.nsamples*2/gfb.samplerate, duration_min
            filtered_hrirset = gfb.process(duration=duration)
            itds = []
            for i in xrange(len(cf)):
                left = filtered_hrirset[:, i]
                right = filtered_hrirset[:, i+len(cf)]
                # This FFT stuff does a correlate(left, right, 'full')
                Lf = fft(hstack((left, zeros(len(left)))))
                Rf = fft(hstack((right[::-1], zeros(len(right)))))
                C = ifft(Lf*Rf).real
                itd = None
                if interpolation:
                    j = argmax(C)
                    x = (arange(j-1, j+2)+1-left.shape[0])/gfb.samplerate
                    y = C[j-1:j+2]
                    p = polyfit(x, y, 2)
                    itd = -0.5*p[1]/p[0]
                if itd is None:
                    itd = (argmax(C)+1-left.shape[0])/gfb.samplerate
                itds.append(itd)
            itds = array(itds)
            if continuous:
                a = itds*2*pi*cf
                a[0] = log(exp(1j*a[0])).imag
                a = hstack((0, cumsum(imag(log(exp(1j*diff(a)))))))+a[0]
                itds = a/(2*pi*cf)
            if phasemultiply:
                cfout, itds, I = phase_multiply(cf, itds, itdmax)
            else:
                I = arange(len(cf))
            all_cfind.append(I)
            all_itds.append(itds)
        result = (all_itds, all_cfind)
        pickle.dump(result, open(fname, 'wb'), -1)
        return result


def apply_adaptive_windowing(hrtfset, aw, reverse, ramptime, cutting=True):
    startmin = hrtfset.data.shape[2]
    startmax = 0
    for i in xrange(2):
        for j in xrange(hrtfset.num_indices):
            D = hrtfset.data[i, j, :]
            onset = min((abs(D)>amax(abs(D))*.15).nonzero()[0])-reverse
            if onset<startmin:
                startmin = onset
            if onset>startmax:
                startmax = onset
            D[:onset] = 0
            D[onset+aw:] = 0
            #D[onset:onset+ramptime] *= hann(2*ramptime)[:ramptime]
            D[onset+aw-ramptime:onset+aw] *= hann(2*ramptime)[ramptime:]
            hrtfset.data[i, j, :] = D
    if cutting:
        hrtfset.data = hrtfset.data[:, :, startmin:startmax+aw].copy()
    hrtfset.hrtf = []
    for i in xrange(hrtfset.num_indices):
        l = Sound(hrtfset.data[0, i, :], samplerate=hrtfset.samplerate)
        r = Sound(hrtfset.data[1, i, :], samplerate=hrtfset.samplerate)
        hrtfset.hrtf.append(HRTF(l, r))
    return hrtfset


class IRCAMHRTFSpatial(HRTFSpatial):
    def __init__(self, subject=1070, itd_max=1*ms, azim_offset=0,
                 windowing=None, cut=None, _elev=0, windowing2=None,
                 adaptive_windowing=None, cutting=True, keep_base=False):
        ircam = get_ircam()
        hrtfset = ircam.load_subject(subject)
        hrtfset.coordinates['azim'] = (hrtfset.coordinates['azim']+azim_offset)%360
        self.windowing = windowing
        self.windowing2 = windowing2
        self.cut = cut
        self.adaptive_windowing = adaptive_windowing
        self._elev = _elev
        if self.cut is not None:
            start, end = cut
            #hrtfset.data = hrtfset.data[:, :, start:end]
            hrtfset.data[:, :, :start] = 0
            hrtfset.data[:, :, end:] = 0           
            hrtfset.hrtf = []
            for i in xrange(hrtfset.num_indices):
                l = Sound(hrtfset.data[0, i, :], samplerate=hrtfset.samplerate)
                r = Sound(hrtfset.data[1, i, :], samplerate=hrtfset.samplerate)
                hrtfset.hrtf.append(HRTF(l, r))
        if self.windowing is not None:
            start, end = windowing
            cutoff = end-start
            # shape (2, num_indices, num_samples)
            window = reshape(hann(cutoff), (1, 1, cutoff))
            #hrtfset.data = hrtfset.data[:, :, start:end]*window
            hrtfset.data[:, :, start:end] = hrtfset.data[:, :, start:end]*window
            hrtfset.data[:, :, :start] = 0
            hrtfset.data[:, :, end:] = 0           
            hrtfset.hrtf = []
            for i in xrange(hrtfset.num_indices):
                l = Sound(hrtfset.data[0, i, :], samplerate=hrtfset.samplerate)
                r = Sound(hrtfset.data[1, i, :], samplerate=hrtfset.samplerate)
                hrtfset.hrtf.append(HRTF(l, r))
        if windowing2 is not None:
            a, b, c, d = windowing2
            hrtfset.data[:, :, a:b] *= reshape(hann((b-a)*2)[:(b-a)], (1, 1, (b-a)))
            hrtfset.data[:, :, c:d] *= reshape(hann((d-c)*2)[(d-c):], (1, 1, (d-c)))
            hrtfset.data[:, :, :a] = 0
            hrtfset.data[:, :, d:] = 0           
            hrtfset.hrtf = []
            for i in xrange(hrtfset.num_indices):
                l = Sound(hrtfset.data[0, i, :], samplerate=hrtfset.samplerate)
                r = Sound(hrtfset.data[1, i, :], samplerate=hrtfset.samplerate)
                hrtfset.hrtf.append(HRTF(l, r))
        if adaptive_windowing is not None:
            aw, reverse, ramptime = adaptive_windowing
            hrtfset = apply_adaptive_windowing(hrtfset, aw, reverse, ramptime,
                                               cutting=cutting)
#            startmin = hrtfset.data.shape[2]
#            startmax = 0
#            for i in xrange(2):
#                for j in xrange(hrtfset.num_indices):
#                    D = hrtfset.data[i, j, :]
#                    onset = min((abs(D)>amax(abs(D))*.15).nonzero()[0])-reverse
#                    if onset<startmin:
#                        startmin = onset
#                    if onset>startmax:
#                        startmax = onset
#                    D[:onset] = 0
#                    D[onset+aw:] = 0
#                    #D[onset:onset+ramptime] *= hann(2*ramptime)[:ramptime]
#                    D[onset+aw-ramptime:onset+aw] *= hann(2*ramptime)[ramptime:]
#                    hrtfset.data[i, j, :] = D
#            if cutting:
#                hrtfset.data = hrtfset.data[:, :, startmin:startmax+aw].copy()
#            hrtfset.hrtf = []
#            for i in xrange(hrtfset.num_indices):
#                l = Sound(hrtfset.data[0, i, :], samplerate=hrtfset.samplerate)
#                r = Sound(hrtfset.data[1, i, :], samplerate=hrtfset.samplerate)
#                hrtfset.hrtf.append(HRTF(l, r))
        if keep_base:
            self.base_hrtfset = hrtfset
        hrtfset = hrtfset.subset(lambda azim, elev: ((azim<=90)|(azim>=270))&(elev==0))
        self.subject = subject
        self.azim_offset = azim_offset
        self.cutting = cutting
        HRTFSpatial.__init__(self, hrtfset, itd_max)
    def __str__(self):
        extraargs = ''
        if self.azim_offset!=0:
            extraargs += ','+str(self.azim_offset)
        if self.windowing is not None:
            extraargs += ','+str(self.windowing)
        if self.windowing2 is not None:
            extraargs += ','+str(self.windowing2)
        if self.cut is not None:
            extraargs += ','+str(self.cut)
        if self._elev is not None:
            extraargs += ','+str(self._elev)
        if self.adaptive_windowing is not None:
            extraargs += ','+str(self.adaptive_windowing)
        if self.cutting:
            extraargs += ',cutting=True'
        return 'IRCAMHRTFSpatial(%d,%s%s)'%(self.subject, str(self.itd_max),
                                            extraargs)


class TollinCatHRTFSpatial(HRTFSpatial):
    def __init__(self, itd_max=440*usecond, adaptive_windowing=None,
                 cutting=True):
        db = TollinCatDatabase()
        hrtfset = db.load_subject()
        self.adaptive_windowing = adaptive_windowing
        if adaptive_windowing is not None:
            aw, reverse, ramptime = adaptive_windowing
            hrtfset = apply_adaptive_windowing(hrtfset, aw, reverse, ramptime,
                                               cutting=cutting)
        HRTFSpatial.__init__(self, hrtfset, itd_max)
    def __str__(self):
        return '%s(%s,%s)'%(self.__class__.__name__, str(self.itd_max),
                            str(self.adaptive_windowing))


class WagnerOwlHRTFSpatial(HRTFSpatial):
    def __init__(self, itd_max=350*usecond, adaptive_windowing=None,
                 cutting=True):
        db = WagnerOwlDatabase()
        hrtfset = db.load_subject('Scottie')
        hrtfset = hrtfset.subset(lambda azim, elev: (elev==0)&((azim<=90)|(azim>=270)))
        self.adaptive_windowing = adaptive_windowing
        if adaptive_windowing is not None:
            aw, reverse, ramptime = adaptive_windowing
            hrtfset = apply_adaptive_windowing(hrtfset, aw, reverse, ramptime,
                                               cutting=cutting)
        HRTFSpatial.__init__(self, hrtfset, itd_max)
    def __str__(self):
        return '%s(%s,%s)'%(self.__class__.__name__, str(self.itd_max),
                            str(self.adaptive_windowing))
        

class SphericalHeadSpatial(Spatial):
    def __init__(self, radius, distance, azim_ear, elev_ear):
        self.radius = radius
        self.distance = distance
        self.azim_ear = azim_ear
        self.elev_ear = elev_ear
        self.itd_max = 2*radius/(343.2*metre/second)
    def get_random_itd(self):
        theta = rand()*pi-pi/2
        itd = sin(theta)*self.itd_max
        return itd
    def apply(self, sound, itd):
        azim = float(arcsin(itd/self.itd_max))
        ir = hrir_sphere(self.radius, self.distance, self.azim_ear, self.elev_ear,
                         azim, 4096, 1024)
        hrtf = HRTF(ir)
        return hrtf.apply(sound)
    def __str__(self):
        return 'SphericalHeadSpatial(%s,%s,%s,%s)'%(self.radius,
                                                    self.distance,
                                                    self.azim_ear,
                                                    self.elev_ear)

if __name__=='__main__':
    if 0:
        space = IRCAMHRTFSpatial(subject=3012, #azim_offset=-25,
                                 #windowing=(5570-1024, 5570+1024)
                                 adaptive_windowing=(764, 16, 512),
                                 cutting=True,
                                 )
        for i in xrange(space.hrtfset.num_indices):
        #for i in [15, 31]:
            ir = space.hrtfset[i].impulse_response.left
            plot(ir)
        show()
    if 1:
        #windowlen = (764, 16, 512)
        windowlen = (512, 16, 256)
        #windowlen = None
        title(str(windowlen))
        space = WagnerOwlHRTFSpatial(adaptive_windowing=windowlen)
        #space = TollinCatHRTFSpatial(adaptive_windowing=windowlen)
        #space = IRCAMHRTFSpatial(subject=3012, azim_offset=-25,
        #                         adaptive_windowing=windowlen,
        #                         #windowing=(5570-2048*2, 5570+2048*2)
        #                         )
        cf = erbspace(400*Hz, 1.5*kHz, 50)
        fm = Parameters(centre_frequencies=hstack((cf, cf)),
                        gammatone_params={},
                        )
        phasemultiply = False
        all_bitds, all_cfind = space.get_bitds(fm,
                                               itdmax=1.0*ms,
                                               continuous=False,
                                               interpolation=True,
                                               phasemultiply=phasemultiply)
        azims = space.hrtfset.coordinates['azim']
        for azim, bitd, cfind in zip(azims, all_bitds, all_cfind):
            if azim>=270: azim -= 360
            c = (azim+90)/180.
            c = (c, 1-c, 0)
            #if azim==0:
            if azim==0 or azim==30 or azim==60 or azim==90:
                c = (0, 0, 1)
            if phasemultiply:
                plot(cf[cfind], -bitd/usecond, '.', color=c)
            else:
                plot(cf, -bitd/usecond)#, color=c, label=str(azim))
#            for n in [-1, 1]:
#                plot(cf, (bitd+n/cf)/usecond, '.k')
        axhline(0, ls='--', color='k')
        show()
        
    if 0:
        for i, spatial in enumerate([
                        IRCAMHRTFSpatial(subject=1002),
                        IRCAMHRTFSpatial(subject=1070),
                        IRCAMHRTFSpatial(subject=2034),
                        IRCAMHRTFSpatial(subject=3010),
                        ]):
            for theta in [90]:#[0, 30, 60, 90]:
                ir = spatial.hrtfset(azim=theta, elev=0).impulse_response
                ir = Sound(ascontiguousarray(ir), samplerate=ir.samplerate)
                print spatial.subject, ir.shape
                left = ir[:, 0].flatten()
                right = ir[:, 1].flatten()
                freq = fftfreq(ir.shape[0], float(1/ir.samplerate))
                I = (freq>50)&(freq<2000)
                fL = fft(left)[I]
                fR = fft(right)[I]
                freq = freq[I]
                subplot(221)
                plot(freq, smooth(10*log10(abs(fL)/abs(fR))),
                     label=str(spatial.subject))
                subplot(222)
                plot(freq, log(fL/fR).imag, '-')                
                subplot(224)
                def docorrect(itd, cf):
                    itd = array(itd)
                    a = itd*2*pi*cf
                    a[0] = log(exp(1j*a[0])).imag
                    a = hstack((0, cumsum(imag(log(exp(1j*diff(a)))))))+a[0]
                    return a/(2*pi*cf)
                itd = log(fL/fR).imag/(2*pi*freq)
                itd = docorrect(itd, freq)
                plot(freq, itd/usecond, '-')
                subplot(223)
                ir = Sound(Cascade(ir, LowPass(ir, 100*Hz), 8).process(),
                           samplerate=ir.samplerate)
                plot(ir.times/ms, ir[:, 0]/amax(abs(ir[:, 0])))
        subplot(221)
        legend(loc='lower right')
        ylabel('ILD (dB)')
        subplot(222)
        ylabel('IPD (rad)')
        subplot(223)
        title('Impulse response (<100 Hz)')
        xlabel('Time (ms)')
        yticks([])
        subplot(224)
        ylabel('ITD (us)')
        xlabel('Frequency (Hz)')            
        show()
    if 0:
        from spatializer.dsp.binauralcues import *
        # IRCAM numbers may not be meaningful because the samplerates are
        # 192 kHz and itd_continuous appears to assume they are 44.1*kHz
        #spatial = IRCAMHRTFSpatial(subject=3012) # guinea pig
        #spatial = IRCAMHRTFSpatial(subject=3014) # cat
        #spatial = IRCAMHRTFSpatial(subject=3016) # macaque
        # Sampling rate problem with the 192 kHz files, the two things below
        # give very different results and th 44.1 kHz ones look right.
        #spatial = IRCAMHRTFSpatial(subject=2034) # rat @ 44.1 kHz
        #spatial = IRCAMHRTFSpatial(subject=3010) # rat @ 192 kHz
        for i, spatial in enumerate([
                        IRCAMHRTFSpatial(subject=2034),
                        #IRCAMHRTFSpatial(subject=3010),
                        #IRCAMHRTFSpatial(subject=3010),
                        ]):
            for theta in [90]:#[0, 30, 60, 90]:
                #ir = hrir_sphere(5*cmetre/2, 2*metre, 90*pi/180, theta*pi/180, 0*pi/180, 4096, 1024)
                ir = spatial.hrtfset(azim=theta, elev=0).impulse_response
                ir = Sound(ascontiguousarray(ir), samplerate=ir.samplerate)
                #ir = Sound(Butterworth(ir, 2, 4, 400*kHz).process(),
                #           samplerate=ir.samplerate)
#                ir = Sound(IIRFilterbank(ir, 2,
#                                         #passband, stopband,
#                                         10*kHz, 20*kHz,
#                                         #gpass, gstop,
#                                         5*dB, 40*dB,
#                                         #btype, ftype
#                                         'high', 'butter'
#                                         ).process(duration=ir.duration+100*ms),
#                           samplerate=ir.samplerate)
                # n=8 appears to be sufficient to get rid of aliasing problems
                ir = Sound(Cascade(ir, LowPass(ir, 20*kHz), 8).process(),
                           samplerate=ir.samplerate)
                #ir.samplerate = 192*kHz
                #if i==1:
                #    ir = ir[::4, :]
                #    ir.samplerate /= 4
                #if ir.samplerate>50*kHz:
                #    ir = ir[::4, :]
                #    ir.samplerate /= 4
                #    print 'woop'
                    #ir[ir.shape[0]/2:, :] = 0
                #print ir.samplerate
                #ir = 3*clip(ir, 0, Inf)**(1./3)
                ir = ImpulseResponse(ir, coordinates=(theta, 0), binaural=True,
                                     samplerate=ir.samplerate)
                #print ir.samplerate
                cf = linspace(50*Hz, 2*kHz, 100)
                itd = itd_continuous(ir, cf, window=1)
                #itd = itd*(44.1*kHz/(192*kHz)) # correct for IRCAM animals
                plot(cf, itd/usecond)
                #semilogx(cf, itd)
        show()
    if 0:
        from spatializer.dsp.binauralcues import *
#        spatial = IRCAMHRTFSpatial(subject=2034) # rat @ 44.1 kHz
#        #spatial = IRCAMHRTFSpatial(subject=3010) # rat @ 192 kHz
#        ir = spatial.hrtfset(azim=90, elev=0).impulse_response
#        #ir = spatial.hrtfset(azim=90, elev=0).apply(whitenoise(200*ms, samplerate=ir.samplerate))
#        #ir = Sound(ascontiguousarray(ir), samplerate=ir.samplerate)
#        #ir = Sound(Butterworth(ir, 2, 4, 20*kHz).process(),
#        #           samplerate=ir.samplerate)
#        #ir = ir[::2, :]
#        #ir.samplerate /= 2
        azim = 90*pi/180
        azim_ear = 90*pi/180
        elev_ear = 0*pi/180
        ir = hrir_sphere(28*mmetre, 2*metre, azim_ear, elev_ear, azim, 4096, 1024)
        cf = linspace(50*Hz, 2000*Hz, 100)
        C = gammatone_correlate(ir, ir.samplerate, cf)
        allex = []
        for i in xrange(len(cf)):
            xmax, ymax = local_maximas(C[:, i], window=int(samplerate/cf[i]), n=3)
            ex, ey, q = interpolated_extremum(xmax, ymax)
            allex.append(ex)
        allex = (array(allex)+1-ir.shape[0])/samplerate
        C = C/reshape(amax(C, axis=0), (1, C.shape[1]))
        imshow(C.T, origin='lower left', aspect='auto', interpolation='nearest',
               extent=(float((1-ir.shape[0])/ir.samplerate/usecond),
                       float((2+ir.shape[0])/ir.samplerate/usecond),
                       0,100))
        plot((argmax(C, axis=0)+1-ir.shape[0])/ir.samplerate/usecond, arange(100), '+w')
        #plot(allex/usecond, arange(100), 'xw')
        show()
    if 0:
        from spatializer.dsp.binauralcues import *
        spatial = IRCAMHRTFSpatial(subject=3010) # rat @ 192 kHz
        ir = spatial.hrtfset(azim=90, elev=0).impulse_response
        ir = Sound(ascontiguousarray(ir), samplerate=ir.samplerate)
        ir = Sound(Butterworth(ir, 2, 12, 20*kHz).process(),
                   samplerate=ir.samplerate)
        ir = 3*clip(ir, 0, Inf)**(1./3)
        ir = ImpulseResponse(ir, coordinates=(90, 0), binaural=True,
                             samplerate=ir.samplerate)
        sound = ir.listen()
        sound.left.spectrogram()
        show()
    if 0:
        from spatializer.dsp.binauralcues import *
#        spatial = IRCAMHRTFSpatial(subject=3010) # rat @ 192 kHz
#        ir = spatial.hrtfset(azim=90, elev=0).impulse_response
        radius = 20*mmetre
        title('radius = %s'%str(radius))
        azim = -90*pi/180
        azim_ear = 110*pi/180
        elev_ear = 30*pi/180
        ir = hrir_sphere(radius, 2*metre, azim_ear, elev_ear, azim, 4096, 1024)
        cf = linspace(50*Hz, 5000*Hz, 200)
        
        ir = ImpulseResponse(ir, coordinates=(90, 0), binaural=True,
                             samplerate=ir.samplerate)
        itd = itd_continuous(ir, cf, window=10)
        plot(cf, itd/usecond)
        show()
    if 0:
        #S = IRCAMHRTFSpatial()
        S = SphericalHeadSpatial(20*cmetre, 2*metre, 100*pi/180)
        itd = S.get_random_itd()
        print itd
        snd = whitenoise(10*ms)
        LR = S.apply(snd, itd)
        L, R = LR.left, LR.right
        f = 5*kHz
        L = Sound(Gammatone(L, f).process())
        R = Sound(Gammatone(R, f).process())
        plot(L.times, L/amax(L))
        plot(R.times, R/amax(R))
        show()
    if 0:
        ir = hrir_sphere(20*cm, 2*metre, 100*pi/180, 30*pi/180, 0)
        sound = whitenoise(100*ms)
        hrtf = HRTF(ir)
        newsound2 = hrtf.apply(sound)
        subplot(311)
        plot(sound.times, sound)
        subplot(313)
        plot(ir.times, ir)
        subplot(312)
        plot(newsound2.times, newsound2)
        show()
    if 0:
        S = Spatial(120*usecond)
        itd = S.get_random_itd()
        print itd
        snd = whitenoise(10*ms)
        LR = S.apply(snd, itd)
        L, R = LR.left, LR.right
        f = 2*kHz
        L = Sound(Gammatone(L, f).process())
        R = Sound(Gammatone(R, f).process())
        plot(L.times, L)
        plot(R.times, R)
        if itd>0:
            plot(L.times+itd, L)
        else:
            plot(R.times-itd, R)
        show()
