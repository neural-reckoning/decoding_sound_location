from base import *
from matplotlib import cm

tollin_path = search_for_path(dan_tollin_cat_hrtf_locations)

class TollinCatDatabase(HRTFDatabase):
    def __init__(self, path=tollin_path):
        self.path = path
    def load_subject(self, subject=None):
        path = os.path.join(self.path, 'Raw HRTFs')
        data = zeros((2, 36, 1999))
        for i, azimindex in enumerate(xrange(0, 360, 10)):
            fname = os.path.join(path, 'nov24_'+str(azimindex)+'_7.txt')
            hrtfdata = loadtxt(fname, usecols=(1, 2)) # shape (nsamples, 2)
            hrtfdata -= reshape(mean(hrtfdata, axis=0), (1, 2))
            data[:, i, :] = hrtfdata.T
        samplerate = 97656.25*Hz
        azim = array(hstack((arange(270, 360, 10), arange(0, 270, 10))))
        coords = make_coordinates(azim=azim)
        return HRTFSet(data, samplerate, coords)
            
if __name__=='__main__':
    #from spatializer.dsp.binauralcues import *
    from scipy.signal import hann
    db = TollinCatDatabase(adaptive_windowing=(512, 16, 256))
    hrtfset = db.load_subject()
    print hrtfset.data.shape
    #exit()
    subplot(211)
    for cutoff_factor in [1]:#[1, 2, 3, 4]:
        for theta in arange(0, 360, 10):
            #subplot(2, 2, cutoff_factor)
            #figure()
            #suptitle(str(theta))
            ir = hrtfset(azim=theta).impulse_response
            ir = Sound(ascontiguousarray(ir), samplerate=ir.samplerate)
            if cutoff_factor>0:
                title('duration = %s'%str(ir.duration/cutoff_factor))
                cutoff = ir.shape[0]/cutoff_factor
                ir[cutoff:, :] = 0
                ir[:cutoff, :] *= reshape(hann(2*cutoff)[cutoff:], (cutoff, 1))    
    
            left = ir[:, 0].flatten()
            right = ir[:, 1].flatten()
            freq = fftfreq(ir.shape[0], float(1/ir.samplerate))
            I = (freq>300)&(freq<1500)
            fL = fft(left)[I]
            fR = fft(right)[I]
            freq = freq[I]
            def docorrect(itd, cf):
                itd = array(itd)
                a = itd*2*pi*cf
                a[0] = log(exp(1j*a[0])).imag
                a = hstack((0, cumsum(imag(log(exp(1j*diff(a)))))))+a[0]
                return a/(2*pi*cf)
            itd = log(fL/fR).imag/(2*pi*freq)
            itd = docorrect(itd, freq)
            if theta==90:
                subplot(212)
            if theta==270:
                subplot(211)
            plot(freq, itd/usecond, '-o', label=str(theta),
                 color=cm.hsv(theta/360.))
#            for n in [-1, 1]:
#                plot(freq, (itd+n/freq)/usecond, '--k')
    
    show()
