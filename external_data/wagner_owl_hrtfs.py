from base import *
import numpy as np

__all__ = ['wagnerHRIR',
           'WagnerOwlDatabase',
           ]

owls = ['Zhaki', 'Pavarotti', 'Hurinfine', 'Hurin', 'Goldi', 'Scottie']

wagner_hrtf_path = search_for_path(wagner_owl_locations)

def wagnerHRIR(name='Scottie'):
    '''
    This function loads the deconvolved HRTFs of the barn owls measured by
    Herman Wagner.

    Returns the irs, the coordinates array and the samplerate

    Data is arranged so that hrir for coordinates coordinates[i] is
    (left) data[:,i] and (right) data[:,i+data.shape[1]/2]
    '''
    f = open(os.path.join(wagner_hrtf_path, 'owl_'+name+'_meta.pkl'), 'rb')
    meta = pickle.load(f)
    f.close()
        
    f = open(os.path.join(wagner_hrtf_path, 'owl_'+name+'_data.npy'), 'rb')
    data = np.load(f)
    f.close()
    
    return data, meta['coordinates'], meta['samplerate']


class WagnerOwlDatabase(HRTFDatabase):
    def __init__(self, path=wagner_hrtf_path):
        self.path = path
        self.name = 'WagnerOwl'
    def load_subject(self, subject='Scottie'):
        rawdata, coords, samplerate = wagnerHRIR(subject)
        nsamples = rawdata.shape[0]
        ncoords = len(coords)/2
        data = zeros((2, ncoords, nsamples))
        for i in xrange(ncoords):
            data[0, i, :] = rawdata[:, i]
            data[1, i, :] = rawdata[:, i+ncoords]
        #azim = array(hstack((arange(270, 360, 10), arange(0, 270, 10))))
        azim = coords['azim'][:ncoords]
        azim[azim<0] += 360
        coords = make_coordinates(azim=azim, elev=coords['elev'][:ncoords])
        return HRTFSet(data, samplerate, coords)


if __name__=='__main__':
    data, coords, samplerate = wagnerHRIR()
    coords = coords[:len(coords)/2]
    I, = ((coords['azim']==0.0)&(coords['elev']==0.0)).nonzero()
    print I
#    print amax(abs(data[:, I[0]]-data[:, I[1]]))
    plot(data[:, I[0]], '-r')
#    plot(data[:, I[1]], '-g')
    plot(data[:, I[0]+len(coords)], '--r')
#    plot(data[:, I[1]+len(coords)], '--g')
#    plot(coords['azim']+randn(len(coords)),
#         coords['elev']+randn(len(coords)),
#         '.',
#         )
    show()
    exit()
    print data.shape
    print coords.dtype
    print len(coords)
    print samplerate
    db = WagnerOwlDatabase()
    hrtfset = db.load_subject()
    hrtfset = hrtfset.subset(lambda azim, elev: (elev==0) & ((azim<=90) | (azim>=270)))
    print len(hrtfset)
    print hrtfset.coordinates['azim']
    for i, azim in enumerate([270, 0, 90]):
        subplot(3, 1, i+1)
        print azim
        plot(hrtfset(azim=azim, elev=0).left)
        plot(hrtfset(azim=azim, elev=0).right)
    show()
    