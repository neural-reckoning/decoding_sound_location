'''
Code for the spherical head model is from:

ftp://ftp.phon.ucl.ac.uk/pub/stuart/StaceyTrainingShifts/spherical_head_model.pdf

From the paper of Duda and Martens 1998.
'''

from shared import *

def Hsphere(a, r, theta, f, c=343.2*metre/second, threshold=1e-20):
    '''
    a: radius of the sphere (m)
    r: distance from the center of the sphere to the source (m)
    theta: angle of incidence (rad)
    f: frequency (Hz)
    c: ambient speed of sound (m/s)
    
    returns:
    
    H: head-related transfer function relative to free field
    '''
    a = float(a)
    r = float(r)
    f = float(f)
    if abs(f)<1e-20:
        return 1.
    if f<0:
        f = -f
        negf = True
    else:
        negf = False
    c = float(c)
    x = cos(theta)
    mu = (2 * pi * f * a) / c
    rho = r / a
    i = 1j
    zr = 1. / (i * mu * rho)
    za = 1. / (i * mu)
    Qr2 = zr
    Qr1 = zr * (1. - zr)
    Qa2 = za
    Qa1 = za * (1. - za)
    P2 = 1.
    P1 = x
    sum = 0.
    term = zr / (za * (za - 1.))
    sum = sum + term
    term = (3 * x * zr * (zr - 1) ) / (za * (2 * za**2 - 2 * za + 1) )
    sum = sum + term
    oldratio = 1.
    newratio = abs(term)/abs(sum)
    m = 2.
    while (oldratio > threshold) or (newratio > threshold):
        Qr = - (2 * m - 1) * zr * Qr1 + Qr2
        Qa = - (2 * m - 1) * za * Qa1 + Qa2
        P = ( (2 * m - 1) * x * P1 - (m - 1) * P2) / m
        term = ( (2 * m + 1) * P * Qr) / ( (m + 1) * za * Qa - Qa1)
        sum = sum + term
        m = m + 1
        Qr2 = Qr1
        Qr1 = Qr
        Qa2 = Qa1
        Qa1 = Qa
        P2 = P1
        P1 = P
        oldratio = newratio
        newratio = abs(term)/abs(sum)
    H = (rho * exp(- i * mu) * sum) / (i * mu)
    if not negf:
        H = conj(H)
    return H

def single_hrir_sphere(radius, distance, azim_ear, elev_ear, azim_sound,
                       N=512, advance=256):
    #theta = arccos(sin(azim)*sin(theta_ear)+cos(azim)*cos(theta_ear)*cos(elev))
    theta = arccos(cos(elev_ear)*cos(azim_ear-azim_sound))
    freq = fftfreq(N, 1/samplerate)
    H = array([Hsphere(radius, distance, theta, f) for f in freq])
    h = ifft(H).real
    h = roll(h, advance)
    hp = hstack((0, h[:-1]))
    hm = hstack((h[1:], 0))
    h = 0.5*(2*h+hp+hm)
    return h

def hrir_sphere(radius, distance, azim_ear, elev_ear, azim_sound, N=512, advance=256):
    ir_left = single_hrir_sphere(radius, distance, azim_ear, elev_ear, azim_sound,
                                 N=N, advance=advance)
    ir_right = single_hrir_sphere(radius, distance, -azim_ear, elev_ear, azim_sound,
                                 N=N, advance=advance)
    return Sound((ir_left, ir_right))

if __name__=='__main__':    
    for theta in linspace(0, 360, 25):
        freq = fftfreq(256, 1/(44.1*kHz))
        times = arange(len(freq))/(44.1*kHz)
        H = [Hsphere(50*cmetre, 2*metre, theta*pi/180, f) for f in freq]
        H = array(H)
        h = ifft(H).real
        h = roll(h, 128)
        hp = hstack((0, h[:-1]))
        hm = hstack((h[1:], 0))
        h = 0.5*(2*h+hp+hm)
        subplot(221)
        semilogx(freq, 20*log10(abs(H)))
        subplot(223)
        semilogx(freq, log(H).imag)
        subplot(122)
        plot(times, h+theta/90, '-k')
    show()
