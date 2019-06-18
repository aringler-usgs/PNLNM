#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import read, UTCDateTime
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import matplotlib as mpl
from matplotlib.mlab import csd
from obspy.signal.invsim import evalresp
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

z1,z2 = 162., 23.

vp, vs =3500., 2200.

sta = 'DWPF'

f = np.linspace(0.001,10,1000.)

dz1 = 2.*z1
dz2 = 2.*z2
tp1 = dz1/vp
ts1 = dz1/vs
tp2 = dz2/vp
ts2 = dz2/vs

pamp1 = 0.5 + 0.5*np.cos(2.*np.pi*f*tp1)
samp1 = 0.5 + 0.5*np.cos(2.*np.pi*f*ts1)
pamp2 = 0.5 + 0.5*np.cos(2.*np.pi*f*tp2)
samp2 = 0.5 + 0.5*np.cos(2.*np.pi*f*ts2)

print(pamp1)

fig = plt.figure(1, figsize=(14,14))
plt.subplot(3,1,1)
plt.plot(1./f, pamp1, label='P-wave 00 Amplitude')
plt.plot(1./f, pamp2, label='P-wave 10 Amplitude')
plt.text(0.15, 1.00, '(a)', fontsize=28)
plt.semilogx(1./f, samp1, label='S-wave 00 Amplitude')
plt.semilogx(1./f, samp2, label='S-wave 10 Amplitude')
plt.legend()
plt.xlim((0.25, 20))
plt.ylabel('Amplitude (Normalized)')
plt.xlabel('Period (s)')
plt.subplot(3,1,2)
plt.semilogx(1./f, 20.*np.log10(pamp1/pamp2), label='P-wave 00 to 10 Ratio (dB)')
plt.semilogx(1./f, 20.*np.log10(samp1/samp2), label='S-wave 00 to 10 Ratio (dB)')
plt.ylabel('Power (dB)')
plt.text(0.15, 0.2, '(b)', fontsize=28)
plt.legend()
plt.xlim((0.25, 20))
plt.xlabel('Period (s)')
#plt.show()

stime = UTCDateTime('2019-101T23:30:00')
st = read('/msd/IU_' + sta + '/2019/' + str(stime.julday).zfill(3) + '/*_B*Z*')
st.trim(stime, stime+60.*60.*5)
st.sort(reverse=True)

for tr in st:
    if tr.stats.location == 'XX':
        print('In the mud')
        seedresp = {'filename': 'RESP.' + tr.stats.network + '.' + tr.stats.station + '.' + tr.stats.location + '.' + tr.stats.channel, 'date': stime, 'units':'ACC'}
    else:
        seedresp = {'filename': '/APPS/metadata/RESPS/RESP.' + tr.stats.network + '.' + tr.stats.station + '.' + tr.stats.location + '.' + tr.stats.channel, 'date': stime, 'units':'ACC'}
    #tr.detrend('constant')
    
    #tr.simulate(paz_remove=None, seedresp=seedresp)
    #tr.filter('bandpass',freqmin=1., freqmax =5.)
    #tr.taper(0.05)



lenfft= 2*512
#fig = plt.figure(1, figsize=(12,8))
plt.subplot(3,1,3)
per, NLNM = get_nlnm()
per, NHNM = get_nhnm()
for tr in st:
    power,freq = csd(tr.data,tr.data, NFFT = lenfft, 
                noverlap = int(lenfft*.5), Fs = 1./tr.stats.delta, 
                            scale_by_freq = True)
    power = np.abs(power[1:])
    freq = freq[1:]
    resp = evalresp(t_samp = tr.stats.delta, nfft = lenfft, filename= '/APPS/metadata/RESPS/RESP.' + tr.stats.network + '.' + tr.stats.station + '.' + tr.stats.location + '.' + tr.stats.channel,  
                    date = tr.stats.starttime, units = 'ACC') 
    resp = resp[1:]

    power = 10.*np.log10(power/np.absolute(resp*np.conjugate(resp)))   
    plt.semilogx(1./freq,power, label=tr.id)
plt.semilogx(per, NLNM,'k', label="NHNM/NLNM")
plt.semilogx(per, NHNM,'k')
plt.ylabel('Power (dB rel. 1 $m/(s^2)^2/Hz$)')  
plt.xlabel('Period (s)')
plt.text(0.137, -110., '(c)', fontsize=28)
plt.legend(loc=4)
plt.xlim((.25, 20.))
plt.ylim((-161., -109.))

plt.savefig('figures/figure7.jpg', format='JPEG', dpi= 400)
#plt.show()
