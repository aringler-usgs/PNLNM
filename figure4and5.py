#!/usr/bin/env python
from obspy.core import read
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy import signal
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)


st = read('/msd/IU_QSPA/2019/001/00_LHZ.512.seed')
st += read('/msd/IU_QSPA/2019/001/40_LFZ.512.seed')
st.detrend('constant')


nfft = 2**12
f, Cxy = signal.coherence(st[0].data, st[1].data, fs = 1., nperseg = nfft, noverlap = int(0.5*nfft))


f = f[1:]
Cxy = Cxy[1:]


st.filter('bandpass', freqmin=0.001, freqmax=0.01)
st.normalize()
t = np.arange(st[0].stats.npts)/(60.*60.)

fig =plt.figure(1, figsize=(12,12))
plt.subplot(2,1,1)
plt.plot(t, st[0].data, label=(st[0].id).replace('.',' '), alpha=.7)
plt.ylabel('Amplitude (Normalized Cnts.)')

plt.plot(t, st[1].data , label=(st[1].id).replace('.',' '), alpha=.7)
plt.ylabel('Amplitude (Normalized Cnts.)')
plt.xlabel('Time (Hours)')
plt.legend()
plt.xlim((min(t), max(t)))
plt.text(-3.2,1.,'(a)')

plt.subplot(2,1,2)
plt.semilogx(1./f, Cxy, label='00 LHZ and 40 LFZ')
plt.xlim((min(1./f),max(1./f)))
plt.ylim((0.,1.))
plt.legend()
plt.xlabel('Period (s)')
plt.ylabel('Coherence ($\gamma^2$)')
plt.text(0.7, 1., '(b)')









plt.savefig('figures/figure4.jpg', format='JPEG',dpi= 400)
plt.clf()
#plt.subplot(3,1,2)



##########################################################################
st = read('/msd/IU_CCM/2019/10*/10_LH1.512.seed')
st += read('/msd/IU_CCM/2019/10*/10_LH2.512.seed')
st += read('/msd/IU_CCM/2019/10*/30_LDO.512.seed')
st.merge()
st.detrend('constant')
st.sort(['location'], reverse=True)
st.filter('bandpass', freqmin=0.001, freqmax=0.02)
st.normalize()
t = np.arange(st[0].stats.npts)/(60.*60.)

fig =plt.figure(1, figsize=(12,12))
for idx, tr in enumerate(st):
    plt.subplot(3,1,idx+1)
    plt.plot(t, tr.data, label=(tr.id).replace('.',' '), alpha=.7)
    if idx == 1:
        plt.ylabel('Amplitude (Normalized Cnts.)')
        plt.text(-37., 1., '(b)', fontsize=28)
    elif idx == 2:
        plt.text(-37., 1., '(c)', fontsize=28)
    elif idx == 0:
        plt.text(-37., 1., '(a)', fontsize=28)
    plt.legend(loc=1)
    plt.xlim((min(t), max(t)))

plt.xlabel('Time (Hours)')

#plt.show()
plt.savefig('figures/figure5.jpg', format='JPEG',dpi= 400)
plt.clf()

##########################################################################
