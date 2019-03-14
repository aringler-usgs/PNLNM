#!/usr/bin/env python
from obspy.core import read, UTCDateTime
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory
from obspy.signal.invsim import cosine_taper, paz_to_freq_resp
from matplotlib.mlab import csd
from math import pi
import sys
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import scipy.io as sio
import matplotlib as mpl
#Set font parameters using matplotlib
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

fig = plt.figure(1,figsize=(10,10))



debug = True
NFFT= 2*2048

def cp(tr1,tr2,lenfft,lenol,delta):
	sr = 1./delta
	cpval,fre = csd(tr1.data,tr2.data,NFFT=lenfft,Fs=sr,noverlap=lenol,scale_by_freq=True)
	fre = fre[1:]
	cpval = cpval[1:]
	return cpval, fre




st = read('/TEST_ARCHIVE/XX_ENG1/2012/2012_116/*BH*')
st.trim(starttime=UTCDateTime('2012-116T15:13:00.0'), endtime=UTCDateTime('2012-116T15:31:00.0'))

if debug:
	print(st)
fig = plt.figure(1,figsize=(10,10))
per_nlnm, pow_nlnm = get_nlnm()
plt.semilogx(per_nlnm,pow_nlnm, linewidth=2, color='k')
per_nhnm, pow_nhnm = get_nhnm()
plt.semilogx(per_nhnm,pow_nhnm, linewidth=2, color='k', label='NLNM')


inv = read_inventory('/APPS/metadata/RESPS/RESP.IU.TUC.00.BHZ')
sts6paz = inv.get_response('IU.TUC.00.BHZ', UTCDateTime('2019-001T00:00:00.0'))
sts6paz = sts6paz.get_paz()
sts6paz.zeros= sts6paz.zeros[1:]

print(sts6paz)

inv = read_inventory('/APPS/metadata/RESPS/RESP.IU.MAJO.10.BHZ')
sts2paz = inv.get_response('IU.MAJO.10.BHZ', UTCDateTime('2019-001T00:00:00.0'))
sts2paz = sts2paz.get_paz()
sts2paz.zeros = sts2paz.zeros[1:]

print(sts2paz)


for tr in st.select(location='00'):
	f, p2 = signal.welch(tr.data, fs = tr.stats.sampling_rate, nperseg=NFFT, noverlap=1024)
	f = f[1:]
	p2 = p2[1:]
	if 'p' in vars():
		p = np.vstack((p, p2))
   	else:
   		p = p2
	


p = np.mean(p, axis=0)
Tfnom, fnom = paz_to_freq_resp(sts6paz.poles, sts6paz.zeros, sts6paz.stage_gain*sts6paz.normalization_factor, 1./st[0].stats.sampling_rate, NFFT, freq=True)
Tfnom = Tfnom[1:]*(2**26)/40.
plt.semilogx(1./f,10.*np.log10(p/np.abs(Tfnom)**2), label='Q330HR Noise with STS-6 Response')
Tfnom, fnom = paz_to_freq_resp(sts2paz.poles, sts2paz.zeros, sts2paz.stage_gain*sts2paz.normalization_factor, 1./st[0].stats.sampling_rate, NFFT, freq=True)
Tfnom = Tfnom[1:]*(2**26)/40.
plt.semilogx(1./f,10.*np.log10(p/np.abs(Tfnom)**2), label='Q330HR Noise with STS-2HG Response')


#plt.show()


###################### Add STS-2HG self-noise


instresp = np.abs(Tfnom)**2

st = read('/TEST_ARCHIVE/XX_ENG7/2016/2016_180/00_BH0*')
st += read('/TEST_ARCHIVE/XX_ENG4/2016/2016_180/10_BH0*')
st += read('/TEST_ARCHIVE/XX_ENG5/2016/2016_180/00_BH0*')


st.merge()
st.sort(['station'])
print(st)

st[0].data = st[0].data.astype(np.float64)
st[0].data *= 4.

(p11, fre1) = cp(st[0],st[0],NFFT,1024,st[0].stats.delta)
(p22, fre1) = cp(st[1],st[1],NFFT,1024,st[0].stats.delta)
(p33, fre1) = cp(st[2],st[2],NFFT,1024,st[0].stats.delta)

(p21, fre1) = cp(st[1],st[0],NFFT,1024,st[0].stats.delta)
(p13, fre1) = cp(st[0],st[2],NFFT,1024,st[0].stats.delta)
(p23, fre1) = cp(st[1],st[2],NFFT,1024,st[0].stats.delta)



n11 = (p11 - p21*p13/p23)/instresp
n22 = (p22 - np.conjugate(p23)*p21/np.conjugate(p13))/instresp
n33 = (p33 - p23*np.conjugate(p13)/p21)/instresp


nn = (n11 + n22)/2.

#plt.semilogx(fre1, 10.*np.log10(p11/(instresp)), label=st[0].id)
#plt.semilogx(fre1, 10.*np.log10(p22/instresp), label=st[1].id)
#plt.semilogx(fre1, 10.*np.log10(p33/(instresp)), label=st[2].id)
plt.semilogx(1./fre1, 10.*np.log10(nn), label='STS-2HG Self-Noise')
#plt.semilogx(fre1, 10.*np.log10(n11), label='STS-2HG Self-Noise')
#plt.semilogx(fre1, 10.*np.log10(n22), label='STS-2HG Self-Noise')
#plt.semilogx(fre1, 10.*np.log10(n33), label='STS-2HG Self-Noise')

#plt.show()
##################### Add STS-6 self-noise

Tfnom, fnom = paz_to_freq_resp(sts6paz.poles, sts6paz.zeros, sts6paz.stage_gain*sts6paz.normalization_factor, 1./st[0].stats.sampling_rate, NFFT, freq=True)
Tfnom = Tfnom[1:]*(2**26)/40.

instresp = np.abs(Tfnom)**2

#st = read('/tr1/telemetry_days/II_BFO/2019/2019_002/50_BHZ*')
#st += read('/tr1/telemetry_days/II_BFO/2019/2019_002/55_BHZ*')
#st += read('/tr1/telemetry_days/II_BFO/2019/2019_002/60_BHZ*')

st = read('/tr1/telemetry_days/XX_CTAO/2019/2019_070/*BHZ*')
st += read('/tr1/telemetry_days/XX_GSN1/2019/2019_070/00_BHZ*')

# 00 is an STS-6 with a pre-amp
# 10 is an STS-2.5
# 00 GSN1 is an STS-6

Tfnom, fnom = paz_to_freq_resp(sts6paz.poles, sts6paz.zeros, sts6paz.stage_gain*sts6paz.normalization_factor, 1./st[0].stats.sampling_rate, NFFT, freq=True)
Tfnomsts6 = Tfnom[1:]*(2**26)/40.
Tfnomsts6pa = Tfnomsts6*20.

inv = read_inventory('/APPS/metadata/RESPS/RESP.IU.TUC.10.BHZ')
sts25paz = inv.get_response('IU.TUC.10.BHZ', UTCDateTime('2019-001T00:00:00.0'))
sts25paz = sts25paz.get_paz()
sts25paz.zeros= sts25paz.zeros[1:]
Tfnom, fnom = paz_to_freq_resp(sts25paz.poles, sts25paz.zeros, sts25paz.stage_gain*sts25paz.normalization_factor, 1./st[0].stats.sampling_rate, NFFT, freq=True)
Tfnomsts25 = Tfnom[1:]*(2**26)/40.

instresp1 = np.abs(Tfnomsts6pa)**2
instresp2 = np.abs(Tfnomsts25)**2
instresp3 = np.abs(Tfnomsts6)**2

st.merge()
st.sort()
if debug:
	print(st)


(p11, fre1) = cp(st[0],st[0],NFFT,1024,st[0].stats.delta)
(p22, fre1) = cp(st[1],st[1],NFFT,1024,st[0].stats.delta)
(p33, fre1) = cp(st[2],st[2],NFFT,1024,st[0].stats.delta)

(p21, fre1) = cp(st[1],st[0],NFFT,1024,st[0].stats.delta)
(p13, fre1) = cp(st[0],st[2],NFFT,1024,st[0].stats.delta)
(p23, fre1) = cp(st[1],st[2],NFFT,1024,st[0].stats.delta)


n11 = (p11 - p21*p13/p23)/instresp1
n22 = (p22 - np.conjugate(p23)*p21/np.conjugate(p13))/instresp2
n33 = (p33 - p23*np.conjugate(p13)/p21)/instresp3

#nn = (n11 + n22 + n33)/3. 

per = [0.833333333, 0.555555556, 0.25, 0.2, 0.1, 0.05, 0.033333333]


LNM=[-165.87, -168.83, -169.22, -168.35, -168.07, -167.00, -162.97]	

plt.semilogx(per, LNM, label='Lajitas (1984)', color='C5', linewidth=2)

plt.semilogx(1./fre1, 10.*np.log10(n11), label='STS-6 Self-Noise')
#plt.semilogx(1./fre1, 10.*np.log10(n22), label='STS-2.5 Self-Noise')
#plt.semilogx(1./fre1, 10.*np.log10(n33), label='STS-6H Self-Noise')
#plt.semilogx(1./fre1, 10.*np.log10(p11/instresp1), label='STS-6 Noise')
#plt.semilogx(1./fre1, 10.*np.log10(p22/instresp2), label='STS-2.5 Noise')
#plt.semilogx(1./fre1, 10.*np.log10(p33/instresp3), label='STS-6H Self-Noise')
plt.xlim((1./20.,1.))
plt.ylim((-210, -140.))
plt.xlabel('Period (s)')
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz$)')
plt.legend()
plt.savefig('figures/figure1.jpg', format='JPEG', dpi= 400)
plt.savefig('figures/figure1.pdf', format='PDF', dpi= 400)
#plt.show()
