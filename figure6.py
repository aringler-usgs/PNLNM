#!/usr/bin/env python
from obspy.core import read, UTCDateTime, Stream
from scipy.signal import coherence
import sys
from scipy.optimize import fmin
import numpy as np
from matplotlib.mlab import csd
from math import pi
import sys
import math
from obspy.signal.invsim import evalresp
import matplotlib.pyplot as plt
from obspy.signal.rotate import rotate2zne
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing

import matplotlib as mpl
#Set font parameters using matplotlib
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)


def octavesmooth(x, octaveFraction = 1.0/3.0 ):
    y = []
    rightFactor = math.pow(2.0, octaveFraction / 2.0 )
    leftFactor =  1.0 / rightFactor
    for n in range(len(x)):
        print('On ' + str(n) + ' of ' + str(len(x)))
        left = long(n * leftFactor)
        right = long(n * rightFactor)
        y.append( sum( x[left : right + 1] ) / (right - left + 1) )
		
    return y

def rotme(st, ang):
    azis = ang[0:3]
    dips = ang[3:6]
    data = rotate2zne(st[0].data, azis[0], dips[0], st[1].data, azis[1], dips[1], st[2].data, azis[2], dips[2])
    st[0].data, st[1].data, st[2].data = data[0], data[1], data[2]
    return st

def rotData2(angle, stR, st, nol, nwin, debug = False):
    if debug:
        print('Here is number overlap: ' + str(nol))
        print('Here is the number in the window: ' + str(nwin))
    minper = 2.
    maxper = 10. 
    stTemp = st.copy()    
    stTemp = rotme(stTemp, angle)
    if debug:
        print('Here we are')
        print(stTemp)
        print(stR)
    newcoh = 0.
    if debug:
        print('Here is noverlap: ' + str(int(nol*nwin)))
        print('Here is FFT: ' + str(nwin))
        print('Here is nperseg: ' + str(nwin))
    for idx in range(3):
        (fre, cxy) =  coherence(stR[idx].data,stTemp[idx].data, \
            fs = 1./st[idx].stats.delta, noverlap = int(nol*nwin), nfft = nwin, nperseg = nwin)
        fre = fre[1:]
        cxy = cxy[1:] 
        mask = (fre <= 1./minper) & (fre >= 1./maxper)
        newcoh += np.abs(np.mean(cxy[mask])-1.)
    if debug:
        print(newcoh)

    return newcoh  

def cp(tr1, tr2, lenfft, lenol, delta):
    # Cross-power function
	sr = 1./float(delta)
	cpval,fre = csd(tr1.data, tr2.data, NFFT=lenfft, Fs=sr, noverlap=int(lenol*lenfft), scale_by_freq=True)
	fre = fre[1:]
	cpval = cpval[1:]
	return cpval, fre

def selfnoise(st, length, overlap):
    (p11, f) = cp(st[0],st[0],length,overlap,st[0].stats.delta)
    (p22, f) = cp(st[1],st[1],length,overlap,st[0].stats.delta)
    (p33, f) = cp(st[2],st[2],length,overlap,st[0].stats.delta)

    (p21, f) = cp(st[1],st[0],length,overlap,st[0].stats.delta)
    (p13, f) = cp(st[0],st[2],length,overlap,st[0].stats.delta)
    (p23, f) = cp(st[1],st[2],length,overlap,st[0].stats.delta)


    n={}
    n['1'] = (p11 - p21*p13/p23)
    n['2'] = (p22 - np.conjugate(p23)*p21/np.conjugate(p13))
    n['3'] = (p33 - p23*np.conjugate(p13)/p21)

    p={}
    p['1'] = p11
    p['2'] = p22
    p['3'] = p33
    
    return n, p, f


debug = True
st = Stream()
for day in range(8, 15):
    st += read('/tr1/telemetry_days/II_BFO/2019/2019_' + str(day).zfill(3) + '/*LH*')
#st.trim(endtime = UTCDateTime('2019-009T18:59:59.0'))
for tr in st:
    if tr.stats.channel == 'LH1':
        tr.stats.channel = 'LHN'
    if tr.stats.channel == 'LH2':
        tr.stats.channel = 'LHE'
st.merge()
st.sort(reverse=True)
# Data is in ZNE format
if debug:
    print(st)

comp = 'Z'
length = 400000
overlap = 0.5
#Treat data1 as the reference and rotate 2 and three to match 1
#angVec = [0., 0., 1.]


# For ENZ data
azis = [0., 0., 90.]
dips = [-90., 0., 0.]
angVec = azis + dips


if debug:
    print('Calculating coherence')
#ang1 = fmin(rotData2, angVec, args=(st.select( location='50'), st.select(location='55'), overlap, length))
#ang2 = fmin(rotData2, angVec, args=(st.select(location='50'), st.select(location='60'), overlap, length))

ang1 = [ -8.09887253e-01, 1.04545136e+00, 9.13021324e+01, -9.01370145e+01, -1.74007786e-01, 5.09743818e-02]
ang2 = [ -4.60261201, 3.96006846, 94.10540385, -89.54792117, -0.56071495, -0.39495279]

if debug:
    #print('Here is cohereVal: ' + str(cohereVal1))
    #print('Here is cohereVal: ' + str(cohereVal2))
    print('Here is angle: ' + str(ang1))
    print('Here is angle: ' + str(ang2))

stgood = st.select(component=comp).copy()
for tr in stgood:
    if tr.stats.location == '00':
        stgood.remove(tr)
stgood55= st.select(location='55')
stgood60 = st.select(location='60')
stgood55= rotme(st.select(location='55'), ang1)
stgood60 = rotme(st.select(location='60'), ang2)
stgood = stgood55 + stgood60 + st.select(location='50')
stgood = stgood.select(component=comp)
print(stgood)
#stgood.select(component=comp, location='50')[0].data = st.select(location='50')[0].data*ang1[0] + st.select(location='50')[1].data*ang1[1] + st.select(location='50')[2].data*ang1[2]
#stgood.select(component=comp, location='55')[0].data = st.select(location='55')[0].data*ang2[0] + st.select(location='55')[1].data*ang2[1] + st.select(location='55')[2].data*ang2[2]

if debug:
    print(stgood)

n, p, fre1 = selfnoise(stgood, length, overlap)

fig = plt.figure(1,figsize=(18,18))
plt.subplot(2,1,1)

nm = (n['1'] + n['2'] + n['3'])/3.

for idx in range(1,4):

    tr = stgood[idx-1]
    resp = evalresp(t_samp = tr.stats.delta, nfft = length, filename = '/APPS/metadata/RESPS/RESP.' + tr.id,  
                date = tr.stats.starttime, station = tr.stats.station,
                channel = tr.stats.channel, network = tr.stats.network, 
                locid = tr.stats.location, units = 'ACC')
    
    n[str(idx)] /= np.abs(resp[1:])**2
    p[str(idx)] = 10.*np.log10(p[str(idx)]/np.abs(resp[1:])**2)
    p[str(idx)] = octavesmooth(p[str(idx)], 1./6.)
    n[str(idx)] = 10.*np.log10(n[str(idx)])
    n[str(idx)] = octavesmooth(n[str(idx)], 1./6.)
    plt.semilogx(1./fre1, p[str(idx)], label='PSD ' + (tr.id).replace('.', ' '), alpha=.7)
    plt.semilogx(1./fre1, n[str(idx)], label='Self-Noise ' + (tr.id).replace('.',' '), alpha=.7)
nm /= np.abs(resp[1:])**2
nm = np.abs(nm)
#N=  5
#nm = np.convolve(nm, np.ones((N,))/N, mode='same')

#nm2 = konno_ohmachi_smoothing(nm, fre1)
plt.semilogx(1./fre1, 10.*np.log10(nm), label='Self-Noise Mean')
per_nlnm, pow_nlnm = get_nlnm()
plt.semilogx(per_nlnm,pow_nlnm, linewidth=2, color='k')
per_nhnm, pow_nhnm = get_nhnm()
plt.semilogx(per_nhnm,pow_nhnm, linewidth=2, color='k', label='NLNM/NHNM')

#plt.semilogx(1./fre1, 10.*np.log10(nm2), label='Self-Noise Mean')
plt.xlabel('Period (s)')
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz$)')
#plt.title('Start Time: ' + str(stgood[0].stats.starttime.format_seed()))
plt.legend(loc=9, ncol=4)
plt.ylim((-225,-30))
plt.xlim((2., 500000))
plt.text(0.5, -30., '(a)', fontsize= 28)
#plt.savefig('SELFNOISE_TUC.jpg',format='JPEG', dpi= 400)
#plt.show()
plt.subplot(2,1,2)

########################### Now we can do it with TUC
st = Stream()
for day in range(15, 30):
    st += read('/tr1/telemetry_days/IU_TUC/2019/2019_' + str(day).zfill(3) + '/*LHZ*')
#st.trim(endtime = UTCDateTime('2019-009T18:59:59.0'))
#for tr in st:
    #if tr.stats.channel == 'LH1':
        #tr.stats.channel = 'LHN'
    #if tr.stats.channel == 'LH2':
        #tr.stats.channel = 'LHE'
st.merge()
st.sort(reverse=True)

n, p, fre1 = selfnoise(st, length, overlap)

for idx in range(1,4):

    tr = st[idx-1]
    resp = evalresp(t_samp = tr.stats.delta, nfft = length, filename = '/APPS/metadata/RESPS/RESP.' + tr.id,  
                date = tr.stats.starttime, station = tr.stats.station,
                channel = tr.stats.channel, network = tr.stats.network, 
                locid = tr.stats.location, units = 'ACC')
    
    n[str(idx)] /= np.abs(resp[1:])**2
    p[str(idx)] = 10.*np.log10(p[str(idx)]/np.abs(resp[1:])**2)
    p[str(idx)] = octavesmooth(p[str(idx)], 1./6.)
    n[str(idx)] = 10.*np.log10(n[str(idx)])
    n[str(idx)] = octavesmooth(n[str(idx)], 1./6.)
    plt.semilogx(1./fre1, p[str(idx)], label='PSD ' + (tr.id).replace('.', ' '), alpha=.7)
    plt.semilogx(1./fre1, n[str(idx)], label='Self-Noise ' + (tr.id).replace('.',' '), alpha=.7)
#plt.semilogx(1./fre1, 10.*np.log10(nm), label='Self-Noise Mean')
per_nlnm, pow_nlnm = get_nlnm()
plt.semilogx(per_nlnm,pow_nlnm, linewidth=2, color='k')
per_nhnm, pow_nhnm = get_nhnm()
plt.semilogx(per_nhnm,pow_nhnm, linewidth=2, color='k', label='NLNM/NHNM')

#plt.semilogx(1./fre1, 10.*np.log10(nm2), label='Self-Noise Mean')
plt.xlabel('Period (s)')
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz$)')
#plt.title('Start Time: ' + str(stgood[0].stats.starttime.format_seed()))
plt.legend(loc=9, ncol=4)
plt.ylim((-225,-30))
plt.xlim((2., 500000))
plt.text(0.5, -30., '(b)', fontsize= 28)

plt.savefig('figures/figure6.jpg',format='JPEG', dpi=400)
plt.show()
