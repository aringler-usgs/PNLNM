#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import scipy.io as sio
import matplotlib as mpl
#Set font parameters using matplotlib
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

fig = plt.figure(1,figsize=(10,10))
ax = plt.subplot(111)
# Plot 5th percentile from IRIS
f = open('data/percentile5IU','r')
freqs, amps =[], []
for line in f:
    line = line.split(',')
    freqs.append(float(line[0]))
    amps.append(float(line[1]))
    
f.close()

freqs =  np.asarray(freqs)
amps = np.asarray(amps)

idx = np.argsort(freqs)
freqs = freqs[idx]
amps = amps[idx]

amps = amps[(freqs < 8.)]
freqs = freqs[(freqs < 8.)]


ax.semilogx(1./freqs, amps, color='C0', label='5th Percentile IU BHZ', linewidth=2)

# Plot NLNM and NHNM
per, nlnm = get_nlnm()
per, nhnm = get_nhnm()

ax.semilogx(per, nhnm, color='k', linewidth=3)
ax.semilogx(per, nlnm, color='k', label='NLNM/NHNM', linewidth=3)


# plot Nyquist
ax.semilogx([1./100., 1./100.], [-250, 100], label='200 sps Nyquist', color='C1', linewidth=2)


# Plot the ALNM/AHNM
ahnm=[-91.5, -91.5, -97.41, -110.5, -120, -98, -96.5, -101., -105., -91.25]
ahnmp=[0.01, .1, .22, .32, .8, 3.80, 4.6, 6.3, 7.1, 150]
alnm=[-135, -135, -130, -118.25]
alnmp=[0.01, 1, 10, 150]


ax.semilogx(alnmp, alnm, color='C2', linestyle='--', linewidth=2)
ax.semilogx(ahnmp, ahnm, color='C2', label='ALNM/AHNM', linestyle='--', linewidth=2)

# Plot 1st percentile of QSPA
f = open('data/percentile1QSPA','r')
freqs, amps =[], []
for line in f:
    line = line.split(',')
    freqs.append(float(line[0]))
    amps.append(float(line[1]))
    
f.close()

freqs =  np.asarray(freqs)
amps = np.asarray(amps)

idx = np.argsort(freqs)
freqs = freqs[idx]
amps = amps[idx]


freqs = np.asarray(freqs)
amps = np.asarray(amps)
amps = amps[(freqs <30.)]
freqs = freqs[(freqs < 30.)]

ax.semilogx(1./freqs, amps, color='C3', label='1st Percentile QSPA HHZ', linewidth=2)

# Plot the GS13 model

matgs13 = sio.loadmat('data/lnm_gs13Vertical.mat')
print(matgs13['lnm'])
freq = []
power = []
for ele in matgs13['lnm']:
    freq.append(ele[1])
    power.append(ele[0])

freq = np.asarray(freq)
power = np.asarray(power)

print(power)
ax.semilogx(1./freq, power, label='GS-13 Self-Noise', color='C4', linewidth=2)

# Add the GSN noise model

f = open('data/GSN_noisemodel.txt','r')
freqs, amps =[], []
for line in f:
    line = ' '.join(line.split())
    line = line.split(' ')
    freqs.append(1./float(line[0]))
    amps.append(float(line[1]))
    
    

freqs = np.asarray(freqs)
ax.semilogx(1./freqs, amps, label='GSN Vertical Noise Model', color='.5', linewidth=2)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

per = [0.833333333, 0.555555556, 0.25, 0.2, 0.1, 0.05, 0.033333333]


LNM=[-165.87, -168.83, -169.22, -168.35, -168.07, -167.00, -162.97]	

ax.semilogx(per, LNM, label='Lajitas (1984)', color='C5', linewidth=2)


#plt.tight_layout()
plt.ylim((-195., -80.))
plt.xlim((1.5*10**-3, 10.))
ax.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5, -0.25))
#plt.subplots_adjust(top = 0.97)
plt.xlabel('Period (s)')
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz$)')
#plt.show()
plt.savefig('figures/figure2.pdf',format='PDF', dpi=400)
plt.savefig('figures/figure2.jpg',format='JPEG', dpi=400)
