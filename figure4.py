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

fig = plt.figure(1,figsize=(12,12))
ax = plt.subplot(111)
# Plot 5th percentile from IRIS
f = open('data/percentile5IU','r')
freqs, amps =[], []
for line in f:
    
    line=line.replace('\r\n','')
    print(line)
    line = line.split(',')
    
    freqs.append(float(line[0]))
    amps.append(float(line[1]))
    
f.close()

freqs =  np.asarray(freqs)
amps = np.asarray(amps)

idx = np.argsort(freqs)
freqs = freqs[idx]
amps = amps[idx]

ax.semilogx(1./freqs, amps, color='C0', label='5th Percentile IU BHZ', linewidth=2)

# Plot NLNM and NHNM
per, nlnm = get_nlnm()
per, nhnm = get_nhnm()

ax.semilogx(per, nhnm, color='k', linewidth=3)
ax.semilogx(per, nlnm, color='k', label='NLNM/NHNM', linewidth=3)









# Plot the sts1 model

f = open('data/lnm_sts1', 'r')

freq = []
power = []
for line in f:
    line = ' '.join(line.split())
    line = line.split(' ')
    print(line)
    freq.append(float(line[1]))
    power.append(float(line[0]))

freq = np.asarray(freq)
power = np.asarray(power)

ax.semilogx(1./freq, power, label='STS-1 Vertical Self-Noise', color='C4', linewidth=2)

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

# Add GSN horizontal noise model 

f = open('data/GSNHmodel2.csv','r')
freqs, amps =[], []
for line in f:
    line = ' '.join(line.split())
    line = line.split(',')
    freqs.append(1./float(line[0]))
    amps.append(float(line[1]))
    
freqs = np.asarray(freqs)
ax.semilogx(1./freqs, amps, label='GSN Horizontal Noise Model', color='C6', linewidth=2)


# Add EG model
#egper = [0.022, 0.035, 0.1, 0.17, 0.4, 0.8, 1.24, 2.4, 4.3, 5., 6., 10., 12., 15.6, 21.9, 31.6, 45, 70, 101, 154, 328, 600, 2000, 10000, 125800]
 
#egLNM = [-163.502, -156.303, -165.413, -169.355, -170.338, -170.525, -167.92, -163.157, -158.826, -156.486, -156.985, -172.776, -175.865, -171.89, -187.047, -192.317, -194.541, -230.774, -226.233, -226.62, -226.996, -233.76, -222.514, -187.802, -157.964]
 
#egHNM = [-106.66, -91.1902, -125.713, -121.844, -122.503, -114.256, -113.551, -120.538, -120.834, -118.592, -124.617, -124.951, -127.319, -127.997, -132.902, -138.665, -145.424, -144.762, -144.49, -142.993, -131.296, -130.05, -123.389, -96.372, -41.3127]

#plt.semilogx(egper, egLNM, color='C5', linewidth=2, label='ENSN NLNM/NHNM', linestyle=':')
#plt.semilogx(egper, egHNM, color='C5', linewidth=2, linestyle=':')


box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])


 


#plt.tight_layout()
plt.ylim((-250., -90.))
plt.xlim((1., 10.**4))
ax.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5, -0.25))
#plt.subplots_adjust(top = 0.97)
plt.xlabel('Period (s)')
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz$)')
#plt.show()
#plt.savefig('figures/figure4.pdf',format='PDF', dpi=400)
plt.savefig('figures/figure4.jpg',format='JPEG', dpi=400)
