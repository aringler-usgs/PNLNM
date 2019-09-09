#!/usr/bin/env python
from obspy.signal import PPSD
import glob
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import numpy as np

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

pernlnm, nlnm = get_nlnm()
pernhnm, nhnm = get_nhnm()

N= 5
path = '/home/aringler/inst_test/sts6AGU/horizontalfigure/'
fig = plt.figure(1,figsize=(14,9))
plt.subplots_adjust(hspace=0.001)
models = {}
for idx, sen in enumerate(['I','S']):
    ax1 = fig.add_subplot(3,1,idx+1)
    plt.semilogx(pernlnm,nlnm,color='k',linewidth=2.)
    plt.semilogx(pernhnm,nhnm,color='k',linewidth=2.)
    chans = ['LHH', 'LHZ']
    for chan in chans:
        if chan == 'LHH':
            files = glob.glob(path + sen + '*LH1.npz') + glob.glob(path + sen + '*LH2.npz')
        else:
            files = glob.glob(path + sen + '*LHZ.npz')
        for curfile in files:
            ppsd = PPSD.load_npz(curfile)
            per, percent = ppsd.get_percentile(10.)
            
            if np.mean(percent) <= -150.:
                print(np.mean(percent))
                continue
            percent = np.convolve(percent, np.ones((N,))/N, mode='same')
            

            
            if 'A' not in vars():
                A = percent
            else:
                A = np.vstack([A, percent])
            if chan == 'LHZ':
                color='C1'
                
            else:
                color='C0'
            plt.semilogx(per,percent,alpha=.1, color=color)
            plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
            if idx == 1:
                plt.ylabel('Power (dB  rel. 1 $(m/s^2)^2/Hz)$)')

        
        #plt.xlabel('Period (s)')
        plt.xlim((2.5,1200.))
        plt.ylim((-195.,-90.))
        plt.yticks([-190, -150, -110])
        if idx == 0:
            plt.text(3, -193, 'KS-54000')
        else:
            plt.text(3, -193, 'STS-6')
        models[sen + chan] = np.percentile(A,10,axis=0)

ax1 = fig.add_subplot(3,1,1)
plt.semilogx(per, models['ILHH'],color='C0', linewidth=2)
plt.semilogx(per, models['ILHZ'],color='C1', linewidth=2)
plt.text(1, -100, '(a)', fontsize=28)
ax1 = fig.add_subplot(3,1,2)
plt.semilogx(per, models['SLHH'],color='C0', linewidth=2)
plt.semilogx(per, models['SLHZ'],color='C1', linewidth=2)
plt.text(1, -100, '(b)', fontsize=28)

ax1 = fig.add_subplot(3,1,3)
plt.text(1, 6, '(c)', fontsize=28)
ax1.semilogx(per, models['ILHH']-models['SLHH'],color='C0',label='LHH')
ax1.semilogx(per, models['ILHZ']-models['SLHZ'],color='C1', label='LHZ')
ax1.semilogx(pernlnm,nlnm,color='k',linewidth=2., label='NHNM/NLNM')
plt.text(3,-6.82,'Improvement of STS-6')
plt.xlabel('Period (s)')
plt.xlim((2.5,1200.))
plt.ylim((-7.,7.))

ax1 = plt.gca()
hand, lab = ax1.get_legend_handles_labels()

handles = hand
    
labels = lab
print(handles)
print(labels)
leg = fig.legend(handles, labels, loc = 'lower center', ncol = 3, fontsize = 15)
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)
plt.subplots_adjust(bottom = 0.14)
#plt.show()
plt.savefig('figures/figures13.jpg',format='JPEG',dpi=400)
