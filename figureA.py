#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
# In displacement
HGLPf = np.asarray([1./40., 1./100.])
HGLPa = np.asarray([((2.*np.pi/40.))*10.**-8, ((2.*np.pi/100.))*5.*10**-8])
HGLPdy = 46.


wwssnf = np.asarray([ 1./40., 1./100.])
wwssna = np.asarray([((2.*np.pi/40.))*9*10.**-7, ((2.*np.pi/100.))*2.*10**-6])
wwssndy = 46.

srof =np.asarray([1./40., 1./100.])
sroa = np.asarray([ ((2.*np.pi/40.))*1.5*10**-9, ((2.*np.pi/100.))*1.5*10**-8])
srody = 120.

NLNMf = np.asarray([1./40., 1./100.])
NLNMa = np.asarray([-187.5, -185.])

fig = plt.figure(1)
for idx in range(2):
    if idx ==  0:
        marker = 'o'
    else:
        marker = 'v'
    plt.plot(20.*np.log10(HGLPa[idx]), HGLPdy, marker, color='C0')
    plt.plot( 20.*np.log10(wwssna[idx]), wwssndy, marker, color='C1')
    plt.plot( 20.*np.log10(sroa[idx]), srody, marker, color='C2')
    plt.plot(NLNMa[idx], 120., marker, color='C3')
plt.ylabel('Dynamic Range (dB)')
plt.xlabel('Lowest Noise')
plt.show()
