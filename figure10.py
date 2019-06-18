#!/usrbing/env python

from obspy.core import read, UTCDateTime, Stream
import matplotlib.pyplot as plt
from scipy.signal import periodogram, hilbert
from obspy.signal.invsim import evalresp
import numpy as np
from scipy.optimize import fmin
import math
import glob
import sys
from obspy.clients.fdsn import Client 
from scipy import signal
debug = True

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

# network
net ='IU'
# station
sta ='ANMO'
# location can be wild carded
loc='00'

# Event time

stime = UTCDateTime('2018-08-19T00:19:40')

stime2 = stime
# min and max frequency in milliHz
mf= 0.5
Mf= 1.3

# In days
length = 60

diresyn = 'data'

# If we have a pressure sensor we can do a pressure correction using pcorr as true
pcorr = False

# Hours of data you want analyzed for sub 1 mHz you want to get 3 or more days
hours = 70
hours = int((1000.*1.1*195./0.76593)/(60.*60.))
#hours = 70
print('Using ' + str(hours) + ' hours of data')

# List of reference modes to plot (only good for sub 1 mHz)
modes = [[0.8140, '0S0'], [0.3092, '0S2'], [0.46884, '0S3'], [0.6469, '0S4'], [0.67934, '1S2'], 
            [0.84022, '0S5'], [0.93908, '1S3'], [0.94434, '3S1'], [1.0379, '0S6'], [1.23133, '0S7'], [1.4135,'0S8'], 
            [1.5783, '0S9'], [1.7265,'0S10'],[1.1720, '1S4'], [1.37,'1S5'], [1.2412, '2S3'],[1.3797,'2S4'], [1.1055,'3S2'],
            [1.4116, '4S1'],
            [0.37889, '0T2'], [0.5856, '0T3'], [0.7649, '0T4'], [0.9272, '0T5'],[1.0775,'0T6'],[1.219,'0T7'],
            [1.35660, '0T8'], [1.31912,'1T2'], [1.43835, '1T3'],[1.48707, '0T9']]
resppath = '/APPS/metadata/RESPS/RESP.'


st = Stream()

client = Client("IRIS")
inventory = client.get_stations(network=net, station=sta, starttime=stime2, endtime = stime2 + 24.*60.*60., level="response")

days = int(float(hours)/24.) + 1
if days < 1:
    days = 1

for day in range(days):
    ctime = stime2 + 24.*60.*60.*day
    string = '/msd/' + net + '_' + sta + '/' + str(ctime.year) + '/' + str(ctime.julday).zfill(3) + '/' 
    #string = '/tr1/telemetry_days/' + net + '_' + sta + '/' + str(ctime.year) + '/*' + str(ctime.julday).zfill(3) + '/' 
    st += read(string + loc + '_LH*.seed')
    if pcorr:
        st += read(string + '30_LDO*.seed')


st.trim(starttime=stime2, endtime=stime2+hours*60.*60.)
st.merge()
st.detrend('linear')
st.detrend('constant')

print(st)

st.merge()
print(st)

#st.decimate(10)
#st.decimate(10)
st.rotate(method="->ZNE", inventory=inventory)
for tr in st:
    if tr.stats.channel == 'LHN':
        tr.stats.channel = "LH1"
    if tr.stats.channel == "LHE":
        tr.stats.channel = "LH2"

# decimate to 1 sample per 20 s
st.decimate(5)
st.decimate(2)
st.decimate(2)

for tr in st:
    win = signal.get_window(('kaiser', 2.*np.pi), tr.stats.npts)
    tr.data *= win

for tr in st:
    tr.data = np.pad(tr.data, (int((length*60*60*24/20. -tr.stats.npts)/2.),int((length*60*60*24/20. -tr.stats.npts)/2.)), 'edge') 


NFFT=2**(math.ceil(math.log(st[0].stats.npts, 2)))

if pcorr:
    trP = st.select(channel='LDO')[0]
if debug:
    if pcorr:
        print(trP)

# We now have the data tie to compute spectra
fig = plt.figure(1, figsize=(16,12))
for idx, tr in enumerate(st.select(channel='LH*')):
    if debug:
        print(tr)
    f,p= periodogram(tr.data, fs=tr.stats.sampling_rate, nfft= NFFT, scaling='spectrum')
    p = p[1:]
    f= f[1:]
    resppathALL = resppath + str(tr.stats.network) + '.' + str(tr.stats.station) + '.' + str(tr.stats.location) + '.' + str(tr.stats.channel)
    if debug:
        print(resppathALL)
    resp = evalresp(t_samp = tr.stats.delta, nfft= NFFT, filename = resppathALL,
                    date = tr.stats.starttime, station = tr.stats.station, channel= tr.stats.channel,
                    network= tr.stats.network, locid= tr.stats.location, units='ACC')
    resp = resp[1:]
    # Convert units to nm/s/s
    p = np.sqrt(p/(np.abs(resp)**2))*10**9
    if debug:
        print(len(f))
    # Now have p in nm/s/s switch f to mHz
    f *= 1000.
    p = p[(f >= mf) & (f <= Mf)]
    if debug:
        print('Here is length of p' + str(len(p)) )
    f = f[(f >= mf) & (f <= Mf)]
    
    
    # Make a static correction
    def pressure_correct(x):
        a = x[0]
        b = x[1]
        # Correction based on a and b
        tr_correct = tr.data - a*trP.data - b*np.imag(hilbert(trP.data)) 
        
        # Now compute the residual
        f, trp = periodogram(tr_correct,  fs=tr.stats.sampling_rate, nfft= NFFT, scaling='spectrum')
        trp = trp[1:]
        f = f[1:]*1000.
        trp = np.sqrt(trp/(np.abs(resp)**2))*10**9
        # Now we have the spectra so we want to minimize the regions outside the modes
        resi = 0.
        bands = [[0.2,0.29], [0.32,0.36], [0.43, 0.45]]
        for band in bands:
            resi += sum(trp[(f>= band[0]) & (f <= band[1])]**2)
        print(resi)
        return resi
    
    # So we want to minimize the correction to identify a and b
    if pcorr:
        bf = fmin(pressure_correct, [0., 0.], maxiter=20)
    #bf =[0., 0.]
    if debug:
        if pcorr:
            print('Here is bf' + str(bf))
    
    if pcorr:
        tr_correct = tr.data - bf[0]*trP.data - bf[1]*np.imag(hilbert(trP.data))
    
    if debug:
        print(len(f))
        if pcorr:
            print(len(tr_correct))
    if pcorr:
        f, pGood =  periodogram(tr_correct , fs=tr.stats.sampling_rate, nfft=NFFT, scaling='spectrum')
        f= f[1:]*1000.
        pGood = pGood[1:]
    
        pGood = np.sqrt(pGood/(np.abs(resp)**2))*10**9
        pGood = pGood[(f >= mf) & (f <= Mf)]
        print('Here is the length of f' + str(len(f)))
        f = f[(f >= mf) & (f <= Mf)]
    
    if debug:
        print('Here is the length of p' + str(len(p)))
        print('Here is the length of f' + str(len(f)))
    if tr.stats.channel == 'LHZ':
        pidx = 1
        label='LHZ'
    elif tr.stats.channel in ['LH1','LHN']:
        pidx = 2
        label = 'LHN'
    elif tr.stats.channel in ['LH2', 'LHE']:
        pidx = 3
        label='LHE'
    plt.subplot(3,1,pidx)
    plt.plot(f,p, label=tr.stats.location + ' ' + label, alpha=0.7)
    plt.fill_between(f, 0, p, alpha=0.7)
    if pcorr:
        plt.plot(f,pGood, label=tr.stats.location + ' ' + tr.stats.channel+ ' Pressure Corrected')
    plt.legend(loc='upper left', fontsize=12)
    plt.xlim((min(f),max(f)))
    if pidx == 3:
        plt.xlabel('Frequency ($mHz$)')
    elif pidx == 2:
        plt.ylabel('Amplitude ($nm/s^2$)')
    if pidx == 1:
        plt.title(tr.stats.network + ' ' + tr.stats.station + ' Start Time: ' +  str(tr.stats.starttime.year)
                    + ' ' + str(tr.stats.starttime.julday).zfill(3) + ' ' + str(tr.stats.starttime.hour).zfill(2) + ':' + 
                    str(tr.stats.starttime.minute).zfill(2) + ':' + str(tr.stats.starttime.second) + ' Duration: ' +
                    str(hours) + ' Hours')
    
    if pidx == 1:
        plt.text(.4, 0.007, '(a)', fontsize=28)
    if pidx == 2:
        plt.text(.4, 0.007, '(b)', fontsize=28)
    if pidx == 3:
        plt.text(.4, 0.007, '(c)', fontsize=28)
    
    
    
    # Plot modes
    if idx < 3:
        locp = max(p)*1.1
        for nidx, mode in enumerate(modes):
            if (mode[0] <= max(f)) and (mode[0] >= min(f)):
                label = mode[1]
                label = '_' + label[0] + label[1] + '_' + label[2]
                if len(label) == 4:
                    label += label[3]
                if debug:
                    print(label)
                plt.plot((mode[0], mode[0]), (-10., 500.), color='grey', alpha=0.5)
                
                plt.text(mode[0] + .003, locp, '$' + label + '$' , fontsize=14)
             
                locp -= max(p)*0.15
                if locp < 0.75*max(p):
                    locp = 1.1*max(p)
    #plt.ylim((-.0, max(p)*1.1))
    plt.ylim((0., 0.007))
# Time to plot a synthetic

#st = read(diresyn + '/' + sta + "*.modes.sac")
#st.trim(starttime=stime, endtime=stime+hours*60.*60.)

#st.detrend('linear')
#st.detrend('constant')
#for tr in st:
    #win = signal.get_window(('kaiser', 2.*np.pi), tr.stats.npts)
    #tr.data *= win

#for tr in st:
    #tr.data = np.pad(tr.data, (int((length*60*60*24/20 -tr.stats.npts)/2.),int((length*60*60*24/20 -tr.stats.npts)/2.)), 'edge')  
#st.merge()



#if debug:
    #print(st)

#for tr in st:
    #f, p =  periodogram(tr.data, fs=tr.stats.sampling_rate, nfft=NFFT, scaling='spectrum')
    #f= f[1:]*1000.
    #p = np.sqrt(p[1:])
    #if tr.stats.channel == 'LHZ':
        #pidx = 1
    #elif tr.stats.channel in ['LH1','LHN']:
        #pidx = 2
    #elif tr.stats.channel in ['LH2', 'LHE']:
        #pidx = 3
    #plt.subplot(3,1,pidx)
    #plt.plot(f,p, label=tr.stats.location + ' ' + (tr.stats.channel).replace('H','X'), alpha=0.7)
    #plt.legend(loc='upper left', fontsize=12)
##sys.exit()



plt.savefig('figures/figure10.jpg', format='JPEG', dpi=400)  


