#!/usr/bin/env python
from obspy.core import UTCDateTime, read, Stream
from obspy.signal.invsim import evalresp
from matplotlib import mlab, rcParams
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.util import smooth
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)
debug = True

sta = 'ANMO'
net = 'IU'
loc = '00'
nfft = 24*60*60
perM= 3000.
perm = 3.
stime = UTCDateTime('2018-195T00:00:00.0')
#stime += 10.*24*60.*60.
#etime = stime+ 10.*60.*60.*24.
etime = UTCDateTime('2019-060T00:00:00.0')
specfinal = {}
chans = ['LH1', 'LH2', 'LHZ']
st = Stream()
ctime = stime
while ctime < etime:
    for chan in chans:
        st += read('/msd/' + net + '_' + sta +'/2018/' + str(ctime.julday).zfill(3) + '/00_' + chan +'*')
    ctime += 24.*60.*60.
    st.merge(fill_value=0)
    st.trim(starttime=stime, endtime=etime)

#st.detrend('constant')
#for tr in st:
#    print('Here we are')
#    tr.data= np.correlate(tr.data, tr.data, "same")

print('done')
specs = {}
for tr in st:
    if debug:
        print(tr)
    #olp = int(0.25 * nfft)   
    olp=20.*60.*60.
    specgram, freq, time = mlab.specgram(tr.data, Fs = tr.stats.sampling_rate, NFFT = nfft,
                                        pad_to = nfft * 2, noverlap = olp, scale_by_freq = True, mode='psd')	
    #still  need to remove response
    specgram = specgram[1:,:]
    freq = freq[1:]

    # Grab the response
    resppath = '/APPS/metadata/RESPS/'

    # Return the response
    resp = evalresp(t_samp = tr.stats.delta, nfft = 2 * nfft, filename = resppath + 'RESP.' + net + '.' + \
                                sta + '.' + loc + '.' + chan,  date = tr.stats.starttime, station = sta,
                                channel = chan, network = net, locid = loc, units = 'ACC') 

    # Remove the 0 frequency
    resp = resp[1:]
    # Correct for the response
    specgram = (specgram.T/(np.abs(resp*np.conjugate(resp)))).T
    
    specgram = 10. * np.log10(specgram)
    #for idx in range(len(specgram[1,:])):
    #    specgram[:,idx] = smooth(specgram[:,idx],10)
    specs[str(tr.stats.channel)] = specgram

print(specs)
# We have calculated all the PSDs for the time period now we want to do things to them

def removebad(specgram, time, debug=True):
    vmean = np.mean(specgram[:,:],axis=1)[(freq >= 1./perM) & (freq <= 1./perm)]
    vstd = np.std(specgram[:,:], axis=1)[(freq >= 1./perM) & (freq <= 1./perm)]
    chop1 = np.mean(vmean) +20.
    chop2 = np.mean(vmean) -20.
    if debug:
        print(chop1)
        print(np.mean(vstd))
    badidx =[]
    for idx in range(len(specgram[1,:])):
        curmean = np.mean(specgram[:,idx][(freq >= 1./perM) & (freq <= 1./perm)])
        if (curmean >= chop1) or (curmean <= chop2):
            if debug:
                print(curmean)
                print('Bad index:' + str(idx))
            #Time to remove
            badidx.append(idx)
    # Now that we have the badidx we should remove them 
    newspecgram = np.delete(specgram,badidx,1)
    newtime = np.delete(time, badidx)
    return newspecgram, newtime

# Here we have earthquake removed horizontal hum
times = {}
specsN={}
timesN={}
for chan in chans:
    specsN[chan],timesN[chan] = removebad(specs[chan],time)



#modetypes = ['S','T','R']
modetypes =[]
modes = {}
for mode in modetypes:
    with open('modes_' + mode + '.eigen','r') as f:
        next(f)
        modes[mode] = []
        for line in f:
            line = ' '.join(line.split())
            if int(line.split(' ')[0]) == 0:
                line = line.split(' ')[4]
                modes[mode].append(1./float(line))

#print(modes['T'])


print('Here are the number of spectra:' + str(len(specsN[chan][1,:])))

if True:
    pers, nlnm = get_nlnm()
    fig = plt.figure(1,figsize=(14,9))
    plt.subplots_adjust(hspace=0.001)
    for idx, chan in enumerate(chans):
        ax1 = fig.add_subplot(len(chans),1,idx+1)
        for pidx in range(len(specsN[chan][1,:])):
            ax1.semilogx(freq, specsN[chan][:,pidx], color = 'C0', alpha=.01)

        #if chan == 'LHZ':
        percent=10.
        minsp = np.percentile(specsN[chan][:,:],percent, axis=1)
        #minsp = np.average(specsN[chan][:,:], axis=1)
        
        #for idx2, mode in enumerate(modes['S']):
            #mode *= 1000.
            #if idx2 == 0:
                #ax1.plot([mode, mode],[-296.,-170.],'C2',":", label=r'${}_{0}S_{l}$', alpha=.5)
            #else:
                #ax1.plot([mode, mode],[-296.,-170.],'C2',":", alpha=.5)
        #for idx3, mode in enumerate(modes['T']):
            #mode *= 1000.
            #if idx3 == 0:
                #ax1.plot([mode, mode],[-296.,-170.],'C3', label=r'${}_{0}T_{l}$', alpha=.5) 
            #else:
                #ax1.plot([mode, mode],[-296.,-170.],'C3', alpha=.5) 
        N=10
        minsp = np.convolve(minsp, np.ones((N,))/N, mode='same')
        ax1.semilogx(freq, minsp, color='C1', label=str(int(percent)) + 'th Percentile')
        plt.xlim(((1./perM),(1./perm)))
        plt.plot((1./pers), nlnm,'k',label='NLNM')
        plt.ylim((-201.,-171.))
        plt.yticks([-195., -190.,-185.,-180., -175.,],[-195, -190, -185, -180, -175])
        if (idx + 1) < len(chans):
            plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
        if (idx+1) == len(chans):
            plt.xlabel('Frequency (Hz)')
            
        if idx == 1:
            plt.ylabel('Power (dB  rel. 1 $(m/s^2)^2/Hz)$)')
        
        plt.text(1./2500.,-195., chan)
    
    ax1 = plt.gca()
    hand, lab = ax1.get_legend_handles_labels()
    #del hand[2]
    #del lab[2]
    handles = hand
    
    labels = lab
    print(handles)
    print(labels)
    leg = fig.legend(handles, labels, loc = 'lower center', ncol = 5, fontsize = 15)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.subplots_adjust(bottom = 0.14)
    
    plt.savefig('figures/figure6.jpg',format='JPEG', dpi =400)
    plt.savefig('figures/figure6.pdf',format='pdf', dpi =400)
    plt.show()



