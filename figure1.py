#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

networks = ['GSN', 'International' ,  'FDSN' , 'Earthscope', 'Portable' , 'US Regional' , 'Engineering' ,'Other']
vals = np.asarray([125./25., 110./25.,  130./55., 70./42., 195./135., 155./150., 5./13., 1./14.])



ships = [125., 110.,   130., 70., 195., 155., 5., 1.]
data= [25., 25., 14., 55., 42., 135., 150., 13.]
inds = range(len(vals))

fig = plt.figure(1, figsize=(12,12))
plt.subplot(2,1,1)
ind = np.arange(len(vals))/float(len(vals))
cols = cm.viridis(np.flip(ind))



plt.scatter(inds,vals*100., s=vals*300., c=cols)
for pair in zip(networks, inds, vals*100.):
    plt.text(pair[1] +.0, pair[2], pair[0])
plt.xlim((-1.,len(inds)+1))
plt.xticks([])
plt.ylim((0, 600.))
plt.title('2018 Shipments per byte at IRIS DMC')
plt.ylabel('Turnover (\%)')
plt.text(-2., 650., '(a)', fontsize=26)

plt.subplot(2,2,3)

plt.pie(ships, colors=cols, labels=networks, autopct="%1.1f%%")
plt.title('Shipped Data (TB)')
plt.text(-1.5, 1.4, '(b)', fontsize=26)
plt.subplot(2,2,4)
plt.pie(data, colors=cols, labels=networks, autopct="%1.1f%%", startangle=20.)
plt.title('Data Holdings (TB)')
plt.text(-1.5, 1.4, '(c)', fontsize=26)
plt.savefig('figures/figure1.jpg', format='JPEG', dpi=400)

#plt.show()
