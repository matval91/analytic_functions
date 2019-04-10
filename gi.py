# Script to plot Gi and Ge
# Gi = Ec/E0 int[0, E0/Ec]dy/(1+y**1.5)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)    
plt.rc('figure', facecolor='white')

x = np.linspace(0.1, 10, 1000) # this is Ec/E0
gi = np.zeros(len(x), dtype=float)

for i, el in enumerate(x):
    x_array = np.linspace(0,1/el, num=100, dtype=float)
    gi[i] = el*np.trapz(1/(1+x_array**1.5), x_array)
xplot=1/x #we want to plot E0/Ec

f = plt.figure()
ax = f.add_subplot(111)
ax2 = ax.twiny()
ax.semilogx(xplot, gi, 'k', lw=2.3, label='Ions')
ax.semilogx(xplot, 1-gi, 'r', lw=2.3, label='Electrons')
ax.grid('on')
ax.set_xlabel(r'$E_0/E_c$')
ax.set_ylabel(r'G')

ind = np.argmin(gi-0.5<0)
plt.axvline(2.41024, ymax=0.5, color='k', linestyle='-', linewidth=1.5) #Gi=Ge

plt.axvline(1.3, color='m', linestyle='--', linewidth=2.3)#TCV
plt.axvline(2.1, color='g', linestyle='--', linewidth=2.3)#ITER
plt.axvline(3.7, color='b', linestyle='--', linewidth=2.3)#SA

ax.xaxis.set_ticks([0.1, 1., 2.41024, 10.])
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

ax2.set_xscale('log')
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks([1.3, 2.1, 3.7])
ax2.set_xticklabels(['TCV', 'ITER', 'JT-60SA'])
plt.setp(ax2.get_xticklabels(), rotation='45', fontsize=16)
[t.set_color(i) for (i,t) in zip(['magenta','green', 'blue'],ax2.xaxis.get_ticklabels())]

plt.setp(ax.get_yticklabels()[0], visible=False) 

ax.legend(loc='center left')
f.tight_layout()
plt.show()
