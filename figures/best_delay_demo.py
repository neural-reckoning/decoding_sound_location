import os
from base import *
from analyse_model import *
from cells import *
from binauraldisplay import *
from colours import *

from models.mcalpine_guinea_pig import *

physrange = 300
markersize = 10
plot_kwds = {'ms':markersize}

seed(32409822)
curbd = generate_random_mcalpine_et_al_2001_bds(cf)

figure(figsize=(15,5))

def common_formatting():
    axis('tight')
    xlabel('Frequency (Hz)')
    ylabel('Best delay (us)')
    
subplot(121)
I = abs(curbd/usecond)>physrange
plot(cf[-I], curbd[-I]/usecond, '.k', mec='k', **plot_kwds)
plot(cf[I], curbd[I]/usecond, '.r', mec='r', **plot_kwds)
title('Peak decoder\n%d%% of cells are useless' % int(sum(I)*100.0/len(I)))
fill([amin(cf), amax(cf), amax(cf), amin(cf)],
     [-physrange, -physrange, physrange, physrange], color=(0.7,)*3)
common_formatting()

subplot(122)
bp = pi/4/(2*pi*cf)
bestphase = 2*pi*cf*bd
#I = abs(curbd/usecond)>physrange
I = -((abs(bestphase)>pi/8) * (abs(bestphase)<3*pi/8)) 
plot(cf[-I], curbd[-I]/usecond, '.k', mec='k', **plot_kwds)
plot(cf[I], curbd[I]/usecond, '.r', mec='r', **plot_kwds)
for f in [1, -1]:#, 0.5, -0.5, 1.5, -1.5]:
    ls = '-' if int(f)==f else '--'
    plot(cf, f*bp/usecond, ls, lw=1+2.0*(int(f)==f), c='b')
fill_between(cf, bp*0.5/usecond, bp*1.5/usecond, color=(0.7,)*3)
fill_between(cf, -bp*0.5/usecond, -bp*1.5/usecond, color=(0.7,)*3)
title('Hemispheric decoder\n%d%% of cells degrade the code' % int(sum(I)*100.0/len(I)))
common_formatting()

show()
