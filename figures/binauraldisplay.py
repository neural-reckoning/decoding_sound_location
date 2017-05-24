from base import *
from scipy.spatial import cKDTree

def best_delay_distribution(cf, bd, usephase=False, tightaxis=True,
                            markersize=4):
    if usephase:
        y = 2*pi*cf*bd
        ylabel('Best phase (rad)')
    else:
        y = bd/msecond
        ylabel('Best delay (ms)')
    plot(cf/kHz, y, '.k', ms=markersize)
    xlabel('Best frequency (kHz)')
    if tightaxis: axis('tight')

def display_response(cf, bd, response, itd=None, use_bestphase=False,
                     tightaxis=True,
                     min_marker_size=5, max_marker_size=20,
                     mincol=0.2, maxcol=0.8,
                     show_xlabel=True,
                     show_ylabel=True,
                     coloured=True,
                     formatting=None,
                     bginterp='cells', # None, 'cells', 'griddata'
                     dmax=0.04, npixels=600, # for cells background
                     bgcol=(1.0,)*3,
                     cmap=cm.jet_r,
                     colour=None,
                     horiz_stretch=1.0,
                     ):
    if formatting is None:
        formatting = {
            'edgecolors':'none',
            'cmap':cmap,
            }
        if colour is None:
            colour = (mincol, mincol, mincol)
        if not coloured:
            formatting['color'] = colour
    force_ylims = None
    I = argsort(response)
    if isinstance(show_ylabel, str):
        prepend = show_ylabel+'\n'
    else:
        prepend = ''
    if use_bestphase:
        x = cf[I]/kHz
        y = 2*pi*cf[I]*bd[I]
        if show_ylabel:
            ylabel(prepend+'Best phase (rad)')
    else:
        x = cf[I]/kHz
        y = bd[I]/msecond
        if show_ylabel:
            ylabel(prepend+'Best delay (ms)')
    nresponse = response[I]*1.0/amax(response)
    if coloured:
        formatting['c'] = (1-nresponse)*(maxcol-mincol)+mincol
    if bginterp is not None:
        if bginterp=='griddata':
            z = (1-nresponse)*(maxcol-mincol)+mincol
            xi = linspace(amin(x), amax(x), 400)
            yi = linspace(amin(y), amax(y), 400)
            # use delaunay triangulation to interpolate
            zi = griddata(x, y, z, xi, yi)
            # use radial basis functions
    #        from scipy.interpolate import Rbf
    #        rbf = Rbf(x, y, z, epsilon=0.1)
    #        Xi, Yi = meshgrid(xi, yi)
    #        zi = rbf(Xi, Yi)
            imshow(zi, origin='lower left',
                   extent=(amin(x), amax(x), amin(y), amax(y)),
                   cmap=formatting['cmap'])
        elif bginterp=='cells':
            xscale = amax(x)-amin(x)
            yscale = amax(y)-amin(y)
            xscale = xscale*horiz_stretch
            k = cKDTree(array([x/xscale, y/yscale]).T)
            u, v = meshgrid(linspace(amin(x), amax(x), npixels),
                            linspace(amin(y), amax(y), npixels))
            u.shape = u.size
            v.shape = v.size
            d, i = k.query(array([u/xscale, v/yscale]).T)
            d.shape = (npixels, npixels)
            i.shape = d.shape
            s = (1-nresponse)*(maxcol-mincol)+mincol
            C = s[i]
            C[d>dmax] = nan
            imshow(C, origin='lower left', aspect='auto',
                   extent=(amin(x), amax(x), amin(y), amax(y)),
                   cmap=formatting['cmap'])
            gca().set_axis_bgcolor(bgcol)
            force_ylims = amin(y), amax(y)
            
    # sort plot order by largest last
    s = nresponse*(max_marker_size-min_marker_size)+min_marker_size
    I = argsort(s)
    x = x[I]
    y = y[I]
    s = s[I]
    scatter(x, y, s=s, vmin=0, vmax=1, marker='o', **formatting)
    if itd is not None:
        if use_bestphase:
            plot(cf/kHz, 2*pi*cf*itd, '--k', lw=2)
        else:
            axhline(itd/msecond, ls='--', color='k', lw=2)                
    if tightaxis: axis('tight')
    if force_ylims is not None:
        ylim(*force_ylims)
    elif not use_bestphase:
        ylim(ymin=-1, ymax=1)
    xlim(xmin=amin(x), xmax=amax(x))
    if show_xlabel:
        if isinstance(show_xlabel, str):
            post = '\n'+show_xlabel
        else:
            post = ''
        xlabel('Best frequency (kHz)'+post)
