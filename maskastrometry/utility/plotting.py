
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.

"""Provide some plotting functions useful for distortion solutions."""


import numpy as np
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator

from maskastrometry.distortions.polynomials import PolyVectorField, Legendre
from maskastrometry.distortions.polyfit import polyfit_svd


def set_fontsize(small=8, medium=10, bigger=12):
    """Set font size of ticks, axes labels and titles."""
    plt.rc('font', size=medium)          # controls default text sizes
    plt.rc('axes', titlesize=medium)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=medium)    # legend fontsize
    plt.rc('figure', titlesize=bigger)   # fontsize of the figure title


def plot_distortions(ax, positions, distortions,
                     inliers=None, limits=(-1.1, 1.1),
                     positionunits=(1, "normalized"), positionunits2=None,
                     arrowunits="arcsec", keylen=None,
                     keypos=(0.7, 1.05), dotsize=2):
    """Make arrow plot of distortions.

    Parameters
    ----------
    ax : matplotlib.Axes
        The target axes used for plotting.
    positions : (N, 2)-shaped array
        The grid of positions in normalised coordinates.
    distortions : (N, 2)-shaped array
        The distortions.
    inliers : (N,)-shaped bool array or None
        Optional array that masks outliers.  `True` for inliers.
    limits : 2-tuple
        Lower and upper limit for plot axes in normalised coordinates.
        Same in x and y direction.
    positionunits : 2-tuple of (float, str)
        Scale of normalised units to another unit and the other unit's name.
        Use e.g. the plate scale and `"arcsec"`.
    positionunits2 : 2-tuple of (float, str)
        Like `positionunits` for a second scale,
        currently unsupported.
    arrowunits : str
        The unit name for distortion data, e.g. `"uas"`.
    keylen : float
        Length of key for arrows.
        If set to `None` the RMS value is used.
    keypos : 2-tuple of floats
        The position of the key in axes coordinates.
    dotsize : float
        Diameter of outlier dots in pt.
    """
    ps = positionunits[0]
    if inliers is None:
        inliers = np.full((positions.shape[0],), True)
    quiver = ax.quiver(positions[inliers, 0] * ps, positions[inliers, 1] * ps,
                       distortions[inliers, 0], distortions[inliers, 1],
                       pivot='middle')
    ax.scatter(positions[~inliers, 0] * ps, positions[~inliers, 1] * ps,
               s=dotsize**2, marker='.', color='red')

    if keylen is None:
        rmslen = np.sqrt(np.sum(distortions[inliers]**2)
                         / distortions[inliers].size)
        keylen = rmslen * 2
        keylen = np.round(keylen, -int(np.floor(np.log10(keylen)))+1)
    ax.quiverkey(quiver, *keypos, keylen,
                 label="{:.2g} {:s}".format(keylen, arrowunits),
                 labelpos='E', coordinates='axes')

    ax.set_xlabel("position / {:s}".format(positionunits[1]))
    ax.set_ylabel("position / {:s}".format(positionunits[1]))
    if positionunits2:
        raise NotImplementedError

    if limits:
        limits = np.array(limits) * ps
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
    ax.set_aspect('equal')


def plot_residuals(fig, positions, residuals, inliers=None, limits=(-1.1, 1.1),
                   positionunits=(1, "normalized"),
                   arrowunits="arcsec", keylen=None,
                   keypos=(1.15, 1.2), dotsize=1, **kwgridspec):
    """Plot residual map with marginalized distributions.

    You pass a figure and keyword arguments to matplotlib.GridSpec to define
    where all five subplots should go.

    Parameters
    ----------
    fig : matplotlib.Figure
        Target figure for plotting.
    positions : (N, 2)-shaped array
        See `plot_distortions()` for doc.
    residuals : (N, 2)-shaped array
    inliers :  (N,)-shaped bool array
    limits : 2-tuple of floats
    positionunits : 2-tuple of (floats, str)
    arrowunits : str
    keylen : float or None
    keypos : 2-tuple of floats
    dotsize : float
    kwgridspec : keyword arguments
        Use keywords `left`, `right`, `bottom`, `top`, `wspace`, `hspace`
        to define the position of the subplots.

    Examples
    --------

    >>> plot_residuals(
    ...     fig, positions, residuals * 1e6,
    ...     positionunits=(posscale, "arcsec"),
    ...     arrowunits="uas",
    ...     left=0.52, right=0.9, bottom=0.52, top=0.9,
    ...     wspace=0.05, hspace=0.05)

    """
    ps = positionunits[0]
    if inliers is None:
        inliers = np.full((positions.shape[0],), True)

    gs = GridSpec(3, 3, width_ratios=[1, 3, 1], height_ratios=[1, 3, 1])
    gs.update(**kwgridspec)

    ax = fig.add_subplot(gs[1, 1])
    plot_distortions(ax, positions, residuals, inliers, limits,
                     positionunits, None, arrowunits, keylen,
                     keypos, dotsize)
    ax.tick_params(axis='both', which='both', right='on', top='on',
                   labelleft='off', labelbottom='off')
    ax.set_xlabel('')
    ax.set_ylabel('')

    axmargxdx = fig.add_subplot(gs[1, 0], sharey=ax)
    axmargxdx.tick_params(axis='both', which='both', left='on', right='off',
                          labelleft='on', labelright='off')
    axmargxdx.scatter(residuals[inliers, 0], positions[inliers, 1] * ps,
                      s=dotsize**2, c='k', alpha=0.2)
    axmargxdx.set_ylabel("position / {:s}".format(positionunits[1]))
    axmargxdx.set_xlabel("res/{:s}".format(arrowunits))

    axmargxdy = fig.add_subplot(gs[1, 2], sharey=ax)
    axmargxdy.tick_params(axis='both', which='both', left='off', right='on',
                          labelleft='off', labelright='on',
                          top='on', bottom='off',
                          labeltop='on', labelbottom='off')
    axmargxdy.scatter(residuals[inliers, 1], positions[inliers, 1] * ps,
                      s=dotsize**2, c='k', alpha=0.2)

    axmargydx = fig.add_subplot(gs[0, 1], sharex=ax)
    axmargydx.tick_params(axis='both', which='both', left='on', right='off',
                          labelleft='on', labelright='off',
                          top='on', bottom='off',
                          labeltop='on', labelbottom='off')
    axmargydx.scatter(positions[inliers, 0] * ps, residuals[inliers, 0],
                      s=dotsize**2, c='k', alpha=0.2)

    axmargydy = fig.add_subplot(gs[2, 1], sharex=ax)
    axmargydy.tick_params(axis='both', which='both', left='off', right='on',
                          labelleft='off', labelright='on',
                          top='off', bottom='on',
                          labeltop='off', labelbottom='on')
    axmargydy.scatter(positions[inliers, 0] * ps, residuals[inliers, 1],
                      s=dotsize**2, c='k', alpha=0.2)
    axmargydy.set_xlabel("position / {:s}".format(positionunits[1]))

    axmargxdx.text(0.3, 1.18, "x",
                   horizontalalignment='center',
                   transform=axmargxdx.transAxes)
    axmargxdy.text(0.8, -0.2, "y",
                   horizontalalignment='center',
                   transform=axmargxdy.transAxes)


def make_grid_analysis(fig, positions, distortions, posscale, maxorder,
                       name="distortions"):
    """Run analysis of distortion pattern and display in four panels.

    The four panels are:

    1. Arrow plot of distortion pattern itself.
    2. Arrow plot with marginalized distributions of residuals of
       Legendre-fit at `maxorder`.
    3. RMS residuals over maximum order of polynomials.
       Maximum order from 1 to `maxorder`.
    4. RMS power in each degree of polynomial of
       fit at `maxorder`.

    Some information is printed to the console.

    Parameters
    ----------
    fig : matplotlib.Axes
        Target figure.
    positions : (N,2)-shaped array
        Position grid in normalized coordinates,
        i.e. [-1, 1].
    distortions : (N,2)-shaped array
        Distortions in arcseconds.
    posscale : float
        Scale of normalized positions to arcseconds.
    maxorder : int
        Maximum order for Legendre fits.
    """
    # RMS residuals over max order
    # Last iteration at maxorder
    orders = list(range(1, maxorder+1))
    rmsperorder = []
    for order in orders:
        print("Fitting {:d}. order Legendre polynomial".format(order))
        vf = PolyVectorField(Legendre(order))
        params, residuals, resvar = polyfit_svd(vf, positions, distortions)
        rmsperorder.append(resvar)
    rmsperorder = np.array(rmsperorder)

    print("In {:d}. order fit:".format(order))
    vf.set_params(params)
    indivorder = list(orders)
    rmspowerinorder = []
    for order in indivorder:
        subvf = vf.make_single_degree_subpoly(order)
        model = subvf.model(positions)
        rmspower = np.sqrt(np.sum(model**2)/model.size)
        rmspowerinorder.append(rmspower)
        print("  {:8.3g} RMS distortions in {:d}. degree terms".format(
            rmspower, order))
    rmspowerinorder = np.array(rmspowerinorder)

    # Actual plotting.
    gs = GridSpec(2, 2, width_ratios=[2, 2], height_ratios=[2, 2])
    gs.update(left=0.1, right=0.9, bottom=0.07, top=0.9,
              wspace=0.25, hspace=0.2)

    # distortion map
    ax1 = fig.add_subplot(gs[0, 0])
    plot_distortions(
        ax1, positions, distortions * 1e3,
        positionunits=(posscale, "arcsec"),
        arrowunits="mas")
    capname = name[:1].upper() + name[1:]
    fig.text(0.28, 0.94, capname+" pattern",
             horizontalalignment='center')

    # residual map of max order fit
    # (creates its own gridspec)
    plot_residuals(
        fig, positions, residuals * 1e6,
        positionunits=(posscale, "arcsec"),
        arrowunits="uas",
        left=0.52, right=0.9, bottom=0.52, top=0.9,
        wspace=0.05, hspace=0.05)
    fig.text(0.72, 0.94,
             "Residual map of {:d}. order fit".format(maxorder),
             horizontalalignment='center')

    # RMS residuals over max order
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.scatter(orders, rmsperorder * 1e6, marker='D')
    vmin = 0.7 * np.min(rmsperorder * 1e6)
    vmax = 1.3 * np.max(rmsperorder * 1e6)
    ax3.set_ylim(vmin, vmax)
    ax3.set_yscale('log')
    ax3.set_ylabel("RMS residuals / uas")
    ax3.set_xlabel("max order of fit")
    ax3.xaxis.set_major_locator(MaxNLocator(integer=True))

    # RMS per individual order of max order fit
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.scatter(indivorder, rmspowerinorder * 1e6, marker='D')
    vmin = 0.7 * np.min(rmspowerinorder * 1e6)
    vmax = 1.5 * np.max(rmspowerinorder * 1e6)
    # cut off for small values, indicate with downward arrow
    if vmin < 1e-2:
        vmin = 1e0
        vmax = max(vmax, 1e2 * vmin)
    for o, v in zip(indivorder, rmspowerinorder):
        if v * 1e6 < vmin:
            ax4.arrow(o, 2 * vmin, 0, -1*vmin, color='black',
                      width=0.01, head_width=0.1, head_length=0.1,
                      length_includes_head=True)
    ax4.set_ylim(vmin, vmax)
    ax4.set_yscale('log')
    uncapname = name[:1].lower() + name[1:]
    ax4.set_ylabel("RMS power of "+uncapname+" / uas")
    ax4.set_xlabel("total degree of terms")
    ax4.xaxis.set_major_locator(MaxNLocator(integer=True))
