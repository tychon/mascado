
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.

"""Provide some plotting functions useful for distortion solutions."""


import numpy as np
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm

from mascado.distortions.polynomials import PolyVectorField, Legendre
from mascado.distortions.polyfit import polyfit_svd
import mascado.distortions.analysis as analysis
import mascado.distortions.powerspectrum as psd


def set_fontsize(small=8, medium=10, bigger=12):
    """Set font size of ticks, axes labels and titles."""
    plt.rc('font', size=medium)          # controls default text sizes
    plt.rc('axes', titlesize=medium)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)     # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)     # fontsize of the tick labels
    plt.rc('legend', fontsize=medium)    # legend fontsize
    plt.rc('figure', titlesize=bigger)   # fontsize of the figure title


def plot_distortions(ax, positions, distortions,
                     inliers=None, limits=(-1.1, 1.1),
                     positionunits=(1, "normalized"), positionunits2=None,
                     arrowunits="arcsec", keylen=None,
                     keypos=(0.5, 1.04), dotsize=2):
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
        If the ``keylen`` is not given and the x component of ``keypos``
        inside the axes (0 to 1), a note about the RMS is appended to
        the key label.
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

    keylabel = None
    if keylen is None:
        rmslen = np.sqrt(np.sum(distortions[inliers]**2)
                         / distortions[inliers].size)
        if np.isclose(rmslen, 0, rtol=1e-11, atol=1e-11):
            raise ValueError("Distortions too close to zero.")
        keylen = rmslen * 2
        keylen = np.round(keylen, -int(np.floor(np.log10(keylen)))+1)
        if 0 < keypos[0] < 1:
            keylabel = r"{:.2g} {:s} = $2\cdot$RMS".format(keylen, arrowunits)
    if not keylabel:
        keylabel = r"{:.2g} {:s}".format(keylen, arrowunits)
    ax.quiverkey(quiver, *keypos, keylen, label=keylabel,
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
                   keypos=(1.2, 1.2), dotsize=1, **kwgridspec):
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
    >>> import matplotlib.pyplot as plt
    >>> positions = np.array([[0,0],[1,1],[1,0],[0,1]], dtype=float)
    >>> residuals = np.array([[0,0],[1,1],[1,0],[0,1]], dtype=float)
    >>> posscale = 1.
    >>>
    >>> fig=plt.figure()
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

    gs = GridSpec(3, 3, width_ratios=[1, 4, 1], height_ratios=[1, 4, 1])
    gs.update(**kwgridspec)

    ax = fig.add_subplot(gs[1, 1])
    plot_distortions(ax, positions, residuals, inliers, limits,
                     positionunits, None, arrowunits, keylen,
                     keypos, dotsize)
    ax.tick_params(axis='both', which='both', right=True, top=True,
                   labelleft=False, labelbottom=False)
    ax.set_xlabel('')
    ax.set_ylabel('')

    axmargxdx = fig.add_subplot(gs[1, 0], sharey=ax)
    axmargxdx.tick_params(axis='both', which='both', left=True, right=False,
                          labelleft=True, labelright=False)
    axmargxdx.scatter(residuals[inliers, 0], positions[inliers, 1] * ps,
                      s=dotsize**2, c='k', alpha=0.2)
    axmargxdx.set_ylabel("position / {:s}".format(positionunits[1]))
    axmargxdx.set_xlabel("res/{:s}".format(arrowunits))

    axmargxdy = fig.add_subplot(gs[1, 2], sharey=ax)
    axmargxdy.tick_params(axis='both', which='both', left=False, right=True,
                          labelleft=False, labelright=True,
                          top=True, bottom=False,
                          labeltop=True, labelbottom=False)
    axmargxdy.scatter(residuals[inliers, 1], positions[inliers, 1] * ps,
                      s=dotsize**2, c='k', alpha=0.2)

    axmargydx = fig.add_subplot(gs[0, 1], sharex=ax)
    axmargydx.tick_params(axis='both', which='both', left=True, right=False,
                          labelleft=True, labelright=False,
                          top=True, bottom=False,
                          labeltop=True, labelbottom=False)
    axmargydx.scatter(positions[inliers, 0] * ps, residuals[inliers, 0],
                      s=dotsize**2, c='k', alpha=0.2)

    axmargydy = fig.add_subplot(gs[2, 1], sharex=ax)
    axmargydy.tick_params(axis='both', which='both', left=False, right=True,
                          labelleft=False, labelright=True,
                          top=False, bottom=True,
                          labeltop=False, labelbottom=True)
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
                       minorder=0, name="distortions", poly=Legendre,
                       maxcondition=1e2):
    """Run analysis of distortion pattern and display in four panels.

    The four panels are:

    1. Arrow plot of distortion pattern itself.
    2. Arrow plot with marginalized distributions of residuals of
       fit at `maxorder`.
    3. RMS residuals over maximum order of polynomials.
       Maximum order from 1 to `maxorder`.
    4. RMS power in each degree of polynomial of
       fit at `maxorder`.

    Since an affine transformation was applied to calculate the distortions,
    the zeroth and first order terms of the polynomial solution are incomplete
    and are plotted with gray data points.

    Some information is printed to the console.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Target figure.
    positions : (N,2)-shaped array
        Position grid in normalized coordinates,
        i.e. [-1, 1].
    distortions : (N,2)-shaped array
        Distortions in arcseconds.
    posscale : float
        Scale of normalized positions to arcseconds.
    maxorder : int
        Maximum order for fits.
    name : str
        Name of distortions in plots
        (eg. ``"drift"`` for distortions difference).

    Returns
    -------
    vf : :class:`mascado.distortions.polynomials.PolyVectorField`
        Fitted 2D vector field of ``maxorder``.
    """
    resrms = analysis.analyze_residuals_over_order(
        positions, distortions, maxorder, poly=poly, minorder=0,
        maxcondition=maxcondition, info='')

    print()
    print("In {:d}. order {:s} fit:".format(maxorder, poly.__name__))
    vf = PolyVectorField(poly(maxorder))
    params, residuals, _ = polyfit_svd(vf, positions, distortions,
                                       maxcondition=maxcondition)
    vf.set_params(params)
    modelrms = analysis.analyze_contributions_over_order(vf)

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
        left=0.55, right=0.85, bottom=0.52, top=0.9,
        wspace=0.05, hspace=0.05)
    fig.text(0.72, 0.94,
             "Residual map of {:d}. order fit".format(maxorder),
             horizontalalignment='center')

    # RMS residuals over max order
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.scatter(
        resrms.order[resrms.order <= 1],
        resrms.resrms[resrms.order <= 1] * 1e6,
        marker='D', color='silver')
    ax3.scatter(
        resrms.order[resrms.order > 1],
        resrms.resrms[resrms.order > 1] * 1e6,
        marker='D')
    vmin = 0.6 * np.min(resrms.resrms * 1e6)
    vmax = 1.4 * np.max(resrms.resrms * 1e6)
    ax3.set_ylim(vmin, vmax)
    ax3.set_yscale('log')
    ax3.set_ylabel("RMS residuals / uas")
    ax3.set_xlabel("max order of fit")
    ax3.xaxis.set_major_locator(MaxNLocator(integer=True))

    # RMS per individual order of max order fit
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.scatter(
        modelrms.order[modelrms.order <= 1],
        modelrms.modelrms[modelrms.order <= 1] * 1e6,
        marker='D', color='silver')
    ax4.scatter(
        modelrms.order[modelrms.order > 1],
        modelrms.modelrms[modelrms.order > 1] * 1e6,
        marker='D')
    vmin = 0.7 * np.min(modelrms.modelrms * 1e6)
    vmax = 1.5 * np.max(modelrms.modelrms * 1e6)
    # cut off for small values, indicate with downward arrow
    if vmin < 1e-2:
        vmin = 1e0
        vmax = max(vmax, 1e2 * vmin)
    for o, v in zip(modelrms.order, modelrms.modelrms):
        if v * 1e6 < vmin:
            ax4.arrow(o, 2 * vmin, 0, -1*vmin, color='black',
                      width=0.01, head_width=0.1, head_length=0.1,
                      length_includes_head=True)
    ax4.set_ylim(vmin, vmax)
    ax4.set_yscale('log')
    uncapname = name[:1].lower() + name[1:]
    ax4.set_ylabel("RMS of "+uncapname+" model / uas")
    ax4.set_xlabel("total degree of terms in {:d}. order fit".format(maxorder))
    ax4.xaxis.set_major_locator(MaxNLocator(integer=True))

    return vf


def plot_spectrum(ax, spectrum, smin, smax, caption):
    """Plot power spectrum as image with correct frequency scale.

    Parameters
    ----------
    ax : matplotlib.Axes
        Target axes for plotting
    spectrum : (N, N)-shaped float array
        Power spectrum
    smin, smax : float
        Minimum and maximum values for color map.
    caption : str
        Plot title

    Returns
    -------
    im
        matplotlib imshow result
    """
    fmin, fmax = -spectrum.shape[0]//2+1, spectrum.shape[0]//2
    im = ax.imshow(spectrum.T, origin='lower', cmap='Purples',
                   norm=LogNorm(vmin=smin, vmax=smax),
                   extent=(fmin-0.5, fmax+0.5, fmin-0.5, fmax+0.5))
    ax.text(0.5, 1.05, caption, horizontalalignment='center',
            transform=ax.transAxes)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("spatial frequency [1/FOV]")
    ax.set_ylabel("spatial frequency [1/FOV]")
    return im


def make_psd_analysis(fig, vf, posscale, name="distortions"):
    """Run analysis of power spectrum of distortion pattern.

    The figure is filled with four panels and one colorbar.  The upper
    two panels contain the unbinned 2D power spectrum for the x- and
    y-components of the vector field with logarithmic colorbar.

    The lower left panel shows the binned power spectra and the lower
    right panel the cumulative version of the lower left panel.

    Units on the power are unituitive and not clear.
    At least to the present me.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Target figure.
    vf : :class:`mascado.distortions.polynomials.PolyVectorField`
        Vector field to analyze.
    posscale : float
        Scale to convert from normalized coordinates to arcseconds.
    name : str
        Meta-name of data.

    """
    print()
    print("Calculating PSD...")
    sxx, syy = psd.psd_spectrums(vf, posscale)
    fx, fy = psd.psd_histogram(sxx), psd.psd_histogram(syy)
    fxc = psd.psd_histogram_cumulative(sxx)
    fyc = psd.psd_histogram_cumulative(syy)
    fxc = ((fxc - fxc[0]) / (fxc[-1] - fxc[0]))
    fyc = ((fyc - fyc[0]) / (fyc[-1] - fyc[0]))

    gs = GridSpec(2, 3, width_ratios=[2, 2, 0.2], height_ratios=[2, 2])
    gs.update(left=0.08, right=0.92, bottom=0.07, top=0.95,
              wspace=0.3, hspace=0.3)

    smax = max(np.max(sxx), np.max(syy))
    # lower limit: lower limit of spectrum
    # or at least down to 0.1 or at least 3 decades
    smin = max(min(0.1,
                   10**np.floor(np.log10(smax/1.1e2))),
               min(np.min(sxx),
                   np.min(syy)))
    ax1 = fig.add_subplot(gs[0, 0])
    im = plot_spectrum(ax1, sxx, smin, smax,
                       name[:1].upper() + name[1:] + " x-component")
    ax2 = fig.add_subplot(gs[0, 1])
    im = plot_spectrum(ax2, syy, smin, smax,
                       name[:1].upper() + name[1:] + " y-component")
    cax = fig.add_subplot(gs[0, 2])
    color = fig.colorbar(im, cax=cax, orientation='vertical')
    color.set_label("power")

    # histogram
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(fx, '.-', label="x-component")
    ax3.plot(fy, '.-', label="y-component")
    ax3.set_xlim(-0.5, fx.size-0.5)
    ax3.set_ylim(0, max(np.max(fx), np.max(fy))*1.1)
    ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax3.set_xlabel("spatial frequency [1/FOV]")
    ax3.set_ylabel("power")
    ax3.legend(loc='upper right')
    ax3t = ax3.twiny()
    ax3t.set_xlim([-0.5/(2*posscale), (fx.size-0.5)/(2*posscale)])
    ax3t.set_xlabel("spatial frequency [1/arcsec]", labelpad=5)

    # cumulative histogram
    ax4 = fig.add_subplot(gs[1, 1])

    def hline(hpos, text):
        ax4.axhline(hpos, color='gainsboro', zorder=0)
        ax4.text((fxc.size-0.7)/(2*posscale), hpos-5, text,
                 color='gainsboro', horizontalalignment='right')
    hline(50, "50% = -3dB loss")
    hline(75, "25% = -6dB loss")
    hline(90, "10% = -10dB loss")
    # hline(99, "1% = -20dB loss") # don't trust this
    hline(100, "no loss")

    ax4.plot(np.arange(fxc.size)/(2*posscale), fxc * 100, '.-',
             label="x-component")
    ax4.plot(np.arange(fxc.size)/(2*posscale), fyc * 100, '.-',
             label="y-component")
    ax4.set_xlim(0, (fxc.size-0.5)/(2*posscale))
    ax4.set_ylim(0, 102)
    ax4.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax4.set_xlabel("spatial frequency [1/arcsec]")
    ax4.set_ylabel("relative cumulative power [%]")
    ax4.legend(loc='lower right')

    # critically sampling pinhole spacing
    ticks = np.arange(fxc.size)[2::2]/(2*posscale)
    ticklabels = ["{:.3g}".format(0.5/freq) for freq in ticks]
    ax4t = ax4.twiny()
    ax4t.set_xticks(list(ticks)+[0])
    ax4t.set_xticklabels(ticklabels+[r'$\infty$'])
    ax4t.set_xlim(0, (fxc.size-0.5)/(2*posscale))
    ax4t.set_xlabel("critically sampling pinhole spacing [arcsec]",
                    labelpad=5)
