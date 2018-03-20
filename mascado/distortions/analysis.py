
import numpy as np
import pandas as pd

from mascado.distortions.polynomials import PolyVectorField, Legendre
from mascado.distortions.polyfit import polyfit_svd


def analyze_residuals_over_order(
        positions, distortions, maxorder, minorder=0, poly=Legendre,
        distortionunits="arcseconds", info='  '):
    """Get RMS of residuals after nth order fit.

    Parameters
    ----------
    positions : (N, 2)-shaped array
        Normalized positions.
    distortions : (N, 2)-shaped array
        Distortions.
    maxorder : int
        End of iteration.
    minorder : int
        Start of iteration.
    poly : :class:`mascado.distortions.polynomials.Polynomial` subclass
        Class of polynomials to use (not an instance).  By default,
        Legendre polynomials are used.
    distortionunits : str
        Units of distortions.  Used for console output.
    info : str or False
        Indent of console output.
        Set to ``False`` to suppress output.

    Returns
    -------
    pandas.DataFrame
        DataFrame with two columns: ``order`` and ``resrms``.
    """
    subinfo = info + '  ' if info is not None else None
    orders = list(range(minorder, maxorder+1))

    rmsperorder = []
    for order in orders:
        if info is not False:
            print("{:s}Fitting {:d}. order {:s} polynomial".format(
                info, order, poly.__name__))
        vf = PolyVectorField(Legendre(order))
        _, residuals, _ = polyfit_svd(vf, positions, distortions,
                                      info=subinfo)
        resrms = np.sqrt(np.sum(residuals**2)/residuals.size)
        rmsperorder.append(resrms)
        if info is not False:
            print("{:s}Residual RMS: {:.3g} {:s}".format(
                subinfo, resrms, distortionunits))

    tab = pd.DataFrame([orders, rmsperorder], index=['order', 'resrms']).T
    return tab


def analyze_contributions_over_order(
        vf, minorder=0, maxorder=None, gridsize=100,
        units="arcseconds", info='  '):
    r"""Get contributions of the terms of a specific degree.

    The contribution is measured as RMS of the vector field over the
    normalized coordinates :math:`[-1,1]\times[-1,1]`.  The model is
    evaluated on a regular grid to calculate the RMS value.  The grid
    has to be well sampled for the result to be correct.  For higher
    orders, a larger ``gridsize`` is needed.

    Parameters
    ----------
    vf : :class:`mascado.distortions.polynomials.PolyVectorField`
        Vector field with internal parameters.
    minorder : int
    maxorder : int or None
        If None, the order of the vector field is used.
    gridsize : int
        Number of points along each axis in normalized coordinates.
    units : str
        Units of vector field.  Used for console output.
    info : str or False
        Indent of console output.
        Set to ``False`` to suppress output.

    Returns
    -------
    pandas.DataFrame
        DataFrame with two columns: ``order`` and ``modelrms``.
    """
    if maxorder is None:
        maxorder = vf.get_degree()
    orders = list(range(minorder, maxorder+1))

    gridside = np.linspace(-1, 1, gridsize)
    xx, yy = np.meshgrid(gridside, gridside, indexing='ij')
    grid = np.stack([xx.ravel(), yy.ravel()], axis=1)

    rmsperdegree = []
    for order in orders:
        subvf = vf.make_single_degree_subpoly(order)
        model = subvf.model(grid)
        modelrms = np.sqrt(np.sum(model**2)/model.size)
        rmsperdegree.append(modelrms)
        if info is not False:
            print("{:s}{:8.3g} {:s} RMS distortions in {:d}. degree terms"
                  "".format(info, modelrms, units, order))

    tab = pd.DataFrame([orders, rmsperdegree], index=['order', 'modelrms']).T
    return tab
