
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.

"""Helpers for Zemax files."""


import numpy as np
import pandas as pd

import mascado.utility.affine as affine


def load_grid_data(fname, encoding='latin1', footer=7, nonsquareokay=False):
    """Load Zemax Grid Distortion Data from file.

    Parameters
    ----------
    fname : string
        File name.
    encoding : string
        File encoding, defaults to `latin1`.
    footer : int
        Number of non-empty lines that are not data rows at the
        end of the file
    nonsquareokay : bool
        Set to ``True`` to avoid exception if the number of
        points in the data set is not a square number.

    Returns
    -------
    pandas.DataFrame
        DataFrame with 9 columns: i, j, x-field (degree), y-field (degrees),
        r-field (degrees), x-predicted (mm), y-predicted (mm),
        x-real (mm), y-real (mm).

    Raises
    ------
    ValueError
        If the number of points is not a square number and
        ``nonsquareokay`` is ``False``.
    """
    data = np.genfromtxt(
        fname, encoding=encoding,
        skip_header=13, skip_footer=footer,
        usecols=list(range(9)))
    if not nonsquareokay and data.shape[0]**0.5 % 1 != 0:
        raise ValueError("Input file does not contain a square number of"
                         " data points: "+fname)
    df = pd.DataFrame(data, columns=[
        'i', 'j', 'x-field', 'y-field', 'r-field',
        'x-predicted', 'y-predicted', 'x-real', 'y-real'])
    return df


def load_grid_data_B(fname, encoding='latin1', nonsquareokay=False):
    """
    Parameters
    ----------
    fname : string
        File name.
    encoding : string
        File encoding, defaults to ``latin1``.
    nonsquareokay : bol
        Complain if the number of points in the data set is not
        a square number.

    Returns
    -------
    pandas.DataFrame
        DataFrame with 5 columns: n (running number),
        x-field (degree), y-field (degree), x-real (mm), y-real (mm).

    Raises
    ------
    ValueError
        If the number of points is not a square number and
        ``nonsquareokay`` is ``False``.
    """
    data = np.genfromtxt(
        fname, encoding=encoding,
        skip_header=8)
    if nonsquareokay and data.shape[0]**0.5 % 1 != 0:
        raise ValueError("Input file does not contain square number of"
                         " data points: "+fname)
    df = pd.DataFrame(data, columns=[
        'n', 'x-field', 'y-field', 'x-real', 'y-real'])
    return df


def load_grid_data_variant(variantkey, fname, **kwargs):
    """Load distortion data from files in various formats.

    Scripts and macros in Zemax produce different outputs that may be
    added here as input format variants for convenience.  Different
    variants are enumerated by upper case letters, where ``A`` is the
    default Zemax Grid Distortion export file.

    Parameters
    ----------
    variantkey : str
        Upper case letter denoting format variant.
    fname : str
        Path to input file.
    kwargs : keyword arguments
        Additional arguments for format loaders.

    Returns
    -------
    :class:`pandas.DataFrame`
        Which should have at least the columns:
        x-field (degree), y-field (degree), x-real (mm), y-real (mm).
        Check variant loader documentations.
    """
    if variantkey == 'A':
        return load_grid_data(fname, **kwargs)
    if variantkey == 'B':
        return load_grid_data_B(fname, **kwargs)


def have_same_grid(cats):
    """Check if all grid catalogs have the same grid in field coordinates.

    Parameters
    ----------
    cats : list of pandas.DataFrame
        Catalogs as described by doc of ``load_grid_data()``.

    Returns
    -------
    bool
    """
    xfields = [c.loc[:, 'x-field'] for c in cats]
    yfields = [c.loc[:, 'y-field'] for c in cats]
    return all(xfields[0].equals(f) for f in xfields[1:]) \
        and all(yfields[0].equals(f) for f in yfields[1:])


def distortions_on_sky(cats, platescale=None, scale=1):
    r"""Get distortions and normalized positions on-sky.

    By default, an affine transform is used for the transformation
    from focal plane to sky.  The affine transform is the
    least-squares solution for the first catalog and the same trafo is
    applied to all catalogs, so that a drift in linear terms is
    conserved.  If ``platescale`` is passed, no affine transform, but
    a fixed scale is used.

    Normalized positions are calculated by shifting and scaling the
    position catalog into the domain :math:`[-1, 1]\times[-1, 1]`.

    Parameters
    ----------
    cats : list of pandas.DataFrame
        Catalogs as described by doc of ``load_grid_data()``.
    platescale : float or None
        Optional fixed plate scale in arcseconds per mm.
    scale : float
        Additional, dimensionless scale apply applied to every position.

    Returns
    -------
    atrafo : (3, 3)-shaped array
        Affine transformation applied to translate from focal plane
        (real coordinates) to sky (field coordinates).
    posnormscale : float
        Scale in ``arcsecond`` for normalized positions.
    positions : (N, 2)-shaped array
        Dimensionless, normalized positions.
    distortions : list of (N, 2)-shaped arrays
        Distortions of each catalog in arcseconds.

    Raises
    ------
    ValueError
        If not all catalogs have the same reference grid
        (field coordinates).
    """
    if not have_same_grid(cats):
        raise ValueError("All catalogs are required to have"
                         " the same reference grid.")

    # get reference catalog
    index = cats[0].index
    refpos = cats[0].loc[index, ['x-field', 'y-field']].as_matrix()  # degree
    refpos = refpos * scale  # apply additional scale
    refpos = refpos * 3600   # convert degree to arcsec

    # get distortions on sky
    realpos = [df.loc[index, ['x-real', 'y-real']].as_matrix() for df in cats]
    if platescale is not None:
        atrafo = np.array([
            [platescale * scale, 0,                  0],
            [0,                  platescale * scale, 0],
            [0,                  0,                  1]])
    else:
        atrafo = affine.affine_lstsq(realpos[0], refpos)
    skypos = [affine.affine_trafo(rpos, atrafo) for rpos in realpos]
    distortions = [spos - refpos for spos in skypos]

    # normalize positions
    posmin, posmax = np.min(refpos), np.max(refpos)
    posscale = (posmax - posmin) / 2
    posshift = (posmax + posmin) / 2
    positions = (refpos - posshift) / posscale

    return atrafo, posscale, positions, distortions
