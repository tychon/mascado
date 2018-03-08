
"""Helpers for Zemax files."""


import numpy as np
import pandas as pd


def load_grid_data(fname, encoding='latin1'):
    """Load Zemax Grid Distortion Data from file.

    Parameters
    ----------
    fname : string
        File name.
    encoding : string
        File encoding, defaults to `latin1`.

    Returns
    -------
    pandas.DataFrame
        DataFrame with 9 columns: i, j, x-field, y-field, r-field,
        x-predicted, y-predicted, x-real, y-real.
    """
    data = np.genfromtxt(
        fname, encoding=encoding,
        skip_header=13, skip_footer=7,
        usecols=list(range(9)))
    df = pd.DataFrame(data, columns=[
        'i', 'j', 'x-field', 'y-field', 'r-field',
        'x-predicted', 'y-predicted', 'x-real', 'y-real'])
    return df


def have_same_grid(cats):
    """Check if all grid catalogs have the same grid in field coordinates.

    Parameters
    ----------
    cats : list of pandas.DataFrame

    Returns
    -------
    bool
    """
    xfields = [c.loc[:, 'x-field'] for c in cats]
    yfields = [c.loc[:, 'y-field'] for c in cats]
    return all(xfields[0].equals(f) for f in xfields[1:]) \
        and all(yfields[0].equals(f) for f in yfields[1:])
