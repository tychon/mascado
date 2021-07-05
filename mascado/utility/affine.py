
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.

"""Provide 2D affine transformations.

Interface follows a functional paradigm.  Runs on plain NumPy arrays
and functions with suffix `_cat` munch Pandas DataFrames.
"""


import numpy as np
import pandas as pd


def affine_lstsq(origin, target):
    r"""Calculate optimal affine transformation in the least-squares sense.

    Parameters
    ----------
    origin : (n,2)-shaped array
        x,y points in the domain of the transformation.
    target: (n,2)-shaped array
        x,y points in the target of the transformation.

    Returns
    -------
    (3,3)-shaped array
        Matrix representation of affine trafo.

    Notes
    -----
    The affine transformation is expressed as matrix :math:`A` in
    homogeneous coordinates:

    .. math::

        \mathtt{target} &= A\cdot\mathtt{origin} \\
        \mathbf x^\prime &= A\cdot \mathbf x \\
        \begin{bmatrix}
          x' \\ y' \\ 1
        \end{bmatrix}
        &= \begin{bmatrix}
          a_{11} & a_{12} & a_{13} \\
          a_{21} & a_{22} & a_{23} \\
          0      & 0      & 1
        \end{bmatrix}
        \begin{bmatrix}
          x \\ y \\ 1
        \end{bmatrix}

    Rewrite as

    .. math::

        \begin{bmatrix}
          x_1' \\ \vdots \\ x_n' \\
          y_1' \\ \vdots \\ y_n'
        \end{bmatrix}
        = \begin{bmatrix}
          x_1 & y_1        & 1  & {}  & {}         & {} \\
          {}  & \vdots     & {} & {}  & \mathbf{0} & {} \\
          x_n & y_n        & 1  & {}  & {}         & {} \\
          {}  & {}         & {} & x_1 & y_1        & 1 \\
          {}  & \mathbf{0} & {} & {}  & \vdots     & {} \\
          {}  & {}         & {} & x_n & y_n        & 1 \\
        \end{bmatrix}
        \begin{bmatrix}
          a_{11} \\ a_{12} \\ a_{13} \\
          a_{21} \\ a_{22} \\ a_{23} \\
        \end{bmatrix}

    or for short

    .. math:: \mathbf y = X \cdot \mathbf a\,.

    Use pseudo-inverse :math:`X^+` for least-squares solution:

    .. math:: \mathbf a = X^+ \cdot \mathbf y
    """
    n = origin.shape[0]  # number of points
    if target.shape[0] != n:
        raise ValueError("Different number of data points "
                         "in origin and target.")
    assert origin.shape == target.shape == (n, 2)

    y = np.concatenate([target[:, 0], target[:, 1]])
    X = np.zeros(shape=(2*n, 6), dtype=float)
    X[0:n, 0:2] = X[n:2*n, 3:5] = origin
    X[0:n, 2] = X[n:2*n, 5] = 1

    Xinv = np.linalg.pinv(X)
    a = Xinv @ y

    A = np.zeros(shape=(3, 3), dtype=float)
    A[0, :] = a[:3]
    A[1, :] = a[3:]
    A[2, 2] = 1
    return A


def affine_lstsq_cat(origin, target, x='x', y='y', ignore_index=False):
    """Like `affine_lstsq` with pandas DataFrame.

    Parameters
    ----------
    origin : (n,2)-shaped array
    target : (n,2)-shaped array
    x : string
        Column name of x components of data points.
    y : string
        Column name of y components of data points.
    ignore_index : bool
        By default, data points are referenced by their index.
        Set to `True` to associate by row position.

    Raises
    ------
    ValueError
        If the catalogs are not of the same length or
        the indices don't match.
    """
    if ignore_index:
        if origin.shape[0] != target.shape[0]:
            raise ValueError("Cannot merge catalogs of different length.")
        origin = origin.loc[:, [x, y]].values
        target = target.loc[:, [x, y]].values
    else:
        if not origin.index.equals(target.index):
            raise ValueError("Cannot merge on unequal indices.")
        idx = origin.index
        origin = origin.loc[idx, [x, y]].values
        target = target.loc[idx, [x, y]].values
    return affine_lstsq(origin, target)


def affine_trafo(points, trafo):
    """Apply affine transformation.

    Parameters
    ----------
    points : (N,2)-shaped array
        x,y points in the domain of the transformation.
    trafo : (3,3)-shaped array
        Affine transformation in matrix form.

    Returns
    -------
    (N,2)-shaped array
        Transformed points.
    """
    assert len(points.shape) == 2
    assert points.shape[1] == 2
    linear = trafo[:2, :2]  # linear component of trafo
    offset = trafo[:2, 2]   # offset of trafo
    # broadcasts dot product over list of vectors
    return np.dot(linear, points.T).T + offset


def affine_trafo_cat(cat, trafo, x='x', y='y'):
    """Like `affine_trafo` with pandas DataFrame.

    Parameters
    ----------
    cat : pandas.DataFrame
        Catalog with points.
    trafo : (3,3)-shaped array
        Affine transformation in matrix form.
    x : string
        Column name of x components of data points.
    y : string
        Column name of y components of data points.

    Returns
    -------
    pandas.DataFrame
        DataFrame with same index and
        two columns named like input columns.

    Examples
    --------

    >>> trafo = np.array([[1, 0, 3], [2, 1, 0], [0, 0, 1]])
    >>> trafo
    array([[1, 0, 3],
           [2, 1, 0],
           [0, 0, 1]])
    >>> cat = pd.DataFrame(
    ...     [[0, 0, 'Anja'], [1, 0, 'Bert'], [2, 2, 'Chris']],
    ...     index=['A', 'B', 'C'], columns=['x', 'y', 'extra'])
    >>> cat
       x  y  extra
    A  0  0   Anja
    B  1  0   Bert
    C  2  2  Chris
    >>> cat2 = affine_trafo_cat(cat, trafo)
    >>> cat2
       x  y
    A  3  0
    B  4  2
    C  5  6

    To integrate back into the original catalog do

    >>> cat[cat2.columns] = cat2
    >>> cat
       x  y  extra
    A  3  0   Anja
    B  4  2   Bert
    C  5  6  Chris
    """
    points = cat.loc[:, [x, y]].values
    points = affine_trafo(points, trafo)
    return pd.DataFrame(points, index=cat.index, columns=[x, y])
