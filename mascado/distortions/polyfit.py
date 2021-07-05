
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.

"""Functions for least-squares fits of polynomial vector fields.

Interface follows functional paradigm.

Some functions take upper bounds, for example for outlier counts and
condition numbers.  The defaults for these bounds are chosen rather
restrictively so you are aware of potential problems, when you have to
turn them up.

Examples
--------
# TODO the example is broken: Need a lot more positions/vectors to make polynomial overdetermined
#  if you have an understanding of how to use this, please adapt and remove comments to re-enable
# >>> from mascado.distortions.polynomials import (
# ...     PolyVectorField, Legendre)
# >>> vf = PolyVectorField(Legendre(3))
# >>> grid = np.array([[0, 0], [1, 1], [0, 1], [1, 0], [0.5, 0.5]])
# >>> vectors = np.array([[1, 1], [0, 1], [0.5, 3], [-1, -1], [100, 3]])
# >>> params, inliers, resvar = polyfit_svd_iterative(
# ...     vf, grid, vectors, sigmas=0.5, maxoutliers=1, iterinfo='  ')
#   5 2D data points, 6 parameters
#   0 outliers before first iteration, 1 max
#     1. iteration: resvar 7.98e+03, outliers: 1 new, 0 old
#     2. iteration: resvar 1.12, outliers: 0 new, 1 old
#   Iterations: 2, outliers: 1
#   Residual variance: 1.12
# >>> vf.set_params(params)

"""

import warnings
import numpy as np


def polyfit_svd(vf, positions, vectors, sigmas=None,
                maxcondition=1e2, info='  '):
    r"""Fit vector field to data using SVD.

    Does not use or alter internal parameters of vector field.

    Parameters
    ----------
    vf : :class:`mascado.distortions.polynomials.PolyVectorField`
        Model of vector field providing the Vandermonde matrix.
    positions : (N,2)-shaped array
        x,y positions of data points.
    vectors : (N,2)-shaped array
        x,y values at data points.
    sigmas : float or (N, 1) array or (N, 2) array or None
        Uncertainties for vectors.  If one-dimensional, the same
        weight is used for x and y component.
    maxcondition : float
        Function raises a RuntimeError if the condition value of the
        weighted Vandermonde matrix is larger than this argument.
    info : whitespace or False
        Indent for info printed to console.
        If False, nothing is printed.

    Returns
    -------
    params : array
        Parameter list for vector field.
    residuals : (N,2)-shaped array
        Residual vectors.
    resvar : float
        Residual variance, aka. reduced chi squared.
        NaN if the problem is not over-determined.

    Raises
    ------
    RuntimeError
        If the weighted Vandermonde matrix is conditioned badly.

    Notes
    -----
    Formulate least-squares problem as matrix equation using
    Vandermonde matrices :math:`V_x,V_y` of polynomials in
    x and y direction at given data points.
    (Polynomials are linear in their coefficients.)
    :math:`q` are the vector components and :math:`a` the polynomial
    coefficients:

    .. math:: q \approx X \cdot a

    .. math::
        q =& [q_{1x}, ..., q_{Nx}, q_{1y}, ..., q_{Ny}]^T \\
        a =& [a_{0x}, ..., a_{mx}, a_{0y}, ..., a_{my}]^T \\
        X =& \begin{bmatrix}
          V_x & \mathbf{0} \\
          \mathbf{0} & V_y \end{bmatrix} \,.

    Minimize the sum of squared residuals

    .. math:: \lVert r \rVert_2^2 = \lVert q - X \cdot a \rVert_2^2

    using the pseudo inverse of :math:`W \cdot X`, obtained by SVD:

    .. math:: a = X^+ \cdot q \,.

    Check this fine
    `StackExchange answer <https://math.stackexchange.com/a/2173715>`_
    for details.

    In the weighted case, this becomes

    .. math:: W \cdot q \approx W \cdot X \cdot a

    .. math:: \lVert r \rVert_2^2 = \lVert W \cdot q -
          W \cdot X \cdot a \rVert_2^2
    .. math:: a = (W \cdot X)^+ \cdot W \cdot q \,.

    :math:`W` is the weighting matrix

    .. math:: W = \operatorname{diag}\left(
          \frac{1}{\sigma_1}, ..., \frac{1}{\sigma_N}\right)\,.

    The values are not squared, because the matrix is applied
    inside the L2-norm.
    """
    if positions.shape[0] != vectors.shape[0]:
        raise ValueError("Mismatching number positions and vectors.")
    assert positions.shape[1] == 2
    assert vectors.shape[1] == 2

    n = positions.shape[0]
    dof = 2*n - vf.paramcount
    if dof <= 0:
        warnings.warn("Problem is not over-determined.", RuntimeWarning)
    if info is not False:
        determination = "over-determined" if dof > 0 else \
                        "well-determined" if dof == 0 else \
                        "under-determined"
        print("{:s}{:d} 2D data points, {:d} parameters: {:s}".format(
            info, n, vf.paramcount, determination))

    y = np.concatenate([vectors[:, 0], vectors[:, 1]])
    X = vf.vandermonde(positions)

    if sigmas is None:
        W = np.identity(2*n)
    elif np.isscalar(sigmas):
        W = np.diag(np.repeat(1/sigmas, y.size))
    elif sigmas.ndim == 1:
        if sigmas.shape[0] != n:
            raise ValueError("Mismatching number of sigmas and positions.")
        W = np.diag(1/np.concatenate([sigmas, sigmas]))
    else:
        if sigmas.shape != (n, 2):
            raise ValueError("Mismatching number of sigmas and positions.")
        W = np.diag(1/np.concatenate([sigmas[:, 0], sigmas[:, 1]]))

    # `@` is NumPys operator for the matrix multiplication of arrays
    U, S, VT = np.linalg.svd(W @ X, full_matrices=False)

    # check conditioning
    smax, smin = S[0], S[-1]
    condition = smax / smin
    if info is not False:
        print("{:s}Condition number: {:.3g} (max {:.3g})".format(
            info, condition, maxcondition))
    if condition > maxcondition:
        raise RuntimeError("Problem is ill-conditioned: {}".format(condition))

    # compute pseudo inverse and coefficients
    WXinv = VT.T @ np.diag((1/S)) @ U.T
    params = WXinv @ W @ y

    # residual variance
    residuals = y - X @ params
    resvar = np.sum((W @ residuals)**2) / dof if dof > 0 else np.nan
    if info is not False and not np.isnan(resvar):
        print("{:s}Residual variance: {:.3g}".format(
            info, resvar))

    residuals = np.reshape(residuals, (n, 2), order='F')
    return params, residuals, resvar


def polyfit_svd_iterative(vf, positions, vectors, sigmas, maxoutliers,
                          sigmacut=5, initialinliers=None,
                          maxcondition=1e2,
                          info='  ', iterinfo=False):
    """Fit vector field with iterative outlier rejection.

    Outlier rejection is realized through sigma clipping.  Because the
    outliers are introducing systematic errors ("bumps") in the
    preliminary solutions, good data points might be rejected, if we
    would just cut all residuals larger than `sigmacut` times
    `sigmas`.  Therefore the `sigmacut` is scaled with the square root
    of the residual variance, so, initially, the applied `sigmacut` is
    larger than given.  As the residual variance converges towards 1,
    if all outliers are found, the threshold converges to the given
    `sigmacut`.

    If no `sigmas` are passed, the root mean square over all residuals
    is used instead (not scaled by the residual variance, because
    that's already in the RMS of residuals).

    The set of outliers increases with every iteration, until no more
    outliers are found and the algorithm halts.  Outliers from
    previous iterations are not reconsidered in subsequent iterations,
    because that might cause an infinite loop.

    Parameters
    ----------
    vf : :class:`mascado.distortions.polynomials.PolyVectorField`
        2D vector field.  Internal parameters are neither used nor changed.
    positions : (N, 2)-shaped array
    vectors : (N, 2)-shaped array
    sigmas : float or { (N,), (N, 2) }-shaped array or None
    maxoutliers : int
        Function raises RuntimeError if the number of outliers
        increases above this number.
    sigmacut : float
        Multiplier for standard deviation specifying
        the upper bound for inliers.
    initialinliers : (N,)-shaped bool array
        Data points, where this array is `False`, are
        not used in any iteration.
    maxcondition : float
        Passed on to `polyfit_svd()`.
    info : whitespace or False
        Indent for info printed to console.
        If `False`, nothing is printed.
    iterinfo : whitespace or False
        Additional indent for info printed in every iteration.
        Nothing is printed by default.

    Returns
    -------
    params : float array
        Parameter list for vector field.
    inliers : (N,)-shaped bool array
        Mask for inliers.
    residuals : (inlier count, 2)-shaped array
        Residual vectors of inliers.
    resvar : float
        Residual variance of inliers.

    Raises
    ------
    RuntimeError
        If the maximum outlier count is exceeded.
    """
    n = positions.shape[0]
    if 2*n <= vf.paramcount:
        raise ValueError("Outlier rejection not possible"
                         " if not over-determined.")

    # mask accumulates outliers of all runs
    # False for outliers, True for inliers
    if initialinliers is not None:
        assert initialinliers.size == n
        mask = initialinliers.copy()
    else:
        mask = np.full((n,), True)

    if info is not False:
        print("{:s}{:d} 2D data points, {:d} parameters".format(
            info, n, vf.paramcount))
        print("{:s}{:d} outliers before first iteration, {:d} max".format(
            info, np.count_nonzero(~mask), maxoutliers))

    if sigmas is not None:
        if np.isscalar(sigmas):
            sigmas = np.full((n, 2), sigmas)
        elif sigmas.ndim == 1:
            sigmas = np.stack([sigmas, sigmas], axis=1)

    nit = 0
    resvar = np.nan
    while True:
        nit += 1
        params, residuals, resvar = polyfit_svd(
            vf, positions[mask], vectors[mask],
            sigmas[mask] if sigmas is not None else None,
            maxcondition=maxcondition, info=False)
        if sigmas is None:
            resrms = np.sqrt(np.sum(residuals**2)/residuals.size)
            inliers = np.logical_and(
                np.absolute(residuals[:, 0]) < resrms * sigmacut,
                np.absolute(residuals[:, 1]) < resrms * sigmacut)
        else:
            deviation = np.absolute(residuals / sigmas[mask])
            inliers = np.logical_and(
                deviation[:, 0] < np.sqrt(resvar) * sigmacut,
                deviation[:, 1] < np.sqrt(resvar) * sigmacut)

        if iterinfo is not False:
            print("{:s}{:s}{:d}. iteration: resvar {:.3g},"
                  " outliers: {:d} new, {:d} old".format(
                      info, iterinfo, nit, resvar,
                      np.count_nonzero(~inliers), np.count_nonzero(~mask)))

        if np.count_nonzero(~inliers) == 0:
            # no new outliers in this iteration
            # stable solution reached
            break

        # accumulate new outliers in mask
        outlieridxs = np.where(mask)[0][~inliers]
        mask[(outlieridxs,)] = False
        if np.count_nonzero(~mask) > maxoutliers:
            raise RuntimeError(
                "Iterative polyfit failed with too many outliers:"
                " {:d} of {:d}, max {:d}".format(
                    np.count_nonzero(~mask), n, maxoutliers))

    if info is not False:
        print("{:s}Iterations: {:d}, outliers: {:d}".format(
            info, nit, np.count_nonzero(~mask)))
        print("{:s}Residual variance: {:.3g}".format(info, resvar))

    return params, mask, residuals[inliers], resvar
