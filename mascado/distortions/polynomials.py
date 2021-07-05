
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.

"""2D polynomials and vector fields.

Interface follows object oriented paradigm.

Important functions are documented in :class:`Polynomial` and
:class:`PolyVectorField` classes.

Currently supported polynomial types are Cartesian, Legendre, Zernike,
and Chebyshev.

Examples
--------
The standard usage to create a vector field is

>>> vf = PolyVectorField(Legendre(3))

set its parameter list

>>> vf.set_params(np.random.randn(vf.paramcount))

and evaluate it on a grid:

>>> grid = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
>>> values = vf.model(grid)

"""


from abc import ABC, abstractmethod
import numpy as np
from scipy.special import binom, factorial

from numpy.polynomial.polynomial import polyval2d, polyvander2d
from numpy.polynomial.legendre import legval2d, legvander2d
from numpy.polynomial.chebyshev import chebval2d, chebvander2d


class Polynomial (ABC):
    r"""Polynomials mapping from 2D to 1D.

    Parameters
    ----------
    degree : int

    Attributes
    ----------
    degree : int
        Maximum degree :math:`k` of this polynomial.
    paramcount : int
        Number of parameters :math:`m`,
        calculated as :math:`\binom{k+2}{2}`.

    Notes
    -----
    A polynomial is build from a 1D basis set :math:`P_i(x)`
    using the parameters :math:`a_{ij}` like this:

    .. math:: P(x, y) = \sum_{i=0}^k \sum_{j=0}^i a_{ij} P_j(x) P_{i-j}(y)\,.

    A flattened parameter list looks like this:

    .. math:: [a_{00}, a_{01}, a_{02}, ...,
               a_{0k}, a_{10}, a_{11}, ... a_{k0}]\,.
    """

    def __init__(self, degree):
        self.degree = degree
        self.paramcount = int(binom(degree+2, 2))
        self.params = None

    def set_params(self, params):
        """Set flat parameters list, e.g. after optimization.

        Parameters
        ----------
        params : list or array of length m
        """
        self.params = params

    def get_params(self):
        """Returns internal parameter list or `None` if not set.

        Returns
        -------
        list or array or None
            Parameter list of length m.
        """
        return self.params

    def model(self, points):
        """Evaluate model with internal parameters.

        Parameters
        ----------
        points : (N,2)-shaped array
            List of x,y points where to evaluate polynomial.

        Returns
        -------
        (N,)-shaped array
            1D array of values.

        Raises
        ------
        ValueError
            If no parameters where set.  Use `set_params`.
        """
        if self.params is None:
            raise ValueError("No internal parameters set.")
        return self.full_model(points, *self.params)

    @abstractmethod
    def full_model(self, points, *params):
        """Evaluate model with given parameters.

        This is an abstract method and has to be overwritten by subclasses.
        `model` uses this method with internal parameter list.

        Parameters
        ----------
        points : (N,2)-shaped array
            List of x,y points where to evaluate polynomial.
        params : (m,)-shaped array
            Flattened list of polynomial parameters.

        Returns
        -------
        (N,)-shaped array
            1D array of values
        """
        NotImplemented

    @abstractmethod
    def vandermonde(self, points):
        r"""Calculate Vandermonde matrix for given points.

        This is an abstract method.

        Parameters
        ----------
        points : (N,2)-shaped array
            List of x,y points where to evaluate polynomial.

        Returns
        -------
        (N,m)-shaped array
            Vandermonde matrix.
            N is the number of points, m the number of parameters.

        Notes
        -----
        The Vandermonde matrix for the points :math:`(x_j, y_j)` is

        .. math:: V = \begin{bmatrix}
              P_0(x_1) P_0(y_1) & \cdots & P_k(x_1) P_0(y_1) \\
              \vdots                & {}     & \vdots \\
              P_0(x_N) P_0(y_N) & \cdots & P_k(x_N) P_0(y_N)
            \end{bmatrix}

        and evaluation with a given parameters is then a simple
        multiplication

        .. math:: \begin{bmatrix}
              P(x_1, y_1) \\ \vdots \\ P(x_N, y_N)
            \end{bmatrix}
            = V \cdot \begin{bmatrix}
              a_0 \\ \vdots \\ a_m
            \end{bmatrix}\,.
        """
        NotImplemented

    def copy(self):
        """Copy of this polynomial with same parameter list.

        Returns
        -------
        Polynomial
            Copy
        """
        copy = self.__class__(self.degree)
        copy.set_params(self.params)
        return copy

    @abstractmethod
    def make_single_degree_subpoly(self, degree):
        """Make new polynomial containing terms of one degree only.

        This is an abstract method.

        Parameters
        ----------
        degree : int

        Returns
        -------
        New instance of same class with some parameters zeroed out.

        Raises
        ------
        ValueError
            If no internal parameters where set.
        """
        NotImplemented


class PolyVectorField:
    """Polynomial vector field mapping from 2D to 2D.

    Two polynomials make up the two components of the vector field.

    Parameters
    ----------
    xpoly : Polynomial
    ypoly : Polynomial
        Optional.  If left out or set to None,
        a copy of `xpoly` is used.

    Attributes
    ----------
    xpoly : Polynomial
    ypoly : Polynomial
    paramcount : int
         Total number of parameters.
    """

    def __init__(self, xpoly, ypoly=None):
        if ypoly is None:
            ypoly = xpoly.copy()
        self.xpoly, self.ypoly = xpoly, ypoly
        self.paramcount = xpoly.paramcount + ypoly.paramcount

    def get_degree(self):
        """Return degree of vector field.

        Raises
        ------
        ValueError
            If the degree is not the same for both polynomials.
        """
        if self.xpoly.degree != self.ypoly.degree:
            raise ValueError("Different degrees per component.")
        return self.xpoly.degree

    def set_params(self, params):
        """Set parameter list for both polynomials.

        Parameters
        ----------
        params : list
            Concatenated list of parameters for polynomial in x-direction
            and then parameters for the y-direction.
        """
        self.xpoly.set_params(params[:self.xpoly.paramcount])
        self.ypoly.set_params(params[self.xpoly.paramcount:])

    def get_params(self):
        """Get internal parameter list.

        Returns
        -------
        list
            Parameter list.  Returns `None` if one or both of the
            polynomials has no internal parameters.
        """
        xp = self.xpoly.get_params()
        yp = self.ypoly.get_params()
        if xp is None or yp is None:
            return None
        else:
            return list(xp) + list(yp)

    def model(self, points):
        """Evaluate model with internal parameters.

        Parameters
        ----------
        points : (N,2)-shaped array
            x,y points where to evaluate field.

        Returns
        -------
        (N,2)-shaped array
            Values at given points.
        """
        assert len(points.shape) == 2 and points.shape[1] == 2
        res = np.ndarray(shape=points.shape, dtype=float)
        res[:, 0] = self.xpoly.model(points)
        res[:, 1] = self.ypoly.model(points)
        return res

    def full_model(self, points, *params):
        """Evaluate model with given parameters.

        Parameters
        ----------
        points : (N,2)-shaped array
            x,y points where to evaluate field.
        params : list
            Concatenated list of parameters for both polynomials.

        Returns
        -------
        (N,2)-shaped array
            Values at given points.
        """
        assert points.shape[1] == 2
        res = np.ndarray(shape=points.shape, dtype=float)
        res[:, 0] = self.xpoly.full_model(
                        points, *params[:self.xpoly.paramcount])
        res[:, 1] = self.ypoly.full_model(
                        points, *params[self.xpoly.paramcount:])
        return res

    def vandermonde(self, points):
        r"""Calculate Vandermonde matrix.

        Parameters
        ----------
        points : (N,2)-shaped array

        Returns
        -------
        (2N,M)-shaped array
            Vandermonde matrix.  M is the total number of
            parameters for both polynomials.

        Notes
        -----
        For a vector field of a x polynomial with Vandermonde matrix
        :math:`V_x` and a y polynomial with Vandermonde matrix
        :math:`V_y` the combined matrix is

        .. math:: V = \begin{bmatrix}
              V_x        & \mathbf{0} \\
              \mathbf{0} & V_y \\
            \end{bmatrix}\,.
        """
        assert points.shape[1] == 2
        n = points.shape[0]
        pc1, pc2 = self.xpoly.paramcount, self.ypoly.paramcount
        V = np.zeros((2*n, pc1+pc2))
        V[:n, :pc1] = self.xpoly.vandermonde(points)
        V[n:, pc1:] = self.ypoly.vandermonde(points)
        return V

    def copy(self):
        """Copy vector field and polynomials with same parameter lists.

        Returns
        -------
        :class:`mascado.distortions.polynomials.PolyVectorField`
            New vector field.
        """
        return PolyVectorField(self.xpoly.copy(),
                               self.ypoly.copy())

    def make_single_degree_subpoly(self, degree):
        """Make new vector field containing terms of a single degree only.

        Parameters
        ----------
        degree : int

        Returns
        -------
        PolyVectorField
             New vector field of given degree.
        """
        subx = self.xpoly.make_single_degree_subpoly(degree)
        suby = self.ypoly.make_single_degree_subpoly(degree)
        return PolyVectorField(subx, suby)


class CoeffPolynomial (Polynomial):
    r"""Abstract class for 2D Polynomials representated as coefficient matrix.

    Parameters
    ----------
    degree : int
        Maximum degree :math:`k`.

    Notes
    -----
    A coefficient matrix is a square matrix where rows denote the
    degree of the x polynomial and columns the degree of the y
    polynomial.  In the simplest case (Cartesian polynomials)
    :math:`1+3x^2+0.5xy` is represented as

    .. math:: \begin{bmatrix}
          1 & 0 & 0 \\
          0 & 0.5 & 0 \\
          3 & 0 & 0
        \end{bmatrix}

    The coefficient matrix is an upper left triangular matrix, because
    ther lower right triangle represents terms of higher degree.

    """

    def __init__(self, degree):
        super().__init__(degree)

    @staticmethod
    def upper_left_triangular_mask(degree):
        """Make mask that covers the upper left trangle of a matrix.

        Parameters
        ----------
        degree : int

        Returns
        -------
        (degree+1, degree+1)-shaped bool array
        """
        return np.fliplr(np.triu(np.full((degree+1, degree+1), True)))

    @staticmethod
    def upper_left_triangular_indices(degree):
        """
        Returns indices you can use to convert a flattened array
        back to an upper left triangular matrix.

        Parameters
        ----------
        degree : int

        Returns
        -------
        2-tuple of (m,)-shaped arrays
            Indices to upper left triangle of a (degree+1, degree+1)-matrix.

        Examples
        --------
        Convert a flattened parameter list to a coefficient matrix:
        >>> degree = 1
        >>> paramlist = [1, -3, 0.12e-3]
        >>> idxs = CoeffPolynomial.upper_left_triangular_indices(degree)
        >>> mat = np.zeros((degree+1, degree+1))
        >>> mat[idxs] = paramlist

        """
        mask = CoeffPolynomial.upper_left_triangular_mask(degree)
        return np.unravel_index(np.where(np.ravel(mask)), (degree+1, degree+1))

    def params_to_coeff(self, params):
        """Convert parameter list to coefficient matrix."""
        ulidxs = CoeffPolynomial.upper_left_triangular_indices(self.degree)
        coeff = np.zeros((self.degree+1, self.degree+1))
        coeff[ulidxs] = params
        return coeff

    def coeff_to_params(self, coeff):
        """Convert coefficient matrix to parameter list."""
        ulmask = CoeffPolynomial.upper_left_triangular_mask(self.degree)
        return coeff[ulmask]

    def make_single_degree_subpoly(self, degree):
        """Concrete implementation."""
        if self.params is None:
            raise ValueError("No internal parameters set.")
        coeff = self.params_to_coeff(self.params)
        subcoeff = np.zeros((degree+1, degree+1))
        for m in range(degree+1):
            subcoeff[m, degree-m] = coeff[m, degree-m]
        subpoly = self.__class__(degree)
        subpoly.set_params(subpoly.coeff_to_params(subcoeff))
        return subpoly


class Cartesian (CoeffPolynomial):
    """Concrete implementation of :class:`Polynomial`."""

    def full_model(self, points, *params):
        """"""
        return polyval2d(points[:, 0], points[:, 1],
                         self.params_to_coeff(params))

    def vandermonde(self, points):
        """"""
        vander = polyvander2d(points[:, 0], points[:, 1],
                              [self.degree, self.degree])
        ulmask = CoeffPolynomial.upper_left_triangular_mask(self.degree)
        return vander[:, ulmask.ravel()]


class Legendre (CoeffPolynomial):
    """Concrete implementation of :class:`Polynomial`."""

    def full_model(self, points, *params):
        """"""
        return legval2d(points[:, 0], points[:, 1],
                        self.params_to_coeff(params))

    def vandermonde(self, points):
        """"""
        vander = legvander2d(points[:, 0], points[:, 1],
                             [self.degree, self.degree])
        ulmask = CoeffPolynomial.upper_left_triangular_mask(self.degree)
        return vander[:, ulmask.ravel()]


class Chebyshev (CoeffPolynomial):
    """Concrete implementation of :class:`Polynomial`."""

    def full_model(self, points, *params):
        """"""
        return chebval2d(points[:, 0], points[:, 1],
                         self.params_to_coeff(params))

    def vandermonde(self, points):
        """"""
        vander = chebvander2d(points[:, 0], points[:, 1],
                              [self.degree, self.degree])
        ulmask = CoeffPolynomial.upper_left_triangular_mask(self.degree)
        return vander[:, ulmask.ravel()]


class Zernike (CoeffPolynomial):
    """Zernike polynomials.

    Concrete implementation of :class:`Polynomial`.

    Implementation details from *Hedser van Brug (1997):*
    Efficient Cartesian representation of Zernike polynomials
    in computer memory. *Proc. SPIE 3190*

    Parameters
    ----------
    degree : int

    Attributes
    ----------
    zernikes : list of arrays
        Precalculated coefficient matrices for single Zernike
        polynomials. Arranged in the same order as
        flattened parameter list.

    Notes
    -----
    Zernike polynomials are defined on the unit circle.
    """

    def __init__(self, degree):
        super().__init__(degree)
        self.zernikes = [
            self.zernikepoly(i, j)
            for j in range(degree+1)
            for i in range(j, degree+1)]
        assert len(self.zernikes) == self.paramcount

    def zernikepoly(self, n, m):
        r"""Build Cartesian coefficient matrix for :math:`Z_n^m`.

        :math:`m,n>=0` and :math:`m<=n`.

        Returns
        -------
        function(x, y) -> 1D float array
            Function that takes x and y coordinates and calculates result for
            given :math:`Z_n^m`.
            (Wrapped `numpy.polynomials.polynomial.polyval2d()`)
        """
        coeff = np.zeros(shape=(n+1, n+1))
        h = n - 2*m
        if h <= 0:
            p = 0
            q = -h//2 if n % 2 == 0 else (-h-1)//2
        else:
            p = 1
            q = h//2-1 if n % 2 == 0 else (h-1)//2
        h = abs(h)
        m = (n - h) // 2
        for i in range(q+1):
            for j in range(m+1):
                for k in range(m-j+1):
                    factor = 1 if (i+j) % 2 == 0 else -1
                    factor *= binom(h, 2*i+p)
                    factor *= binom(m-j, k)
                    factor *= factorial(n-j) / (factorial(j) * factorial(m-j)
                                                * factorial(n-m-j))
                    ypow = 2 * (i + k) + p
                    xpow = n - 2 * (i + j + k) - p
                    coeff[xpow, ypow] += factor
        return lambda x, y: polyval2d(x, y, coeff)

    def full_model(self, points, *params):
        """"""
        result = np.zeros(shape=(points.shape[0],))
        for i, factor in enumerate(params):
            result += factor * self.zernikes[i](points[:, 0], points[:, 1])
        return result

    def vandermonde(self, points):
        """"""
        vander = np.zeros(shape=(points.shape[0], self.paramcount))
        for i in range(self.paramcount):
            vander[:, i] = self.zernikes[i](points[:, 0], points[:, 1])
        return vander
