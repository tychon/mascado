
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.


from unittest import TestCase

import numpy as np

from mascado.distortions.polynomials import (
    Cartesian, Legendre, Chebyshev, Zernike,
    PolyVectorField)


class PolynomialTestCase (TestCase):
    """
    Results are not externally validated.
    Only change, not correctness is tested.
    """

    def test_cartesian(self):
        poly = Cartesian(3)
        poly.set_params([
            1, 0, 1, 0,
            0, 0, 0,
            0, 0,
            1])
        grid = np.array([[0, 0], [1, 1], [0, -1]])
        values = poly.model(grid)
        self.assertTrue(all(values == np.array([1, 3, 2])))

    def test_cartesian_vandermonde(self):
        poly = Cartesian(3)
        vander = poly.vandermonde(np.array([[0, 0], [0, -1]]))
        self.assertTrue(np.all(vander == np.array(
            [[1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
             [1., -1.,  1., -1.,  0., -0.,  0.,  0., -0.,  0.]]
            )))

    def test_legendre(self):
        poly = Legendre(3)
        poly.set_params([
            1, 0, 1, 0,
            0, 0, 0,
            0, 0,
            1])
        grid = np.array([[0, 0], [1, 1], [0, -1]])
        values = poly.model(grid)
        self.assertTrue(all(values == np.array([0.5, 3, 2])))

    def test_legendre_vandermonde(self):
        poly = Legendre(3)
        vander = poly.vandermonde(np.array([[0, 0], [0, -1]]))
        self.assertTrue(np.all(vander == np.array(
            [[1.,  0., -0.5, -0.,  0.,  0., -0., -0.5,  -0., -0.],
             [1., -1.,   1., -1.,  0., -0.,  0., -0.5,  0.5, -0.]]
            )))

    def test_chebyshev(self):
        poly = Chebyshev(3)
        poly.set_params([
            1, 0, 1, 0,
            0, 0, 0,
            0, 0,
            1])
        grid = np.array([[0, 0], [1, 1], [0, -1]])
        values = poly.model(grid)
        self.assertTrue(all(values == np.array([0, 3, 2])))

    def test_chebyshev_vandermonde(self):
        poly = Chebyshev(3)
        vander = poly.vandermonde(np.array([[0, 0], [0, -1]]))
        self.assertTrue(np.all(vander == np.array(
            [[1.,  0., -1., -0.,  0.,  0., -0., -1., -0., -0.],
             [1., -1.,  1., -1.,  0., -0.,  0., -1.,  1., -0.]]
            )))

    def test_zernike(self):
        poly = Zernike(3)
        poly.set_params([
            1, 0, 1, 0,
            0, 0, 0,
            0, 0,
            1])
        grid = np.array([[0, 0], [1, 1], [0, -1]])
        values = poly.model(grid)
        self.assertTrue(all(values == np.array([1, 1, 1])))

    def test_zernike_vandermonde(self):
        poly = Zernike(3)
        vander = poly.vandermonde(np.array([[0, 0], [0, -1]]))
        self.assertTrue(np.all(vander == np.array(
            [[1.,  0.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  0.],
             [1., -1.,  0.,  1.,  0.,  1., -1., -1.,  0.,  0.]]
            )))

    def test_copy(self):
        vf = PolyVectorField(Legendre(3))
        vf.set_params(np.random.randn(20))
        vf2 = vf.copy()
        self.assertTrue(vf.xpoly != vf2.xpoly)
        self.assertTrue(vf.ypoly != vf2.ypoly)
        self.assertTrue(all(vf.xpoly.get_params() == vf2.xpoly.get_params()))
        self.assertTrue(all(vf.ypoly.get_params() == vf2.ypoly.get_params()))
