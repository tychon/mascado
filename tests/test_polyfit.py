
# Copyright 2018 Hannes Riechert at Max-Planck-Institute for Astronomy.
# Licensed under GPL-3.0-or-later.  See COPYING for details.


from unittest import TestCase

import numpy as np
from mascado.distortions.polynomials import PolyVectorField, Legendre
from mascado.distortions.polyfit import (
    polyfit_svd, polyfit_svd_iterative)


class PolyFitTestCase (TestCase):
    def test_polyfit_randomized(self):
        # prepare random field
        vf = PolyVectorField(Legendre(3))
        vf.set_params(np.random.randn(vf.paramcount))
        xx, yy = np.meshgrid(np.linspace(-1, 1, 10), np.linspace(-1, 1, 10))
        # make artificial data
        grid = np.stack([xx.ravel(), yy.ravel()], axis=1)
        sigmas = np.absolute(np.random.normal(0.1, 0.01, size=(grid.shape[0])))
        data = np.random.normal(vf.model(grid),
                                np.stack([sigmas, sigmas], axis=1))
        # fit field
        params, residuals, resvar = polyfit_svd(vf, grid, data, sigmas,
                                                info=False)
        self.assertTrue(0.6 < resvar < 1.4)
        self.assertTrue(np.allclose(params, vf.get_params(), atol=0.1))

    def test_outlier_rejection_randomized(self):
        # prepare random field
        vf = PolyVectorField(Legendre(3))
        vf.set_params(np.random.randn(vf.paramcount))
        xx, yy = np.meshgrid(np.linspace(-1, 1, 10), np.linspace(-1, 1, 10))
        # make artificial data
        grid = np.stack([xx.ravel(), yy.ravel()], axis=1)
        sigmas = np.absolute(np.random.normal(0.01, 0.001,
                                              size=(grid.shape[0])))
        vectors = np.random.normal(vf.model(grid),
                                   np.stack([sigmas, sigmas], axis=1))
        vectors[15] += sigmas[5] * 10
        vectors[55] -= sigmas[5] * 10
        params, inliers, residuals, resvar = polyfit_svd_iterative(
            vf, grid, vectors, sigmas=sigmas,
            maxoutliers=2, info=False)
        self.assertEqual(np.count_nonzero(~inliers), 2)
        self.assertTrue(0.6 < resvar < 1.4)

    def test_outlier_rejection_maxoutliers_bound(self):
        # prepare random field
        vf = PolyVectorField(Legendre(3))
        vf.set_params(np.random.randn(vf.paramcount))
        xx, yy = np.meshgrid(np.linspace(-1, 1, 10), np.linspace(-1, 1, 10))
        # make artificial data
        grid = np.stack([xx.ravel(), yy.ravel()], axis=1)
        sigmas = np.absolute(np.random.normal(0.01, 0.001,
                                              size=(grid.shape[0])))
        vectors = np.random.normal(vf.model(grid),
                                   np.stack([sigmas, sigmas], axis=1))
        vectors[15] += sigmas[5] * 50
        vectors[55] -= sigmas[5] * 50
        with self.assertRaises(RuntimeError):
            params, inliers, residuals, resvar = polyfit_svd_iterative(
                vf, grid, vectors, sigmas=sigmas,
                maxoutliers=0, info=False)
