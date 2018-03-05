
from unittest import TestCase

import numpy as np
import pandas as pd

from maskastrometry.utility import affine


class AffineTestCase (TestCase):
    def test_plain_transfrom(self):
        mask = np.array([[0, 0], [1, 0], [1, 1]])
        image = np.array([[100, 1], [1, 0], [3, 1]])
        trafo = affine.affine_lstsq(mask, image)
        self.assertTrue(np.allclose(affine.affine_trafo(mask, trafo), image))

    def test_catalog_transform(self):
        mask = pd.DataFrame(
            np.array([[0, 0], [1, 0], [1, 1]]),
            index=['A', 'B', 'C'], columns=['xref', 'yref'])
        image = pd.DataFrame(
            np.array([[100, 1], [1, 0], [3, 1]]),
            index=['A', 'B', 'C'], columns=['xref', 'yref'])
        trafo = affine.affine_lstsq_cat(mask, image, x='xref', y='yref')
        mapped = affine.affine_trafo_cat(mask, trafo, x='xref', y='yref')
        self.assertTrue(np.allclose(mapped, image))

    def test_nonmatching_catalog_transform(self):
        mask = pd.DataFrame(
            np.array([[0, 0], [1, 0], [1, 1]]),
            index=['A', 'B', 'X'], columns=['xref', 'yref'])
        image = pd.DataFrame(
            np.array([[100, 1], [1, 0], [3, 1]]),
            index=['A', 'B', 'Y'], columns=['xref', 'yref'])
        with self.assertRaises(ValueError):
            affine.affine_lstsq_cat(mask, image, x='xref', y='yref')

    def test_nonmatching_ignore_index_transform(self):
        mask = pd.DataFrame(
            np.array([[0, 0], [1, 0], [1, 1]]),
            index=['A', 'B', 'X'], columns=['xref', 'yref'])
        image = pd.DataFrame(
            np.array([[100, 1], [1, 0], [3, 1]]),
            index=['A', 'B', 'Y'], columns=['xref', 'yref'])
        trafo = affine.affine_lstsq_cat(mask, image, x='xref', y='yref',
                                        ignore_index=True)
        mapped = affine.affine_trafo_cat(mask, trafo, x='xref', y='yref')
        self.assertTrue(np.allclose(mapped, image))
