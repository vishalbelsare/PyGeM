
from unittest import TestCase
import unittest
import pygem.ffd_parameters as ffdp
import numpy as np


class TestFFDParameters(TestCase):


	def test_class_members_default_0(self):
		params = ffdp.FFDParameters()
		assert params.Nx == 1
		assert params.Ny == 1
		assert params.Nz == 1


	def test_class_members_default_1(self):
		params = ffdp.FFDParameters()
		assert params.conversion_unit == 1.
		assert params.a == 1.
		assert params.b == 1.
		assert params.c == 1.
		assert params.x0 == 0.
		assert params.y0 == 0.
		assert params.z0 == 0.
		assert params.aX == 0
		assert params.aY == 0
		assert params.aZ == 0


	def test_class_members_default_2(self):
		params = ffdp.FFDParameters()
		np.testing.assert_array_almost_equal(params.muXX, np.zeros((2, 2, 2)))
		np.testing.assert_array_almost_equal(params.muYY, np.zeros((2, 2, 2)))
		np.testing.assert_array_almost_equal(params.muZZ, np.zeros((2, 2, 2)))
		np.testing.assert_array_almost_equal(params.psi, np.diag([1, 1, 1]))
		np.testing.assert_array_almost_equal(params.inv_psi, np.diag([1, 1, 1]))
		np.testing.assert_array_almost_equal(params.rot_mat, np.eye(3))
		np.testing.assert_array_almost_equal(params.P0, np.zeros(3))
		np.testing.assert_array_almost_equal(params.P1, np.zeros(3))
		np.testing.assert_array_almost_equal(params.P2, np.zeros(3))
		np.testing.assert_array_almost_equal(params.P3, np.zeros(3))


	def test_class_members_0(self):
		params = ffdp.FFDParameters([2, 3, 5])
		assert params.Nx == 2
		assert params.Ny == 3
		assert params.Nz == 5
		np.testing.assert_array_almost_equal(params.muXX, np.zeros((3, 4, 6)))
		np.testing.assert_array_almost_equal(params.muYY, np.zeros((3, 4, 6)))
		np.testing.assert_array_almost_equal(params.muZZ, np.zeros((3, 4, 6)))

