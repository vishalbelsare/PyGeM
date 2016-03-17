
from unittest import TestCase
import unittest
import pygem.ffd_parameters as ffdp
import numpy as np


class TestFFDParameters(TestCase):


	def test_class_members_default_0(self):
		params = ffdp.FFDParameters()
		assert params.n_control_points_x == 1
		assert params.n_control_points_y == 1
		assert params.n_control_points_z == 1


	def test_class_members_default_1(self):
		params = ffdp.FFDParameters()
		assert params.conversion_unit == 1.
		assert params.lenght_box_x == 1.
		assert params.lenght_box_y == 1.
		assert params.lenght_box_z == 1.
		assert params.origin_box_x == 0.
		assert params.origin_box_y == 0.
		assert params.origin_box_z == 0.
		assert params.rot_angle_x == 0
		assert params.rot_angle_y == 0
		assert params.rot_angle_z == 0


	def test_class_members_default_2(self):
		params = ffdp.FFDParameters()
		np.testing.assert_array_almost_equal(params.array_mu_x, np.zeros((2, 2, 2)))
		np.testing.assert_array_almost_equal(params.array_mu_y, np.zeros((2, 2, 2)))
		np.testing.assert_array_almost_equal(params.array_mu_z, np.zeros((2, 2, 2)))
		np.testing.assert_array_almost_equal(params.psi_mapping, np.diag([1, 1, 1]))
		np.testing.assert_array_almost_equal(params.inv_psi_mapping, np.diag([1, 1, 1]))
		np.testing.assert_array_almost_equal(params.rot_mat, np.eye(3))
		np.testing.assert_array_almost_equal(params.position_vertex_0, np.zeros(3))
		np.testing.assert_array_almost_equal(params.position_vertex_1, np.zeros(3))
		np.testing.assert_array_almost_equal(params.position_vertex_2, np.zeros(3))
		np.testing.assert_array_almost_equal(params.position_vertex_3, np.zeros(3))


	def test_class_members_0(self):
		params = ffdp.FFDParameters([2, 3, 5])
		assert params.n_control_points_x == 2
		assert params.n_control_points_y == 3
		assert params.n_control_points_z == 5
		np.testing.assert_array_almost_equal(params.array_mu_x, np.zeros((3, 4, 6)))
		np.testing.assert_array_almost_equal(params.array_mu_y, np.zeros((3, 4, 6)))
		np.testing.assert_array_almost_equal(params.array_mu_z, np.zeros((3, 4, 6)))

