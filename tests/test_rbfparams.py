
from unittest import TestCase
import unittest
import pygem.params as rbfp
import numpy as np
import filecmp
import os


class TestFFDParameters(TestCase):


	def test_class_members_default_basis(self):
		params = rbfp.RBFParameters()
		assert params.basis == None


	def test_class_members_default_radius(self):
		params = rbfp.RBFParameters()
		assert params.radius == None


	def test_class_members_default_n_control_points(self):
		params = rbfp.RBFParameters()
		assert params.n_control_points == None


	def test_class_members_default_original_control_points(self):
		params = rbfp.RBFParameters()
		assert params.original_control_points == None


	def test_class_members_default_deformed_control_points(self):
		params = rbfp.RBFParameters()
		assert params.deformed_control_points == None


	def test_read_parameters_basis(self):
		params = rbfp.RBFParameters()
		params.read_parameters('tests/test_datasets/parameters_rbf_default.prm')
		assert params.basis == 'gaussian_spline'


	def test_read_parameters_radius(self):
		params = rbfp.RBFParameters()
		params.read_parameters('tests/test_datasets/parameters_rbf_default.prm')
		assert params.radius == 0.5


	def test_read_parameters_n_control_points(self):
		params = rbfp.RBFParameters()
		params.read_parameters('tests/test_datasets/parameters_rbf_default.prm')
		assert params.n_control_points == 8


	def test_read_parameters_original_control_points(self):
		params = rbfp.RBFParameters()
		params.read_parameters('tests/test_datasets/parameters_rbf_default.prm')
		original_control_points_exact = np.array([0., 0., 0., 0., 0., 1., 0., 1., 0., 1., 0., 0., \
				0., 1., 1., 1., 0., 1., 1., 1., 0., 1., 1., 1.]).reshape((8, 3))
		np.testing.assert_array_almost_equal(params.original_control_points, original_control_points_exact)


	def test_read_parameters_deformed_control_points(self):
		params = rbfp.RBFParameters()
		params.read_parameters('tests/test_datasets/parameters_rbf_default.prm')
		deformed_control_points_exact = np.array([0., 0., 0., 0., 0., 1., 0., 1., 0., 1., 0., 0., \
				0., 1., 1., 1., 0., 1., 1., 1., 0., 1., 1., 1.]).reshape((8, 3))
		np.testing.assert_array_almost_equal(params.deformed_control_points, deformed_control_points_exact)


	def test_read_parameters_failing_filename_type(self):
		params = rbfp.RBFParameters()
		with self.assertRaises(TypeError):
			params.read_parameters(3)


	def test_read_parameters_failing_number_deformed_control_points(self):
		params = rbfp.RBFParameters()
		with self.assertRaises(TypeError):
			params.read_parameters('tests/test_datasets/parameters_rbf_bugged_01.prm')


	def test_print_info(self):
		params = rbfp.RBFParameters()
		params.print_info()
