
from unittest import TestCase
import unittest
import pygem.ffd as ffd
import pygem.ffd_parameters as ffdp
import pygem.file_handler as fh
import numpy as np



class TestFFD(TestCase):


	def test_ffd_identity(self):
		params = ffdp.FFDParameters()
		params.read_parameters_file(filename='tests/test_datasets/parameters_test_ffd_identity.prm')
		mesh_points = np.load('tests/test_datasets/meshpoints_sphere_orig.npy')
		free_form = ffd.FFD(params, mesh_points)
		mesh_points_test = free_form.perform()

		np.testing.assert_array_almost_equal(mesh_points_test, mesh_points)
		

	def test_ffd_sphere_mod(self):
		params = ffdp.FFDParameters()
		params.read_parameters_file(filename='tests/test_datasets/parameters_test_ffd_sphere.prm')
		mesh_points = np.load('tests/test_datasets/meshpoints_sphere_orig.npy')
		mesh_points_ref = np.load('tests/test_datasets/meshpoints_sphere_mod.npy')
		free_form = ffd.FFD(params, mesh_points)
		mesh_points_test = free_form.perform()

		np.testing.assert_array_almost_equal(mesh_points_test, mesh_points_ref)
		

