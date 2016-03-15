
from unittest import TestCase
import unittest
import pygem.affine_trans as at
import numpy as np


class TestAffineTrans(TestCase):


	def test_angles2matrix_rot0(self):
		rotz = 0
		roty = 0
		rotx = 0
		
		mat_exact = np.eye(3)
		mat_test = at.angles2matrix(rotz, roty, rotx)

		np.testing.assert_array_almost_equal(mat_exact, mat_test)

	def test_angles2matrix_rot_x(self):
		rotz = 0
		roty = 0
		rotx = 50*np.pi/180
		
		mat_exact = np.array([1., 0., 0., 0., 0.64278761, -0.76604444, 0., 0.76604444, 0.64278761]).reshape((3,3))
		mat_test = at.angles2matrix(rotz, roty, rotx)

		np.testing.assert_array_almost_equal(mat_exact, mat_test)

	def test_angles2matrix_rot_y(self):
		rotz = 0
		roty = 23*np.pi/180
		rotx = 0
		
		mat_exact = np.array([0.92050485, 0., 0.39073113, 0., 1., 0., -0.39073113, 0., 0.92050485]).reshape((3,3))
		mat_test = at.angles2matrix(rotz, roty, rotx)

		np.testing.assert_array_almost_equal(mat_exact, mat_test)

	def test_angles2matrix_rot_z(self):
		rotz = -57*np.pi/180
		roty = 0
		rotx = 0
		
		mat_exact = np.array([0.54463904, 0.83867057, 0., -0.83867057, 0.54463904, 0., 0., 0., 1.]).reshape((3,3))
		mat_test = at.angles2matrix(rotz, roty, rotx)

		np.testing.assert_array_almost_equal(mat_exact, mat_test)

	def test_angles2matrix_rot_xyz(self):
		rotz = 10*np.pi/180
		roty = 20*np.pi/180
		rotx = 30*np.pi/180
		
		mat_exact = np.array([0.92541658, -0.16317591, 0.34202014, 0.31879578, 0.82317294, -0.46984631, -0.20487413, 0.54383814, 0.81379768]).reshape((3,3))
		mat_test = at.angles2matrix(rotz, roty, rotx)

		np.testing.assert_array_almost_equal(mat_exact, mat_test)


	
