
from unittest import TestCase
import unittest
import pygem.file_handler as fh
import numpy as np
import filecmp
import os


class TestFileHandler(TestCase):


	def test_base_class_filename(self):
		file_handler = fh.FileHandler()
		assert file_handler.filename == None


	def test_base_class_extension(self):
		file_handler = fh.FileHandler()
		assert file_handler.extension == None

	
	def test_base_class_parse(self):
		file_handler = fh.FileHandler()
		with self.assertRaises(NotImplementedError):
			file_handler.parse()


	def test_base_class_write(self):
		file_handler = fh.FileHandler()
		mesh_points = np.zeros((3, 3))
		with self.assertRaises(NotImplementedError):
			file_handler.write(mesh_points, 'output')


	# UNV tests
	def test_unv_filename_member(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		assert unv_handler.filename == 'tests/test_datasets/test_square.unv'
	

	def test_unv_extension_member(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		assert unv_handler.extension == 'unv'


	def test_unv_failing_filename_type(self):
		with self.assertRaises(TypeError):
			unv_handler = fh.UnvHandler(3)


	def test_unv_instantiation(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
	

	def test_unv_failing_check_extension(self):
		with self.assertRaises(ValueError):
			unv_handler = fh.UnvHandler('tests/test_datasets/test_square.iges')


	def test_unv_parse_shape(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		assert mesh_points.shape == (256, 3)


	def test_unv_parse_coords_1(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		np.testing.assert_almost_equal(mesh_points[33][0], 1.0)


	def test_unv_parse_coords_2(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		np.testing.assert_almost_equal(mesh_points[178][1], 0.4)


	def test_unv_parse_coords_3(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		np.testing.assert_almost_equal(mesh_points[100][2], 0.0)


	def test_unv_parse_coords_4(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		np.testing.assert_almost_equal(mesh_points[0][0], 0.0)
		
		
	def test_unv_parse_coords_5(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		np.testing.assert_almost_equal(mesh_points[-1][2], 0.0)		


	def test_unv_write_failing_outfile_type(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		with self.assertRaises(TypeError):
			unv_handler.write(mesh_points, 3)
 

	def test_unv_write_outfile(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		mesh_points = unv_handler.parse()
		mesh_points[0][0] = 2.2
		mesh_points[5][1] = 4.3
		mesh_points[9][2] = 0.5
		mesh_points[45][0] = 7.2
		mesh_points[132][1] = -1.2
		mesh_points[255][2] = -3.6

		outfilename = 'tests/test_datasets/test_square_out.unv'
		outfilename_expected = 'tests/test_datasets/test_square_out_true.unv'

		unv_handler.write(mesh_points, outfilename)
		self.assertTrue(filecmp.cmp(outfilename, outfilename_expected))
		os.remove(outfilename)


	# STL tests
	def test_stl_filename_member(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		assert stl_handler.filename == 'tests/test_datasets/test_sphere.stl'
	

	def test_stl_extension_member(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		assert stl_handler.extension == 'stl'


	def test_stl_failing_filename_type(self):
		with self.assertRaises(TypeError):
			stl_handler = fh.StlHandler(3)


	def test_stl_instantiation(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
	

	def test_stl_failing_check_extension(self):
		with self.assertRaises(ValueError):
			stl_handler = fh.StlHandler('tests/test_datasets/test_square.iges')


	def test_stl_parse_shape(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		assert mesh_points.shape == (7200, 3)


	def test_stl_parse_coords_1(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		np.testing.assert_almost_equal(mesh_points[33][0], -21.31975937)


	def test_stl_parse_coords_2(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		np.testing.assert_almost_equal(mesh_points[1708][1], 2.58431911)


	def test_stl_parse_coords_3(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		np.testing.assert_almost_equal(mesh_points[3527][2], -2.47207999)


	def test_stl_parse_coords_4(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		np.testing.assert_almost_equal(mesh_points[0][0], -21.31975937)	


	def test_stl_parse_coords_5(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		np.testing.assert_almost_equal(mesh_points[-1][2], -39.05963898)		


	def test_stl_write_failing_outfile_type(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		with self.assertRaises(TypeError):
			stl_handler.write(mesh_points, 3)
 

	def test_stl_write_outfile(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()
		mesh_points[0] = [-40.2, -20.5, 60.9]
		mesh_points[1] = [-40.2, -10.5, 60.9]
		mesh_points[2] = [-40.2, -10.5, 60.9]
		mesh_points[2000] = [-40.2, -20.5, 60.9]
		mesh_points[2001] = [-40.2, -10.5, 60.9]
		mesh_points[2002] = [-40.2, -10.5, 60.9]
		mesh_points[6100] = [-40.2, -20.5, 60.9]
		mesh_points[6101] = [-40.2, -10.5, 60.9]
		mesh_points[6102] = [-40.2, -10.5, 60.9]

		outfilename = 'tests/test_datasets/test_sphere_out.stl'
		outfilename_expected = 'tests/test_datasets/test_sphere_out_true.stl'

		stl_handler.write(mesh_points, outfilename)
		self.assertTrue(filecmp.cmp(outfilename, outfilename_expected))
		os.remove(outfilename)


	def test_stl_plot_failing_outfile_type(self):
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		with self.assertRaises(TypeError):
			stl_handler.plot(plot_file=3)
			
	
	# openFOAM tests
	def test_open_foam_filename_member(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		assert open_foam_handler.filename == 'tests/test_datasets/test_openFOAM'


	def test_open_foam_failing_filename_type(self):
		with self.assertRaises(TypeError):
			open_foam_handler = fh.OpenFoamHandler(3)


	def test_open_foam_instantiation(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')


	def test_open_foam_parse_shape(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		assert mesh_points.shape == (21812, 3)


	def test_open_foam_parse_coords_1(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		np.testing.assert_almost_equal(mesh_points[33][0], 1.42254)


	def test_open_foam_parse_coords_2(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		np.testing.assert_almost_equal(mesh_points[1708][1], -3.13059)


	def test_open_foam_parse_coords_3(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		np.testing.assert_almost_equal(mesh_points[3527][2], .0)


	def test_open_foam_parse_coords_4(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		np.testing.assert_almost_equal(mesh_points[0][0], -17.5492)	


	def test_open_foam_parse_coords_5(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		np.testing.assert_almost_equal(mesh_points[-1][2], 0.05)		


	def test_open_foam_write_failing_outfile_type(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		with self.assertRaises(TypeError):
			open_foam_handler.write(mesh_points, 3)
 

	def test_open_foam_write_outfile(self):
		open_foam_handler = fh.OpenFoamHandler('tests/test_datasets/test_openFOAM')
		mesh_points = open_foam_handler.parse()
		mesh_points[0] = [-14.,  1.55, 0.2]
		mesh_points[1] = [-14.3, 2.55, 0.3]
		mesh_points[2] = [-14.3, 2.55, 0.3]
		mesh_points[2000] = [7.8, -42.8, .0]
		mesh_points[2001] = [8.8, -41.8, .1]
		mesh_points[2002] = [9.8, -40.8, .0]
		mesh_points[-3] = [236.3, 183.7, 0.06]
		mesh_points[-2] = [237.3, 183.7, 0.06]
		mesh_points[-1] = [236.3, 185.7, 0.06]

		outfilename = 'tests/test_datasets/test_openFOAM_out'
		outfilename_expected = 'tests/test_datasets/test_openFOAM_out_true'

		open_foam_handler.write(mesh_points, outfilename)
		self.assertTrue(filecmp.cmp(outfilename, outfilename_expected))
		os.remove(outfilename)


