
from unittest import TestCase
import unittest
import pygem.file_handler as fh
import numpy as np
import filecmp
import os


class TestFileHandler(TestCase):


	def test_base_class_infile(self):
		file_handler = fh.FileHandler()
		assert file_handler.infile == None


	def test_base_class_outfile(self):
		file_handler = fh.FileHandler()
		assert file_handler.outfile == None


	def test_base_class_extension(self):
		file_handler = fh.FileHandler()
		assert file_handler.extension == None

	
	def test_base_class_parse(self):
		file_handler = fh.FileHandler()
		with self.assertRaises(NotImplementedError):
			file_handler.parse('input')


	def test_base_class_write(self):
		file_handler = fh.FileHandler()
		mesh_points = np.zeros((3, 3))
		with self.assertRaises(NotImplementedError):
			file_handler.write(mesh_points, 'output')


	# UNV tests
	def test_unv_instantiation(self):
		unv_handler = fh.UnvHandler()
	

	def test_unv_default_infile_member(self):
		unv_handler = fh.UnvHandler()
		assert unv_handler.infile == None


	def test_unv_default_outfile_member(self):
		unv_handler = fh.UnvHandler()
		assert unv_handler.outfile == None


	def test_unv_default_extension_member(self):
		unv_handler = fh.UnvHandler()
		assert unv_handler.extension == 'unv'
	

	def test_unv_parse_failing_filename_type(self):
		unv_handler = fh.UnvHandler()
		with self.assertRaises(TypeError):
			mesh_points = unv_handler.parse(5.2)

	
	def test_unv_parse_failing_check_extension(self):
		unv_handler = fh.UnvHandler()
		with self.assertRaises(ValueError):
			mesh_points = unv_handler.parse('tests/test_datasets/test_square.iges')


	def test_unv_parse_infile(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		assert unv_handler.infile == 'tests/test_datasets/test_square.unv'


	def test_unv_parse_shape(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		assert mesh_points.shape == (256, 3)


	def test_unv_parse_coords_1(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		np.testing.assert_almost_equal(mesh_points[33][0], 1.0)


	def test_unv_parse_coords_2(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		np.testing.assert_almost_equal(mesh_points[178][1], 0.4)


	def test_unv_parse_coords_3(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		np.testing.assert_almost_equal(mesh_points[100][2], 0.0)


	def test_unv_parse_coords_4(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		np.testing.assert_almost_equal(mesh_points[0][0], 0.0)
		
		
	def test_unv_parse_coords_5(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		np.testing.assert_almost_equal(mesh_points[-1][2], 0.0)		


	def test_unv_write_failing_filename_type(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		with self.assertRaises(TypeError):
			unv_handler.write(mesh_points, -2)


	def test_unv_write_failing_check_extension(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		with self.assertRaises(ValueError):
			unv_handler.write(mesh_points, 'tests/test_datasets/test_square.iges')


	def test_unv_write_failing_infile_instantiation(self):
		unv_handler = fh.UnvHandler()
		mesh_points = np.zeros((20, 3))
		with self.assertRaises(RuntimeError):
 			unv_handler.write(mesh_points, 'tests/test_datasets/test_square_out.unv')


	def test_unv_write_outfile(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
		outfilename = 'tests/test_datasets/test_square_out.unv'
		unv_handler.write(mesh_points, outfilename)
		assert unv_handler.outfile == outfilename
		os.remove(outfilename)


	def test_unv_write_comparison(self):
		unv_handler = fh.UnvHandler()
		mesh_points = unv_handler.parse('tests/test_datasets/test_square.unv')
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
		
		
	# VTK tests
	def test_vtk_instantiation(self):
		vtk_handler = fh.VtkHandler()
	

	def test_vtk_default_infile_member(self):
		vtk_handler = fh.VtkHandler()
		assert vtk_handler.infile == None


	def test_vtk_default_outfile_member(self):
		vtk_handler = fh.VtkHandler()
		assert vtk_handler.outfile == None


	def test_vtk_default_extension_member(self):
		vtk_handler = fh.VtkHandler()
		assert vtk_handler.extension == 'vtk'
	

	def test_vtk_parse_failing_filename_type(self):
		vtk_handler = fh.VtkHandler()
		with self.assertRaises(TypeError):
			mesh_points = vtk_handler.parse(5.2)

	
	def test_vtk_parse_failing_check_extension(self):
		vtk_handler = fh.VtkHandler()
		with self.assertRaises(ValueError):
			mesh_points = vtk_handler.parse('tests/test_datasets/test_square.iges')


	def test_vtk_parse_infile(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		assert vtk_handler.infile == 'tests/test_datasets/test_red_blood_cell.vtk'


	def test_vtk_parse_shape(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		assert mesh_points.shape == (500, 3)


	def test_vtk_parse_coords_1(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		np.testing.assert_almost_equal(mesh_points[33][0], -2.2977099)


	def test_vtk_parse_coords_2(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		np.testing.assert_almost_equal(mesh_points[178][1], 0.143506)


	def test_vtk_parse_coords_3(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		np.testing.assert_almost_equal(mesh_points[100][2], 2.3306999)


	def test_vtk_parse_coords_4(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		np.testing.assert_almost_equal(mesh_points[0][0], -3.42499995)
		
		
	def test_vtk_parse_coords_5(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		np.testing.assert_almost_equal(mesh_points[-1][2], -2.8480699)		


	def test_vtk_write_failing_filename_type(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		with self.assertRaises(TypeError):
			vtk_handler.write(mesh_points, -2)


	def test_vtk_write_failing_check_extension(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		with self.assertRaises(ValueError):
			vtk_handler.write(mesh_points, 'tests/test_datasets/test_square.iges')


	def test_vtk_write_failing_infile_instantiation(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = np.zeros((20, 3))
		with self.assertRaises(RuntimeError):
 			vtk_handler.write(mesh_points, 'tests/test_datasets/test_red_blood_cell_out.vtk')


	def test_vtk_write_outfile(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		outfilename = 'tests/test_datasets/test_red_blood_cell_out.vtk'
		vtk_handler.write(mesh_points, outfilename)
		assert vtk_handler.outfile == outfilename
		os.remove(outfilename)


	def test_vtk_write_comparison(self):
		vtk_handler = fh.VtkHandler()
		mesh_points = vtk_handler.parse('tests/test_datasets/test_red_blood_cell.vtk')
		mesh_points[0][0] = 2.2
		mesh_points[5][1] = 4.3
		mesh_points[9][2] = 0.5
		mesh_points[45][0] = 7.2
		mesh_points[132][1] = -1.2
		mesh_points[255][2] = -3.6

		outfilename = 'tests/test_datasets/test_red_blood_cell_out.vtk'
		outfilename_expected = 'tests/test_datasets/test_red_blood_cell_out_true.vtk'

		vtk_handler.write(mesh_points, outfilename)
		self.assertTrue(filecmp.cmp(outfilename, outfilename_expected))
		os.remove(outfilename)


	# STL tests
	def test_stl_instantiation(self):
		stl_handler = fh.StlHandler()
	

	def test_stl_default_infile_member(self):
		stl_handler = fh.StlHandler()
		assert stl_handler.infile == None


	def test_stl_default_outfile_member(self):
		stl_handler = fh.StlHandler()
		assert stl_handler.outfile == None


	def test_stl_default_extension_member(self):
		stl_handler = fh.StlHandler()
		assert stl_handler.extension == 'stl'
	

	def test_stl_parse_failing_filename_type(self):
		stl_handler = fh.StlHandler()
		with self.assertRaises(TypeError):
			mesh_points = stl_handler.parse(5.2)

	
	def test_stl_parse_failing_check_extension(self):
		stl_handler = fh.StlHandler()
		with self.assertRaises(ValueError):
			mesh_points = stl_handler.parse('tests/test_datasets/test_square.iges')


	def test_stl_parse_infile(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		assert stl_handler.infile == 'tests/test_datasets/test_sphere.stl'


	def test_stl_parse_shape(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		assert mesh_points.shape == (7200, 3)


	def test_stl_parse_coords_1(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		np.testing.assert_almost_equal(mesh_points[33][0], -21.31975937)


	def test_stl_parse_coords_2(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		np.testing.assert_almost_equal(mesh_points[1708][1], 2.58431911)


	def test_stl_parse_coords_3(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		np.testing.assert_almost_equal(mesh_points[3527][2], -2.47207999)


	def test_stl_parse_coords_4(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		np.testing.assert_almost_equal(mesh_points[0][0], -21.31975937)	


	def test_stl_parse_coords_5(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		np.testing.assert_almost_equal(mesh_points[-1][2], -39.05963898)		


	def test_stl_write_failing_filename_type(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		with self.assertRaises(TypeError):
			stl_handler.write(mesh_points, 4.)


	def test_stl_write_failing_check_extension(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		with self.assertRaises(ValueError):
			stl_handler.write(mesh_points, 'tests/test_datasets/test_square.iges')


	def test_stl_write_failing_infile_instantiation(self):
		stl_handler = fh.StlHandler()
		mesh_points = np.zeros((40, 3))
		with self.assertRaises(RuntimeError):
 			stl_handler.write(mesh_points, 'tests/test_datasets/test_sphere_out.stl')


	def test_stl_write_outfile(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		outfilename = 'tests/test_datasets/test_sphere_out.stl'
		stl_handler.write(mesh_points, outfilename)
		assert stl_handler.outfile == outfilename
		os.remove(outfilename)


	def test_stl_write_comparison(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
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


	def test_stl_plot_save_fig(self):
		stl_handler = fh.StlHandler()
		mesh_points = stl_handler.parse('tests/test_datasets/test_sphere.stl')
		stl_handler.plot(save_fig=True)
		if not os.path.isfile('tests/test_datasets/test_sphere.png'):
			assert False
		else:
			os.remove('tests/test_datasets/test_sphere.png')


	def test_stl_plot_failing_outfile_type(self):
		stl_handler = fh.StlHandler()
		with self.assertRaises(TypeError):
			stl_handler.plot(plot_file=3)
			
	
	# openFOAM tests
	def test_open_foam_instantiation(self):
		open_foam_handler = fh.OpenFoamHandler()
	

	def test_open_foam_default_infile_member(self):
		open_foam_handler = fh.OpenFoamHandler()
		assert open_foam_handler.infile == None


	def test_open_foam_default_outfile_member(self):
		open_foam_handler = fh.OpenFoamHandler()
		assert open_foam_handler.outfile == None


	def test_open_foam_default_extension_member(self):
		open_foam_handler = fh.OpenFoamHandler()
		assert open_foam_handler.extension == ''
	

	def test_open_foam_parse_failing_filename_type(self):
		open_foam_handler = fh.OpenFoamHandler()
		with self.assertRaises(TypeError):
			mesh_points = open_foam_handler.parse(.2)

	
	def test_open_foam_parse_failing_check_extension(self):
		open_foam_handler = fh.OpenFoamHandler()
		with self.assertRaises(ValueError):
			mesh_points = open_foam_handler.parse('tests/test_datasets/test_square.iges')


	def test_open_foam_parse_infile(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		assert open_foam_handler.infile == 'tests/test_datasets/test_openFOAM'


	def test_open_foam_parse_shape(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		assert mesh_points.shape == (21812, 3)


	def test_open_foam_parse_coords_1(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		np.testing.assert_almost_equal(mesh_points[33][0], 1.42254)


	def test_open_foam_parse_coords_2(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		np.testing.assert_almost_equal(mesh_points[1708][1], -3.13059)


	def test_open_foam_parse_coords_3(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		np.testing.assert_almost_equal(mesh_points[3527][2], .0)


	def test_open_foam_parse_coords_4(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		np.testing.assert_almost_equal(mesh_points[0][0], -17.5492)	


	def test_open_foam_parse_coords_5(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		np.testing.assert_almost_equal(mesh_points[-1][2], 0.05)		


	def test_open_foam_write_failing_filename_type(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		with self.assertRaises(TypeError):
			open_foam_handler.write(mesh_points, -1.)


	def test_open_foam_write_failing_check_extension(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		with self.assertRaises(ValueError):
			open_foam_handler.write(mesh_points, 'tests/test_datasets/test_square.iges')


	def test_open_foam_write_failing_infile_instantiation(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = np.zeros((40, 3))
		with self.assertRaises(RuntimeError):
 			open_foam_handler.write(mesh_points, 'tests/test_datasets/test_openFOAM_out')


	def test_open_foam_write_outfile(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
		outfilename = 'tests/test_datasets/test_openFOAM_out'
		open_foam_handler.write(mesh_points, outfilename)
		assert open_foam_handler.outfile == outfilename
		os.remove(outfilename)


	def test_open_foam_write_comparison(self):
		open_foam_handler = fh.OpenFoamHandler()
		mesh_points = open_foam_handler.parse('tests/test_datasets/test_openFOAM')
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


