
from unittest import TestCase
import unittest
import pygem.ffd as ffd
import pygem.ffd_parameters as ffdp
import pygem.file_handler as fh
import numpy as np



class TestFFD(TestCase):


	def test_class_default(self):
		params = ffdp.FFDParameters()
		params.read_parameters_file(filename='tests/test_datasets/parameters_sphere.prm')
		stl_handler = fh.StlHandler('tests/test_datasets/test_sphere.stl')
		mesh_points = stl_handler.parse()

		free_form = ffd.FFD(params, mesh_points)
		new_mesh_points = free_form.perform()


