
from unittest import TestCase
import unittest
import pygem.file_handler as fh
import numpy as np


class TestFileHandler(TestCase):


	def test_base_class_members(self):
		file_handler = fh.FileHandler()
		assert file_handler.file == None
		assert file_handler.filename == None
		assert file_handler.extension == None

	
	def test_base_class_parse(self):
		works = True
		try:
			file_handler = fh.FileHandler()
			file_handler.parse()
			works = False
		except:
			pass
		assert works


	def test_base_class_write(self):
		works = True
		try:
			file_handler = fh.FileHandler()
			file_handler.write()
			works = False
		except:
			pass
		assert works
	

	def test_unv_file_member(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		assert unv_handler.file == None


	def test_unv_filename_member(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		assert unv_handler.filename == 'tests/test_datasets/test_square.unv'
	

	def test_unv_extension_member(self):
		unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
		assert unv_handler.extension == 'unv'


	def test_unv_passing_filename_type(self):
		works = False
		try:
			unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
			works = True
		except:
			pass
		assert works


	def test_unv_failing_filename_type(self):
		works = True
		try:
			unv_handler = fh.UnvHandler(test_square.iges)
			works = False
		except:
			pass
		assert works


	def test_unv_passing_check_extension(self):
		works = False
		try:
			unv_handler = fh.UnvHandler('tests/test_datasets/test_square.unv')
			works = True
		except:
			pass
		assert works


	def test_unv_failing_check_extension(self):
		works = True
		try:
			unv_handler = fh.UnvHandler('tests/test_datasets/test_square.iges')
			works = False
		except:
			pass
		assert works

