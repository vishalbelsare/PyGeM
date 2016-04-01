
from unittest import TestCase
import unittest
import pkgutil
from os import walk


class TestPackage(TestCase):


	def test_modules_name(self):
		import pygem
		package = pygem
		mod = ['__init__.py']
		for __, modname, __ in pkgutil.iter_modules(package.__path__):
			mod.append(modname + '.py')
		f = []
		for (__, __, filenames) in walk('pygem'):
			f.extend(filenames)
		self.assertItemsEqual(mod, f)


	def test_import_pg_1(self):
		import pygem as pg
		params = pg.params.FFDParameters()


	def test_import_pg_2(self):
		import pygem as pg
		mat = pg.affine.angles2matrix(2)