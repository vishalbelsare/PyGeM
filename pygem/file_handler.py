"""
Utilities for reading and writing different CAD files.
"""
import numpy as np


class FileHandler(object):
	"""
	A base class for file handling.

	:cvar string filename: name of the file to be processed
	:cvar string extension: extension of the file to be processed

	.. todo:: Add virtual methods like parse and write
	"""
	def __init__(self):
		self.filename = None
		self.extension = None


	def _extract_extension(self):
		"""
		It is used in child classes.
		
		.. todo:: DOCS
		"""
		ext = self.filename.split('.')[-1].lower()
		return ext


	def _check_extension(self, extension):
		"""
		Internal function that checks that the member extension is equal to
		a given extension.
		It is used in child classes.

		:param string extension: file extension to check.
		"""
		if not self.extension == extension:
			raise ValueError('The input file does not have the proper extension. \
				It should be %s.' % extension)


	def parse(self):
		"""
		Abstract method to parse a specific file.

		Not implemented, it has to be implemented in subclasses.
		"""
		raise NotImplementedError("Subclass must implement abstract method " \
			+ self.__class__.__name__ + ".parse")


	def write(self, mesh_points, outfile):
		"""
		Abstract method to write a specific file.

		Not implemented, it has to be implemented in subclasses.
		"""
		raise NotImplementedError("Subclass must implement abstract method " \
			+ self.__class__.__name__ + ".write")



class UnvHandler(FileHandler):
	"""
	Universal file handler class

	.. todo:: DOCS
	"""
	def __init__(self, filename):
		super(UnvHandler, self).__init__()

		if not isinstance(filename, basestring):
			raise TypeError("filename must be a string")

		self.filename = filename
		self.extension = self._extract_extension()
		self._check_extension('unv')


	def parse(self):
		"""
		Method to parse the `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: float ndarray

		.. todo::

			- specify when it works
		"""

		input_file = open(self.filename, 'r')
		nline = 0
		while True:
			line = input_file.readline()
			nline += 1
			if len(line) == 0:
				break
			if line.startswith('    -1'):
				section_id = input_file.readline().strip()
				nline += 1
				if section_id == '2411':
					count = 0
					while not input_file.readline().startswith('    -1'):
						count += 1
					start_line = nline + 2
					last_line = start_line + count
				else:
					while not input_file.readline().startswith('    -1'):
						nline += 1

		input_file.close()

		n_points = count/2
		mesh_points = np.zeros(shape=(n_points, 3))

		nline = 0
		i = 0
		with open(self.filename, 'r') as input_file:
			for line in input_file:
				nline += 1
				if nline % 2 == 1 and nline > start_line and nline < last_line:
					line = line.strip()
					j = 0
					for number in line.split():
						mesh_points[i][j] = float(number)
						j += 1
					i += 1

		return mesh_points


	def write(self, mesh_points, outfile):
		"""
		Writes a unv file, called outfile, copying all the lines from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the unv file.

		:param ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh

		.. todo:: DOCS
		"""
		if not isinstance(outfile, basestring):
			raise TypeError("outfile must be a string")

		n_points = mesh_points.shape[0]
		nrow = 0
		i = 0
		with open(self.filename, 'r') as input_file, open(outfile, 'w') as output_file:
			for line in input_file:
				nrow += 1
				if nrow % 2 == 1 and nrow > 20 and nrow <= (20 + n_points * 2):
					for j in range(0, 3):
						output_file.write('   ' + str(mesh_points[i][j]))
					output_file.write('\n')
					i += 1
				elif nrow > 17:
					output_file.write(line)


