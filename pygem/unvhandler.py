"""
Derived module from filehandler.py to handle Universal (unv) files.
"""
import numpy as np
import pygem.filehandler as fh


class UnvHandler(fh.FileHandler):
	"""
	Universal file handler class

	:cvar string infile: name of the input file to be processed.
	:cvar string outfile: name of the output file where to write in.
	:cvar string extension: extension of the input/output files. It is equal to '.unv'.
	"""
	def __init__(self):
		super(UnvHandler, self).__init__()
		self.extension = '.unv'


	def parse(self, filename):
		"""
		Method to parse the file `filename`. It returns a matrix with all the coordinates.
		It reads only the section 2411 of the unv files and it assumes there are only triangles.

		:param string filename: name of the input file.
		
		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh.
		:rtype: numpy.ndarray
		"""
		self._check_filename_type(filename)
		self._check_extension(filename)

		self.infile = filename

		with open(self.infile, 'r') as input_file:
			nline = 0
			while True:
				line = input_file.readline()
				nline += 1
				if not line:
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

		n_points = count/2
		mesh_points = np.zeros(shape=(n_points, 3))

		nline = 0
		i = 0
		with open(self.infile, 'r') as input_file:
			for line in input_file:
				nline += 1
				if nline % 2 == 1 and start_line < nline < last_line:
					line = line.strip()
					j = 0
					for number in line.split():
						mesh_points[i][j] = float(number)
						j += 1
					i += 1

		return mesh_points


	def write(self, mesh_points, filename):
		"""
		Writes a unv file, called filename, copying all the lines from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the unv file.

		:param numpy.ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh
		:param string filename: name of the output file.
		"""
		self._check_filename_type(filename)
		self._check_extension(filename)
		self._check_infile_instantiation(self.infile)

		self.outfile = filename

		n_points = mesh_points.shape[0]
		nrow = 0
		i = 0
		with open(self.infile, 'r') as input_file, open(self.outfile, 'w') as output_file:
			for line in input_file:
				nrow += 1
				if nrow % 2 == 1 and 20 < nrow <= (20 + n_points * 2):
					for j in range(0, 3):
						output_file.write('   ' + str(mesh_points[i][j]))
					output_file.write('\n')
					i += 1
				elif nrow > 17:
					output_file.write(line)
