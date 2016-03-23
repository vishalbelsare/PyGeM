"""
Utilities for reading and writing different CAD files.
"""
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from stl import mesh
#import os


class FileHandler(object):
	"""
	A base class for file handling.

	:cvar string filename: name of the file to be processed
	:cvar string extension: extension of the file to be processed
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
		:rtype: float numpy.ndarray

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
		:param string outfile: name of the output file.

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



class StlHandler(FileHandler):
	"""
	STereoLithography file handler class

	.. todo:: DOCS
	"""
	def __init__(self, filename):
		super(StlHandler, self).__init__()

		if not isinstance(filename, basestring):
			raise TypeError("filename must be a string")

		self.filename = filename
		self.extension = self._extract_extension()
		self._check_extension('stl')


	def parse(self):
		"""
		Method to parse the `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: numpy.ndarray

		.. todo::

			- specify when it works
		"""
		stl_mesh = mesh.Mesh.from_file(self.filename)
		mesh_points = np.array([stl_mesh.x.ravel(), stl_mesh.y.ravel(), stl_mesh.z.ravel()])
		mesh_points = mesh_points.T

		return mesh_points


	def write(self, mesh_points, outfile):
		"""
		Writes a unv file, called outfile, copying all the lines from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the unv file.

		:param numpy.ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh.
		:param string outfile: name of the output file.

		.. todo:: DOCS
		"""
		if not isinstance(outfile, basestring):
			raise TypeError("outfile must be a string")

		n_vertices = mesh_points.shape[0]
		# number of triplets of vertices
		n_triplets = n_vertices/3
		data = np.zeros(n_triplets, dtype=mesh.Mesh.dtype)
		stl_mesh = mesh.Mesh(data, remove_empty_areas=False)

		for i in range(0, n_triplets):
			for j in range(0, 3):
				data['vectors'][i][j] = mesh_points[3*i + j]

		stl_mesh.save(outfile, mode=1, update_normals=True)


	def plot(self, plot_file=None):
		"""
		Method to plot an stl file. If `plot_file` is not given it plots `self.filename`.

		:param string plot_file: the stl filename you want to plot.
		"""
		if plot_file is None:
			plot_file = self.filename
		else:
			if not isinstance(plot_file, basestring):
				raise TypeError("plot_file must be a string")

		# Create a new plot
		figure = pyplot.figure()
		axes = mplot3d.Axes3D(figure)

		# Load the STL files and add the vectors to the plot
		stl_mesh = mesh.Mesh.from_file(plot_file)
		axes.add_collection3d(mplot3d.art3d.Poly3DCollection(stl_mesh.vectors))

		# Auto scale to the mesh size
		scale = stl_mesh.points.flatten(-1)
		axes.auto_scale_xyz(scale, scale, scale)

		# Show the plot to the screen
		pyplot.show()
		
		
class OpenFoamHandler(FileHandler):
	"""
	OpenFOAM mesh file handler class.

	.. todo:: DOCS
	"""
	def __init__(self, filename):
		super(OpenFoamHandler, self).__init__()

		if not isinstance(filename, basestring):
			raise TypeError("filename must be a string")

		self.filename = filename


	def parse(self):
		"""
		Method to parse the `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: numpy.ndarray

		.. todo::

			- specify when it works
		"""
		
		nrow = 0
		i = 0
		with open(self.filename, 'r') as input_file:
			for line in input_file:
				nrow += 1
				if nrow == 19:
					n_points = int(line)
					mesh_points = np.zeros(shape=(n_points,3))
				if nrow > 20 and nrow < 21 + n_points:
					line = line[line.index("(") + 1:line.rindex(")")]
					j = 0
					for number in line.split():
						mesh_points[i][j] = float(number)
						j += 1
					i += 1

		return mesh_points



	def write(self, mesh_points, outfile):
		"""
		Writes a openFOAM file, called outfile, copying all the lines from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the openFOAM file.

		:param numpy.ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh.
		:param string outfile: name of the output file.
		
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
				if nrow > 20 and nrow < 21 + n_points:
					output_file.write('(' + str(mesh_points[i][0]) + ' ' + str(mesh_points[i][1]) + ' ' + str(mesh_points[i][2]) +')')
					output_file.write('\n')
					i += 1
				else:
					output_file.write(line)
					
					
