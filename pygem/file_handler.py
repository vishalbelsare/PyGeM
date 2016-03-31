"""
Utilities for reading and writing different CAD files.
"""
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from stl import mesh
import vtk


class FileHandler(object):
	"""
	A base class for file handling.

	:cvar string infile: name of the input file to be processed
	:cvar string outfile: name of the output file where to write in
	:cvar string extension: extension of the input/output files. It is specific for each
		subclass
	"""
	def __init__(self):
		self.infile = None
		self.outfile = None
		self.extension = None


	def parse(self, infile):
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
	def __init__(self):
		super(UnvHandler, self).__init__()
		self.extension = 'unv'


	def parse(self, filename):
		"""
		Method to parse the file `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: numpy.ndarray

		.. todo::

			- specify when it works
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)

		self.infile = filename

		input_file = open(self.infile, 'r')
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
		with open(self.infile, 'r') as input_file:
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


	def write(self, mesh_points, filename):
		"""
		Writes a unv file, called filename, copying all the lines from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the unv file.

		:param numpy.ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh
		:param string filename: name of the output file.

		.. todo:: DOCS
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)
		_check_infile_instantiation(self.infile)

		self.outfile = filename

		n_points = mesh_points.shape[0]
		nrow = 0
		i = 0
		with open(self.infile, 'r') as input_file, open(self.outfile, 'w') as output_file:
			for line in input_file:
				nrow += 1
				if nrow % 2 == 1 and nrow > 20 and nrow <= (20 + n_points * 2):
					for j in range(0, 3):
						output_file.write('   ' + str(mesh_points[i][j]))
					output_file.write('\n')
					i += 1
				elif nrow > 17:
					output_file.write(line)



class VtkHandler(FileHandler):
	"""
	Vtk file handler class

	.. todo:: DOCS
	"""
	def __init__(self):
		super(VtkHandler, self).__init__()
		self.extension = 'vtk'


	def parse(self, filename):
		"""
		Method to parse the file `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: numpy.ndarray

		.. todo::

			- specify when it works
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)

		self.infile = filename

		reader = vtk.vtkDataSetReader()
		reader.SetFileName(self.infile)
		reader.ReadAllVectorsOn()
		reader.ReadAllScalarsOn()
		reader.Update()
		data = reader.GetOutput()

		n_points = data.GetNumberOfPoints()
		mesh_points = np.zeros([n_points, 3])

		for i in range(n_points):
			mesh_points[i, 0], mesh_points[i, 1], mesh_points[i, 2] = data.GetPoint(i)

		return mesh_points


	def write(self, mesh_points, filename):
		"""
		Writes a vtk file, called filename, copying all the structures from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the vtk file.

		:param numpy.ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh
		:param string filename: name of the output file.

		.. todo:: DOCS
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)
		_check_infile_instantiation(self.infile)

		self.outfile = filename

		reader = vtk.vtkDataSetReader()
		reader.SetFileName(self.infile)
		reader.ReadAllVectorsOn()
		reader.ReadAllScalarsOn()
		reader.Update()
		data = reader.GetOutput()
	
		points = vtk.vtkPoints()
	
		for i in range(data.GetNumberOfPoints()):
			points.InsertNextPoint(mesh_points[i,:])
		
		data.SetPoints(points)
		
		writer = vtk.vtkDataSetWriter()
		writer.SetFileName(self.outfile)
		
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(data)
		else:
			writer.SetInputData(data)
		
		writer.Write()


class StlHandler(FileHandler):
	"""
	STereoLithography file handler class

	.. todo:: DOCS
	"""
	def __init__(self):
		super(StlHandler, self).__init__()
		self.extension = 'stl'


	def parse(self, filename):
		"""
		Method to parse the `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: numpy.ndarray

		.. todo::

			- specify when it works
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)

		self.infile = filename

		stl_mesh = mesh.Mesh.from_file(self.infile)
		mesh_points = np.array([stl_mesh.x.ravel(), stl_mesh.y.ravel(), stl_mesh.z.ravel()])
		mesh_points = mesh_points.T

		return mesh_points


	def write(self, mesh_points, filename):
		"""
		Writes a unv file, called filename, copying all the lines from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the unv file.

		:param numpy.ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh.
		:param string filename: name of the output file.

		.. todo:: DOCS
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)
		_check_infile_instantiation(self.infile)

		self.outfile = filename

		n_vertices = mesh_points.shape[0]
		# number of triplets of vertices
		n_triplets = n_vertices/3
		data = np.zeros(n_triplets, dtype=mesh.Mesh.dtype)
		stl_mesh = mesh.Mesh(data, remove_empty_areas=False)

		for i in range(0, n_triplets):
			for j in range(0, 3):
				data['vectors'][i][j] = mesh_points[3*i + j]

		stl_mesh.save(self.outfile, mode=1, update_normals=True)


	def plot(self, plot_file=None, save_fig=False):
		"""
		Method to plot an stl file. If `plot_file` is not given it plots `self.infile`.

		:param string plot_file: the stl filename you want to plot.
		:param bool save_fig: a flag to save the figure in png or not. If True the
			plot is not shown.
		"""
		if plot_file is None:
			plot_file = self.infile
		else:
			_check_filename_type(plot_file)

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
		if not save_fig:
			pyplot.show()
		else:
			figure.savefig(plot_file.split('.')[0] + '.png')



class OpenFoamHandler(FileHandler):
	"""
	OpenFOAM mesh file handler class.

	.. todo:: DOCS
	"""
	def __init__(self):
		super(OpenFoamHandler, self).__init__()
		self.extension = ''


	def parse(self, filename):
		"""
		Method to parse the `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: numpy.ndarray

		.. todo::

			- specify when it works
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)

		self.infile = filename

		nrow = 0
		i = 0
		with open(self.infile, 'r') as input_file:
			for line in input_file:
				nrow += 1
				if nrow == 19:
					n_points = int(line)
					mesh_points = np.zeros(shape=(n_points, 3))
				if nrow > 20 and nrow < 21 + n_points:
					line = line[line.index("(") + 1:line.rindex(")")]
					j = 0
					for number in line.split():
						mesh_points[i][j] = float(number)
						j += 1
					i += 1

		return mesh_points



	def write(self, mesh_points, filename):
		"""
		Writes a openFOAM file, called filename, copying all the lines from self.filename but
		the coordinates. mesh_points is a matrix that contains the new coordinates to
		write in the openFOAM file.

		:param numpy.ndarray mesh_points: it is a `n_points`-by-3 matrix containing
			the coordinates of the points of the mesh.
		:param string filename: name of the output file.

		.. todo:: DOCS
		"""
		_check_filename_type(filename)
		_check_extension(filename, self.extension)
		_check_infile_instantiation(self.infile)

		self.outfile = filename

		n_points = mesh_points.shape[0]
		nrow = 0
		i = 0
		with open(self.infile, 'r') as input_file, open(self.outfile, 'w') as output_file:
			for line in input_file:
				nrow += 1
				if nrow > 20 and nrow < 21 + n_points:
					output_file.write('(' + str(mesh_points[i][0]) + ' ' + str(mesh_points[i][1]) + \
									  ' ' + str(mesh_points[i][2]) +')')
					output_file.write('\n')
					i += 1
				else:
					output_file.write(line)



def _check_extension(filename, extension):
	"""
	This private method checks if the given `filename` has the proper `extension`.
	If not it raises a ValueError.

	:param string filename: file to check.
	:param string extension: file extension to check.
	"""
	file_ext = filename.split('.')[-1].lower()
	# to manage the case of open foam (no extension) the following check is needed
	if file_ext == filename.lower():
		pass
	elif not file_ext == extension:
		raise ValueError('The input file does not have the proper extension. \
			It should be %s.' % extension)


def _check_filename_type(filename):
	"""
	This private method checks if `filename` is a string. If not it raises a TypeError.

	:param string filename: file to check.
	"""
	if not isinstance(filename, basestring):
		raise TypeError('The given filename (%s) must be a string' % filename)


def _check_infile_instantiation(infile):
	"""
	This private method checks if the input file `infile` is instantiated. If not it means
	that nobody called the parse method, i.e. `self.infile` is None. If the check fails
	it raises a RuntimeError.

	:param string infile: file to check.
	"""
	if not infile:
		raise RuntimeError("You can not write a file without having parsed one.")
