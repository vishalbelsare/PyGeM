"""
Utilities for reading and writing different CAD files.
"""
import numpy as np
import vtk
import pygem.filehandler as fh


class VtkHandler(fh.FileHandler):
	"""
	Vtk file handler class

	:cvar string infile: name of the input file to be processed.
	:cvar string outfile: name of the output file where to write in.
	:cvar string extension: extension of the input/output files. It is equal to '.vtk'.
	"""
	def __init__(self):
		super(VtkHandler, self).__init__()
		self.extension = '.vtk'


	def parse(self, filename):
		"""
		Method to parse the file `filename`. It returns a matrix with all the coordinates.

		:return: mesh_points: it is a `n_points`-by-3 matrix containing the coordinates of
			the points of the mesh
		:rtype: numpy.ndarray

		.. todo::

			- specify when it works
		"""
		self._check_filename_type(filename)
		self._check_extension(filename)

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
		self._check_filename_type(filename)
		self._check_extension(filename)
		self._check_infile_instantiation(self.infile)

		self.outfile = filename

		reader = vtk.vtkDataSetReader()
		reader.SetFileName(self.infile)
		reader.ReadAllVectorsOn()
		reader.ReadAllScalarsOn()
		reader.Update()
		data = reader.GetOutput()

		points = vtk.vtkPoints()

		for i in range(data.GetNumberOfPoints()):
			points.InsertNextPoint(mesh_points[i, :])

		data.SetPoints(points)

		writer = vtk.vtkDataSetWriter()
		writer.SetFileName(self.outfile)

		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(data)
		else:
			writer.SetInputData(data)

		writer.Write()


	def plot(self, plot_file=None):
		"""
		Method to plot an stl file. If `plot_file` is not given it plots `self.infile`.

		:param string plot_file: the stl filename you want to plot.
		"""
		if plot_file is None:
			plot_file = self.infile
		else:
			self._check_filename_type(plot_file)

		# Read the source file.
		reader = vtk.vtkUnstructuredGridReader()
		reader.SetFileName(plot_file)
		reader.Update() # Needed because of GetScalarRange
		output = reader.GetOutput()
		scalar_range = output.GetScalarRange()
		 
		# Create the mapper that corresponds the objects of the vtk file
		# into graphics elements
		mapper = vtk.vtkDataSetMapper()
		if vtk.VTK_MAJOR_VERSION <= 5:
			mapper.SetInput(output)
		else:
			mapper.SetInputData(output)
		mapper.SetScalarRange(scalar_range)
		 
		# Create the Actor
		actor = vtk.vtkActor()
		actor.SetMapper(mapper)
		 
		# Create the Renderer
		renderer = vtk.vtkRenderer()
		renderer.AddActor(actor)
		renderer.SetBackground(20, 20, 20) # Set background color (white is 1, 1, 1)
		 
		# Create the RendererWindow
		renderer_window = vtk.vtkRenderWindow()
		renderer_window.AddRenderer(renderer)
		 
		# Create the RendererWindowInteractor and display the vtk_file
		interactor = vtk.vtkRenderWindowInteractor()
		interactor.SetRenderWindow(renderer_window)
		interactor.Initialize()
		interactor.Start()

		