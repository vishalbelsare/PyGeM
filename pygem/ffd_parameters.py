"""
Utilities for reading and writing parameters files to perform the Free Form Deformation (FFD)
"""
import os
import numpy as np
import pygem.affine_trans as at


class FFDParameters(object):
	"""
	Class that handles the Free Form Deformation parameters in terms of FFD bounding box and
	weight of the FFD control points.

	:param list n_control_points: number of control points in the x, y, and z direction.
		If not provided it is set to [1, 1, 1].

	:cvar float length_box_x: length of the FFD bounding box in the x direction
		(local coordinate system).
	:cvar float length_box_y: length of the FFD bounding box in the y direction
		(local coordinate system).
	:cvar float length_box_z: length of the FFD bounding box in the z direction
		(local coordinate system).

	:cvar numpy.ndarray origin_box: a 3-by-1 vector of float numbers representing the x, y and z
		coordinates of the origin of the FFD bounding box.

	:cvar float rot_angle_x: rotation angle around x axis of the FFD bounding box.
	:cvar float rot_angle_y: rotation angle around y axis of the FFD bounding box.
	:cvar float rot_angle_z: rotation angle around z axis of the FFD bounding box.

	:cvar list n_control_points: list of 3 int representing the number of control points in the
		x, y, and z direction.

	:cvar numpy.ndarray array_mu_x: weights of control points in direction x.
	:cvar numpy.ndarray array_mu_y: weights of control points in direction y.
	:cvar numpy.ndarray array_mu_z: weights of control points in direction z.

	:cvar numpy.ndarray psi_mapping: map from the pysical domain to the reference domain.
	:cvar numpy.ndarray inv_psi_mapping: map from the reference domain to the physical domain.

	:cvar numpy.ndarray rotation_matrix: rotation matrix (according to rot_angle_x,
		rot_angle_y, rot_angle_z).

	:cvar numpy.ndarray position_vertex_0: position of the first vertex of the FFD bounding box.
		It is always equal to the member `origin_box`.
	:cvar numpy.ndarray position_vertex_1: position of the second vertex of the FFD bounding box.
	:cvar numpy.ndarray position_vertex_2: position of the third vertex of the FFD bounding box.
	:cvar numpy.ndarray position_vertex_3: position of the fourth vertex of the FFD bounding box.

	.. note::
		Four vertex (non coplanar) are sufficient to uniquely identify a parallelepiped.
		If the four vertex are coplanar, an assert is thrown when affine_points_fit is used.


	:Example:

	>>> import pygem.ffd_parameters as ffdp

	>>> # Reading an existing file
	>>> params1 = ffdp.FFDParameters()
	>>> params1.read_parameters_file(filename='tests/test_datasets/parameters_test_ffd_identity.prm')

	>>> # Creating a defaul paramters file with the right dimensions (if the file does not exists
	>>> # it is created with that name). So it is possible to manually edit it and read it again.
	>>> params2 = ffdp.FFDParameters(n_control_points=[2, 3, 2])
	>>> params2.read_parameters_file(filename='parameters_test.prm')
    """

	def __init__(self, n_control_points=None):
		self.conversion_unit = 1.

		self.lenght_box_x = 1.
		self.lenght_box_y = 1.
		self.lenght_box_z = 1.

		self.origin_box = np.array([0., 0., 0.])

		self.rot_angle_x = 0.
		self.rot_angle_y = 0.
		self.rot_angle_z = 0.

		if n_control_points is None:
			n_control_points = [1, 1, 1]
		self.n_control_points = n_control_points

		self.array_mu_x = np.zeros((self.n_control_points[0]+1, self.n_control_points[1]+1, \
			self.n_control_points[2]+1))
		self.array_mu_y = np.zeros((self.n_control_points[0]+1, self.n_control_points[1]+1, \
			self.n_control_points[2]+1))
		self.array_mu_z = np.zeros((self.n_control_points[0]+1, self.n_control_points[1]+1, \
			self.n_control_points[2]+1))

		self.psi_mapping = np.diag([1./self.lenght_box_x, 1./self.lenght_box_y, 1./self.lenght_box_z])
		self.inv_psi_mapping = np.diag([self.lenght_box_x, self.lenght_box_y, self.lenght_box_z])

		self.rotation_matrix = np.eye(3)
		self.position_vertex_0 = self.origin_box
		self.position_vertex_1 = np.array([1., 0., 0.])
		self.position_vertex_2 = np.array([0., 1., 0.])
		self.position_vertex_3 = np.array([0., 0., 1.])


	def read_parameters_file(self, filename=None):
		"""
		Reads in the parameters file and fill the self structure.

		:param string filename: parameters file to be read in.

		"""
		if not isinstance(filename, basestring):
			raise TypeError("filename must be a string")

		if filename is None:
			filename = 'parameters.prm'

		# Checks if the parameters file exists. If not it writes the default class into filename.
		if not os.path.isfile(filename):
			self.write_parameters_file(filename)
			return

		found_x = found_y = found_z = False
		with open(filename, 'r') as input_file:
			for line in input_file:
				if found_x and line.strip() != '' and line.strip() != 'mu_y:':
					values = line.split()
					self.array_mu_x[int(values[0])][int(values[1])][int(values[2])] = float(values[3])
				if found_y and line.strip() != '' and line.strip() != 'mu_z:':
					values = line.split()
					self.array_mu_y[int(values[0])][int(values[1])][int(values[2])] = float(values[3])
				if found_z:
					values = line.split()
					self.array_mu_z[int(values[0])][int(values[1])][int(values[2])] = float(values[3])
				if line.strip() == 'mu_x:':
					found_x = True
					self.array_mu_x = np.zeros((self.n_control_points[0]+1, self.n_control_points[1]+1, \
						self.n_control_points[2]+1))
					self.array_mu_y = np.zeros((self.n_control_points[0]+1, self.n_control_points[1]+1, \
						self.n_control_points[2]+1))
					self.array_mu_z = np.zeros((self.n_control_points[0]+1, self.n_control_points[1]+1, \
						self.n_control_points[2]+1))
				elif line.strip() == 'mu_y:':
					found_y = True
					found_x = False
				elif line.strip() == 'mu_z:':
					found_z = True
					found_y = False
				else:
					for spl in line.split():
						if spl == 'n_control_points_x:':
							self.n_control_points[0] = int(line.strip(spl))
						if spl == 'n_control_points_y:':
							self.n_control_points[1] = int(line.strip(spl))
						if spl == 'n_control_points_z:':
							self.n_control_points[2] = int(line.strip(spl))
						if spl == 'lenght_box_x:':
							self.lenght_box_x = float(line.strip(spl))
						if spl == 'lenght_box_y:':
							self.lenght_box_y = float(line.strip(spl))
						if spl == 'lenght_box_z:':
							self.lenght_box_z = float(line.strip(spl))
						if spl == 'origin_box_x:':
							self.origin_box[0] = float(line.strip(spl))
						if spl == 'origin_box_y:':
							self.origin_box[1] = float(line.strip(spl))
						if spl == 'origin_box_z:':
							self.origin_box[2] = float(line.strip(spl))
						if spl == 'rot_angle_x:':
							self.rot_angle_x = float(line.strip(spl))
						if spl == 'rot_angle_y:':
							self.rot_angle_y = float(line.strip(spl))
						if spl == 'rot_angle_z:':
							self.rot_angle_z = float(line.strip(spl))

			self.rotation_matrix = at.angles2matrix(self.rot_angle_z*np.pi/180, \
				self.rot_angle_y*np.pi/180, self.rot_angle_x*np.pi/180)

			self.position_vertex_0 = self.origin_box
			self.position_vertex_1 = self.position_vertex_0 + \
				np.dot(self.rotation_matrix, [self.lenght_box_x, 0, 0])
			self.position_vertex_2 = self.position_vertex_0 + \
				np.dot(self.rotation_matrix, [0, self.lenght_box_y, 0])
			self.position_vertex_3 = self.position_vertex_0 + \
				np.dot(self.rotation_matrix, [0, 0, self.lenght_box_z])

			self.psi_mapping = np.diag([1./self.lenght_box_x, 1./self.lenght_box_y, 1./self.lenght_box_z])
			self.inv_psi_mapping = np.diag([self.lenght_box_x, self.lenght_box_y, self.lenght_box_z])


	def write_parameters_file(self, filename=None):
		"""
		This method writes a parameters file (.prm) called `filename` and fills it with all
		the parameters class members.

		.. todo::
			document better

		"""
		if filename is None:
			filename = 'parameters.prm'

		with open(filename, 'w') as output_file:
			output_file.write('\nn_control_points_x: ' + str(self.n_control_points[0]) + '\n')
			output_file.write('n_control_points_y: ' + str(self.n_control_points[1]) + '\n')
			output_file.write('n_control_points_z: ' + str(self.n_control_points[2]) + '\n')

			output_file.write('\nlenght_box_x: ' + str(self.lenght_box_x) + '\n')
			output_file.write('lenght_box_y: ' + str(self.lenght_box_y) + '\n')
			output_file.write('lenght_box_z: ' + str(self.lenght_box_z) + '\n')

			output_file.write('\norigin_box_x: ' + str(self.origin_box[0]) + '\n')
			output_file.write('origin_box_y: ' + str(self.origin_box[1]) + '\n')
			output_file.write('origin_box_z: ' + str(self.origin_box[2]) + '\n')

			output_file.write('\nrot_angle_x: ' + str(self.rot_angle_x) + '\n')
			output_file.write('rot_angle_y: ' + str(self.rot_angle_y) + '\n')
			output_file.write('rot_angle_z: ' + str(self.rot_angle_z) + '\n')

			output_file.write('\nmu_x:\n')
			for i in range(0, self.n_control_points[0] + 1):
				for j in range(0, self.n_control_points[1] + 1):
					for k in range(0, self.n_control_points[2] + 1):
						output_file.write(str(i) + '   ' + str(j) + '   ' + str(k) + \
						'   ' + str(self.array_mu_x[i][j][k]) + '\n')
			output_file.write('\nmu_y:\n')
			for i in range(0, self.n_control_points[0] + 1):
				for j in range(0, self.n_control_points[1] + 1):
					for k in range(0, self.n_control_points[2] + 1):
						output_file.write(str(i) + '   ' + str(j) + '   ' + str(k) + \
						'   ' + str(self.array_mu_y[i][j][k]) + '\n')
			output_file.write('\nmu_z:\n')
			for i in range(0, self.n_control_points[0] + 1):
				for j in range(0, self.n_control_points[1] + 1):
					for k in range(0, self.n_control_points[2] + 1):
						output_file.write(str(i) + '   ' + str(j) + '   ' + str(k) + \
						'   ' + str(self.array_mu_z[i][j][k]) + '\n')


	def print_info(self):
		"""
		This method prints all the FFD parameters on the screen.
		"""
		print 'conversion_unit = ' + str(self.conversion_unit) + '\n'
		print '(lenght_box_x, lenght_box_y, lenght_box_z) = (' + str(self.lenght_box_x) + \
			', ' + str(self.lenght_box_y) + ', ' + \
		str(self.lenght_box_z) + ')'
		print 'origin_box = ' + str(self.origin_box)
		print 'n_control_points = ' + str(self.n_control_points)
		print '(rot_angle_x, rot_angle_y, rot_angle_z) = (' + str(self.rot_angle_x) + \
			', ' + str(self.rot_angle_y) + ', ' + \
		str(self.rot_angle_z) + ')'
		print '\narray_mu_x ='
		print self.array_mu_x
		print '\narray_mu_y ='
		print self.array_mu_y
		print '\narray_mu_z ='
		print self.array_mu_z
		print '\npsi_mapping ='
		print self.psi_mapping
		print '\ninv_psi_mapping ='
		print self.inv_psi_mapping
		print '\nrotation_matrix ='
		print self.rotation_matrix
		print '\nposition_vertex_0 ='
		print self.position_vertex_0
		print '\nposition_vertex_1 ='
		print self.position_vertex_1
		print '\nposition_vertex_2 ='
		print self.position_vertex_2
		print '\nposition_vertex_3 ='
		print self.position_vertex_3


