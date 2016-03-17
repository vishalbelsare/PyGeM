"""
Utilities for reading and writing parameters.
"""
import numpy as np


class FFDParameters():
	"""
	Class that handle the Free Form Deformation parameters in terms of FFD bounding box and weight of the FFD control points.

	:param list n_control_points: number of control points in direction (x, y, z). Default value (1, 1, 1).
	
	:cvar float length_box_x: length of the FFD bounding box in the x direction (local coordinate system).
	:cvar float length_box_y: length of the FFD bounding box in the y direction (local coordinate system).
	:cvar float length_box_z: length of the FFD bounding box in the z direction (local coordinate system).
	
	:cvar float origin_box_x: x coordinate of the origin of the FFD bounding box.
	:cvar float origin_box_y: y coordinate of the origin of the FFD bounding box.
	:cvar float origin_box_z: z coordinate of the origin of the FFD bounding box.
	
	:cvar float rot_angle_x: rotation angle around x axis of the FFD bounding box.
	:cvar float rot_angle_y: rotation angle around y axis of the FFD bounding box.
	:cvar float rot_angle_z: rotation angle around z axis of the FFD bounding box.
	
	:cvar int n_control_points_x: number of control points in direction x.
	:cvar int n_control_points_y: number of control points in direction y.
	:cvar int n_control_points_z: number of control points in direction z.
	
	:cvar numpy.ndarray array_mu_x: weights of control points in direction x.
	:cvar numpy.ndarray array_mu_y: weights of control points in direction y.
	:cvar numpy.ndarray array_mu_z: weights of control points in direction z.
	
	:cvar numpy.ndarray psi_mapping: map from the pysical domain to the reference domain.
	:cvar numpy.ndarray inv_psi_mapping: map from the reference domain to the physical domain.
	
	:cvar numpy.ndarray rot_mat: rotation matrix (according to rot_angle_x, rot_angle_y, rot_angle_z).
	
	:cvar numpy.ndarray position_vertex_0: position of the first vertex of the FFD bounding box.
	:cvar numpy.ndarray position_vertex_1: position of the second vertex of the FFD bounding box.
	:cvar numpy.ndarray position_vertex_2: position of the third vertex of the FFD bounding box.
	:cvar numpy.ndarray position_vertex_3: position of the fourth vertex of the FFD bounding box.
	
	.. note::
	
		Four vertex (non coplanar) are sufficient to uniquely identify a parallelepiped. If the four vertex are coplanar, an assert is thrown when affine_points_fit is used.
	
    """
    
	def __init__(self, n_control_points=[1, 1, 1]):
		self.conversion_unit = 1.
		
		self.lenght_box_x = 1.
		self.lenght_box_y = 1.
		self.lenght_box_z = 1.
		
		self.origin_box_x = 0.
		self.origin_box_y = 0.
		self.origin_box_z = 0.

		self.rot_angle_x = 0.
		self.rot_angle_y = 0.
		self.rot_angle_z = 0.

		self.n_control_points_x = n_control_points[0]
		self.n_control_points_y = n_control_points[1]
		self.n_control_points_z = n_control_points[2]
		
		self.array_mu_x = np.zeros((self.n_control_points_x+1, self.n_control_points_y+1, self.n_control_points_z+1))
		self.array_mu_y = np.zeros((self.n_control_points_x+1, self.n_control_points_y+1, self.n_control_points_z+1))
		self.array_mu_z = np.zeros((self.n_control_points_x+1, self.n_control_points_y+1, self.n_control_points_z+1))

		self.psi_mapping = np.diag([1./self.lenght_box_x, 1./self.lenght_box_y, 1./self.lenght_box_z])
		self.inv_psi_mapping = np.diag([self.lenght_box_x, self.lenght_box_y, self.lenght_box_z])

		self.rot_mat = np.eye(3)
		self.position_vertex_0 = np.zeros(3)
		self.position_vertex_1 = np.zeros(3)
		self.position_vertex_2 = np.zeros(3)
		self.position_vertex_3 = np.zeros(3)


	def print_info(self):
		"""
		This method prints all the FFD parameters on the screen.
		"""
		print 'conversion_unit = ' + str(self.conversion_unit) + '\n'
		print '(lenght_box_x, lenght_box_y, lenght_box_z) = (' + str(self.lenght_box_x) + ', ' + str(self.lenght_box_y) + ', ' + \
		str(self.lenght_box_z) + ')'
		print '(origin_box_x, origin_box_y, origin_box_z) = (' + str(self.origin_box_x) + ', ' + str(self.origin_box_y) + ', ' + \
		str(self.origin_box_z) + ')'
		print '(n_control_points_x, n_control_points_y, n_control_points_z) = (' + str(self.n_control_points_x) + ', ' + str(self.n_control_points_y) + ', ' + \
		str(self.n_control_points_z) + ')'
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
		print '\nrot_mat ='
		print self.rot_mat
		print '\nposition_vertex_0 ='
		print self.position_vertex_0
		print '\nposition_vertex_1 ='
		print self.position_vertex_1
		print '\nposition_vertex_2 ='
		print self.position_vertex_2
		print '\nposition_vertex_3 ='
		print self.position_vertex_3


