"""
Utilities for performing Free Form Deformation (FFD)
"""
import numpy as np
from scipy import special
import pygem.ffd_parameters as ffdp
import pygem.affine_trans as at


class FFD(object):
	"""
	DOCS
    """

	def __init__(self, ffd_parameters, original_mesh_points):
		self.parameters = ffd_parameters
		self.original_mesh_points = original_mesh_points


	def perform(self):
		"""
		DOCS
		"""
		(cm, cn) = self.original_mesh_points.shape

		## translation and then affine transformation
		translation = np.array([self.parameters.origin_box_x, self.parameters.origin_box_y, self.parameters.origin_box_z])

		v0 = np.array([self.parameters.position_vertex_1 - translation, \
			  self.parameters.position_vertex_2 - translation, \
			  self.parameters.position_vertex_3 - translation, \
			  [0,0,0]])
		v1 = np.array([[1,0,0], [0,1,0], [0,0,1], [0,0,0]])

		transformation = at.affine_points_fit(v0, v1)
		inverse_transformation = at.affine_points_fit(v1, v0)
	
		##
		(n, m, t) = self.parameters.array_mu_x.shape
		bx = np.zeros((cm, n))
		by = np.zeros((cm, m))
		bz = np.zeros((cm, t))
		Cp_tilde = np.zeros((cm, cn))

		mesh_points_mod = self.original_mesh_points - translation	

		Cp0 = self._transform_points(mesh_points_mod, transformation)

		for i in range  (0, n):
			aux1 = np.power((1-Cp0[:,0]),n-1-i)
			aux2 = np.power(Cp0[:,0],i)
			bx[:,i] = special.binom(n-1,i) * np.multiply(aux1, aux2)

		for i in range  (0, m):
			aux1 = np.power((1-Cp0[:,1]),m-1-i)
			aux2 = np.power(Cp0[:,1],i)
			by[:,i] = special.binom(m-1,i) * np.multiply(aux1, aux2)

		for i in range  (0, t):
			aux1 = np.power((1-Cp0[:,2]),t-1-i)
			aux2 = np.power(Cp0[:,2],i)
			bz[:,i] = special.binom(t-1,i) * np.multiply(aux1, aux2)

		for i in range (0, n):
			for j in range (0, m):
				for k in range (0, t):
					aux = np.multiply(bx[:,i], np.multiply(by[:,j], bz[:,k]))
					aux_x = aux * self.parameters.array_mu_x[i, j, k]
					aux_y = aux * self.parameters.array_mu_y[i, j, k]
					aux_z = aux * self.parameters.array_mu_z[i, j, k]
					Cp_tilde[:,0] += aux_x
					Cp_tilde[:,1] += aux_y
					Cp_tilde[:,2] += aux_z

		# Splitting points inside and outside the lattice: TODO not very efficient
		for i in range (0, cm):
			if (Cp0[i,0] < 0) or (Cp0[i,1] < 0) or (Cp0[i,2] < 0) or \
				(Cp0[i,0] > 1) or (Cp0[i,1] > 1) or (Cp0[i,2] > 1):
				Cp_tilde[i,:] = 0

		return self._transform_points(Cp_tilde + Cp0, inverse_transformation) + translation


	def _transform_points(self, mesh_points, transformation):
		"""
		DOCS
		"""
		n, m = mesh_points.shape
		Cp = np.zeros((n, m))

		for i in range(0, n):
			fp = mesh_points[i]
			Cp[i,:] = transformation(fp)

		return Cp
