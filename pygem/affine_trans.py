"""
Utilities for the affine transformations of the bounding box of the Free Form Deformation.
"""
import math
import sys
import numpy as np
import sympy as syp


def angles2matrix(rot_z=0, rot_y=0, rot_x=0):
	"""
	Returns matrix for given rotations around z, y and x axes. The output rotation matrix is
	equal to the composition of the	individual rotations. Rotations are counter-clockwise.
	The default value of the three rotations is zero.

	:param float rot_z: rotation angle (in radians) around z-axis.
	:param float rot_y: rotation angle (in radians) around y-axis.
	:param float rot_x: rotation angle (in radians) around x-axis.

	:return: rot_matrix: rotation matrix for the given angles. The matrix shape is always (3, 3).
	:rtype: float

	.. note::

		- The direction of rotation is given by the right-hand rule.
		- When applying the rotation to a vector, the vector should be column vector
		  to the right of the rotation matrix.

	"""
	rot_matrix = []
	if rot_z:
		cos = math.cos(rot_z)
		sin = math.sin(rot_z)
		rot_matrix.append(np.array([cos, -sin, 0, sin, cos, 0, 0, 0, 1]).reshape((3, 3)))
	if rot_y:
		cos = math.cos(rot_y)
		sin = math.sin(rot_y)
		rot_matrix.append(np.array([cos, 0, sin, 0, 1, 0, -sin, 0, cos]).reshape((3, 3)))
	if rot_x:
		cos = math.cos(rot_x)
		sin = math.sin(rot_x)
		rot_matrix.append(np.array([1, 0, 0, 0, cos, -sin, 0, sin, cos]).reshape((3, 3)))
	if rot_matrix:
		return reduce(np.dot, rot_matrix[::-1])
	return np.eye(3)


def affine_points_fit(points_start, points_end):
	"""
	Fit an affine transformation from starting points to ending points through a
	least square procedure.

	:param numpy.ndarray points_start: set of starting points.
	:param numpy.ndarray points_end: set of ending points.

	:return: transform_vector: function that transforms a vector according to the
			 affine map. It takes a source vector and return a vector transformed
			 by the reduced row echelon form of the map.
	:rtype: function

	:Example:

	>>> import pygem.affine_trans as at
	>>> # Example of a rotation
	>>> p_start = np.array([[1,0,0], [0,1,0], [0,0,1], [0,0,0]])
	>>> p_end = np.array([[0,1,0], [-1,0,0], [0,0,1], [0,0,0]])
	>>> v_test = np.array([1., 2., 3.])
	>>> transformation = at.affine_points_fit(p_start, p_end)
	>>> v_trans = transformation(v_test)
	"""
	if len(points_start) != len(points_end):
		raise RuntimeError("points_start and points_end must be of same size.")

	dim = len(points_start[0])
	if len(points_start) < dim:
		raise RuntimeError("Too few starting points => under-determined system.")

    # Fill an an empty (dim+1) x (dim) matrix
	c = [[0.0 for a in range(dim)] for i in range(dim+1)]
	for j in range(dim):
		for k in range(dim+1):
			for i in range(len(points_start)):
				qt = list(points_start[i]) + [1]
				c[k][j] += qt[k] * points_end[i][j]

	# Fill an an empty (dim+1) x (dim+1) matrix
	Q = [[0.0 for a in range(dim)] + [0] for i in range(dim+1)]
	for qi in points_start:
		qt = list(qi) + [1]
		for i in range(dim+1):
			for j in range(dim+1):
				Q[i][j] += qt[i] * qt[j]

	
	# Ultra simple linear system solver.
	def gauss_jordan(m, eps = 1.0/(10**10)):
		"""
		TO DOC

		Puts given matrix (2D array) into the Reduced Row Echelon Form.
		Returns True if successful, False if 'm' is singular.
		NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
		Written by Jarno Elonen in April 2005, released into Public Domain"""
		(h, w) = (len(m), len(m[0]))
		
		for y in range(0,h):
			maxrow = y
			for y2 in range(y+1, h):    # Find max pivot
				if abs(m[y2][y]) > abs(m[maxrow][y]):
					maxrow = y2
			(m[y], m[maxrow]) = (m[maxrow], m[y])
			if abs(m[y][y]) <= eps:     # Singular?
				raise RuntimeError("Error: singular matrix. Points are probably coplanar.")
			
			for y2 in range(y+1, h):    # Eliminate column y
				c = m[y2][y] / m[y][y]
				for x in range(y, w):
					m[y2][x] -= m[y][x] * c
		
		for y in range(h-1, 0-1, -1): # Backsubstitute
			c  = m[y][y]
			for y2 in range(0,y):
				for x in range(w-1, y-1, -1):
					m[y2][x] -=  m[y][x] * m[y2][y] / c
			m[y][y] /= c
			for x in range(h, w):       # Normalize row y
				m[y][x] /= c
		return m

	
	# Augement Q with c and solve Q * a' = c by Gauss-Jordan
	M = [ Q[i] + c[i] for i in range(dim+1)]
	rref_M = gauss_jordan(M)

	def transform_vector(source):
		"""
		Transform a vector according to the affine map.

		:param numpy.ndarray source: vector to be transformed.

		:return destination: numpy.ndarray representing the transformed vector.
		:rtype: float
		"""
		destination = np.zeros(dim)
		for i in range(dim):
			for j in range(dim):
				destination[j] += source[i] * rref_M[i][j+dim+1]
			# Add the last line of the rref
			destination[i] += rref_M[dim][i+dim+1]
		return destination


	return transform_vector
	

