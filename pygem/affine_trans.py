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
	c = np.zeros((dim+1, dim))

	for j in range(dim):
		for k in range(dim+1):
			for i in range(len(points_start)):
				qt = list(points_start[i]) + [1]
				c[k][j] += qt[k] * points_end[i][j]

	# Fill an an empty (dim+1) x (dim+1) matrix
	Q = np.zeros((dim+1, dim+1))
	for qi in points_start:
		qt = list(qi) + [1]
		for i in range(dim+1):
			for j in range(dim+1):
				Q[i][j] += qt[i] * qt[j]

	# Augement Q with c and solve Q * a' = c and put resulting matrix into the Reduced Row Echelon Form
	M = np.append(Q, c, axis=1)

	if np.linalg.cond(M) < 1/sys.float_info.epsilon:
		M, __ = syp.Matrix(M).rref()
		M = np.array(M)
	else:
		raise RuntimeError("Error: singular matrix. Points are probably coplanar.")


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
				destination[j] += source[i] * M[i][j+dim+1]
			# Add the last line of the rref
			destination[i] += M[dim][i+dim+1]
		return destination


	return transform_vector


