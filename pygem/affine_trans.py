"""
Utilities for the affine transformations of the bounding box of the Free Form Deformation.
"""

import math
import numpy as np
import sys
import sympy as syp

def angles2matrix(z=0, y=0, x=0):
	"""
	Returns matrix for given rotations around z, y and x axes. The output rotation matrix is equal to the composition of the
	individual rotations. Rotations are counter-clockwise. The default value of the three rotations is zero.

	:param float z: rotation angle (in radians) around z-axis.
	:param float y: rotation angle (in radians) around y-axis.
	:param float x: rotation angle (in radians) around x-axis.

	:return: rot_matrix: rotation matrix for the given angles. The matrix shape is always (3,3).
	:rtype: float numpy.ndarray
	
	.. note::
	
		- The direction of rotation is given by the right-hand rule.
		- When applying the rotation to a vector, the vector should be column vector to the right of the rotation matrix.
	
    """

	rot_matrix = []
	if z:
		cos = math.cos(z)
		sin = math.sin(z)
		rot_matrix.append(np.array([cos, -sin, 0, sin, cos, 0, 0, 0, 1]).reshape((3,3)))
	if y:
		cos = math.cos(y)
		sin = math.sin(y)
		rot_matrix.append(np.array([cos, 0, sin, 0, 1, 0, -sin, 0, cos]).reshape((3,3)))
	if x:
		cos = math.cos(x)
		sin = math.sin(x)
		rot_matrix.append(np.array([1, 0, 0, 0, cos, -sin, 0, sin, cos]).reshape((3,3)))
	if rot_matrix:
		return reduce(np.dot, rot_matrix[::-1])
	return np.eye(3)
	
	
def affine_points_fit(points_start, points_end):
	"""
	Fit an affine transformation from starting points to ending points through a least square procedure.
	
	:param numpy.ndarray points_start: set of starting points.
	:param numpy.ndarray points_end: set of ending points.
	
	:return: AffineTrasformation: class containing the affine transformation object and its method to evaluate the transformation.
	:rtype: AffineTransformation
	"""

	q = points_start
	p = points_end
	if len(q) != len(p) or len(q)<1:
		raise RuntimeError("points_start and points_end must be of same size.")

	dim = len(q[0])
	if len(q) < dim:
		raise RuntimeError("Too few points => under-determined system.")

    # Fill an an empty (dim+1) x (dim) matrix
	c = np.zeros((dim+1,dim))

	for j in range(dim):
		for k in range(dim+1):
			for i in range(len(q)):
				qt = list(q[i]) + [1]
				c[k][j] += qt[k] * p[i][j]

	# Fill an an empty (dim+1) x (dim+1) matrix
	Q = np.zeros((dim+1,dim+1))
	for qi in q:
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
	
	# Make a result object
	class AffineTransformation:
		"""
		Result object that represents the transformation from affine fitter.
		"""

		def transform_vector(self, vector):
			"""
			Transform a vector according to the affine map.
		
			:param numpy.ndarray vector: vector to be transformed.
	
			:rtype: numpy.ndarray
			"""
			res = np.zeros(dim)
			for j in range(dim):
				for i in range(dim):
					res[j] += vector[i] * M[i][j+dim+1]
				res[j] += M[dim][j+dim+1]
			return res
		
	return AffineTransformation()
