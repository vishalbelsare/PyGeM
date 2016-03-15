"""
Utilities for the affine transformations of the bounding box of the Free Form Deformation.
"""

import math
import numpy as np

def angles2matrix(z=0, y=0, x=0):
    """
    Returns matrix for given rotations around z, y and x axes. The output rotation matrix is equal to the composition of the
    individual rotations. Rotations are counter-clockwise.

    :param float z: rotation angle (in radians) around z-axis.
    :param float y: rotation angle (in radians) around y-axis.
    :param float x: rotation angle (in radians) around x-axis.

    :return: rot_matrix: rotation matrix for the given angles. The matrix shape is always (3,3).
    :rtype: float numpy array.
    
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
