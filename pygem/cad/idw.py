"""
Module focused on the Inverse Distance Weighting interpolation technique.
The IDW algorithm is an average moving interpolation that is usually applied to
highly variable data. The main idea of this interpolation strategy lies in
fact that it is not desirable to honour local high/low values but rather to look
at a moving average of nearby data points and estimate the local trends.
The node value is calculated by averaging the weighted sum of all the points.
Data points that lie progressively farther from the node inuence much less the
computed value than those lying closer to the node.

:Theoretical Insight:

    This implementation is based on the simplest form of inverse distance
    weighting interpolation, proposed by D. Shepard, A two-dimensional
    interpolation function for irregularly-spaced data, Proceedings of the 23 rd
    ACM National Conference.

    The interpolation value :math:`u` of a given point :math:`\\mathrm{x}`
    from a set of samples :math:`u_k = u(\\mathrm{x}_k)`, with
    :math:`k = 1,2,\\dotsc,\\mathcal{N}`, is given by:

    .. math::
        u(\\mathrm{x}) = \\displaystyle\\sum_{k=1}^\\mathcal{N}
        \\frac{w(\\mathrm{x},\\mathrm{x}_k)}
        {\\displaystyle\\sum_{j=1}^\\mathcal{N} w(\\mathrm{x},\\mathrm{x}_j)}
        u_k

    where, in general, :math:`w(\\mathrm{x}, \\mathrm{x}_i)` represents the
    weighting function:

    .. math::
        w(\\mathrm{x}, \\mathrm{x}_i) = \\| \\mathrm{x} - \\mathrm{x}_i \\|^{-p}

    being :math:`\\| \\mathrm{x} - \\mathrm{x}_i \\|^{-p} \\ge 0` is the
    Euclidean distance between :math:`\\mathrm{x}` and data point
    :math:`\\mathrm{x}_i` and :math:`p` is a power parameter, typically equal to
    2.
"""

import numpy as np
from pygem import IDW as OriginalIDW
from .cad_deformation import CADDeformation

class IDW(CADDeformation, OriginalIDW):
    """
    Class that handles the Free Form Deformation on the mesh points.

    :param FFDParameters ffd_parameters: parameters of the Free Form
        Deformation.
    :param numpy.ndarray original_mesh_points: coordinates of the original
        points of the mesh.

    :param list n_control_points: number of control points in the x, y, and z
        direction. If not provided it is set to [2, 2, 2].

    :cvar numpy.ndarray box_length: dimension of the FFD bounding box, in the
        x, y and z direction (local coordinate system).
    :cvar numpy.ndarray box_origin: the x, y and z coordinates of the origin of
        the FFD bounding box.
    :cvar numpy.ndarray rot_angle: rotation angle around x, y and z axis of the
        FFD bounding box.
    :cvar numpy.ndarray n_control_points: the number of control points in the
        x, y, and z direction.
    :cvar numpy.ndarray array_mu_x: collects the displacements (weights) along
        x, normalized with the box length x.
    :cvar numpy.ndarray array_mu_y: collects the displacements (weights) along
        y, normalized with the box length y.
    :cvar numpy.ndarray array_mu_z: collects the displacements (weights) along
        z, normalized with the box length z.

    :Example:

        >>> from pygem.cad import IDW
        >>> idw = IDW()
        >>> idw.read_parameters(
        >>>        'tests/test_datasets/parameters_test_idw_iges.prm')
        >>> input_cad_file_name = "input.iges"
        >>> modified_cad_file_name = "output.iges"
        >>> idw(input_cad_file_name, modified_cad_file_name)
    """
    def __init__(self,
                 original_control_points=None,
                 deformed_control_points=None,
                 power=2,
                 u_knots_to_add=30,
                 v_knots_to_add=30,
                 t_knots_to_add=30,
                 tolerance=1e-4):
        OriginalIDW.__init__(self,
                             original_control_points=original_control_points,
                             deformed_control_points=deformed_control_points,
                             power=power)
        CADDeformation.__init__(self,
                                u_knots_to_add=u_knots_to_add,
                                v_knots_to_add=v_knots_to_add,
                                t_knots_to_add=t_knots_to_add,
                                tolerance=tolerance)
