"""
Module focused on the implementation of the Radial Basis Functions interpolation
technique.  This technique is still based on the use of a set of parameters, the
so-called control points, as for FFD, but RBF is interpolatory. Another
important key point of RBF strategy relies in the way we can locate the control
points: in fact, instead of FFD where control points need to be placed inside a
regular lattice, with RBF we hano no more limitations. So we have the
possibility to perform localized control points refiniments.
The module is analogous to the freeform one.

:Theoretical Insight:

    As reference please consult M.D. Buhmann, Radial Basis Functions, volume 12
    of Cambridge monographs on applied and computational mathematics. Cambridge
    University Press, UK, 2003.  This implementation follows D. Forti and G.
    Rozza, Efficient geometrical parametrization techniques of interfaces for
    reduced order modelling: application to fluid-structure interaction coupling
    problems, International Journal of Computational Fluid Dynamics.

    RBF shape parametrization technique is based on the definition of a map,
    :math:`\\mathcal{M}(\\boldsymbol{x}) : \\mathbb{R}^n \\rightarrow
    \\mathbb{R}^n`, that allows the possibility of transferring data across
    non-matching grids and facing the dynamic mesh handling. The map introduced
    is defines as follows

    .. math::
        \\mathcal{M}(\\boldsymbol{x}) = p(\\boldsymbol{x}) + 
        \\sum_{i=1}^{\\mathcal{N}_C} \\gamma_i
        \\varphi(\\| \\boldsymbol{x} - \\boldsymbol{x_{C_i}} \\|)

    where :math:`p(\\boldsymbol{x})` is a low_degree polynomial term,
    :math:`\\gamma_i` is the weight, corresponding to the a-priori selected
    :math:`\\mathcal{N}_C` control points, associated to the :math:`i`-th basis
    function, and :math:`\\varphi(\\| \\boldsymbol{x} - \\boldsymbol{x_{C_i}}
    \\|)` a radial function based on the Euclidean distance between the control
    points position :math:`\\boldsymbol{x_{C_i}}` and :math:`\\boldsymbol{x}`.
    A radial basis function, generally, is a real-valued function whose value
    depends only on the distance from the origin, so that
    :math:`\\varphi(\\boldsymbol{x}) = \\tilde{\\varphi}(\\| \\boldsymbol{x}
    \\|)`.

    The matrix version of the formula above is:

    .. math::
        \\mathcal{M}(\\boldsymbol{x}) = \\boldsymbol{c} +
        \\boldsymbol{Q}\\boldsymbol{x} +
        \\boldsymbol{W^T}\\boldsymbol{d}(\\boldsymbol{x})

    The idea is that after the computation of the weights and the polynomial
    terms from the coordinates of the control points before and after the
    deformation, we can deform all the points of the mesh accordingly.  Among
    the most common used radial basis functions for modelling 2D and 3D shapes,
    we consider Gaussian splines, Multi-quadratic biharmonic splines, Inverted
    multi-quadratic biharmonic splines, Thin-plate splines, Beckert and
    Wendland :math:`C^2` basis and Polyharmonic splines all defined and
    implemented below.
"""

import numpy as np
from pygem import RBF as OriginalRBF
from .cad_deformation import CADDeformation

class RBF(CADDeformation, OriginalRBF):
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

        >>> from pygem.cad import RBF
        >>> rbf = RBF()
        >>> rbf.read_parameters(
        >>>        'tests/test_datasets/parameters_test_ffd_iges.prm')
        >>> input_cad_file_name = "input.iges"
        >>> modified_cad_file_name = "output.iges"
        >>> rbf(input_cad_file_name, modified_cad_file_name)
    """
    def __init__(self,
                 original_control_points=None,
                 deformed_control_points=None,
                 func='gaussian_spline',
                 radius=0.5,
                 extra_parameter=None,
                 u_knots_to_add=30,
                 v_knots_to_add=30,
                 t_knots_to_add=30,
                 tolerance=1e-4):
        OriginalRBF.__init__(self, 
                             original_control_points=original_control_points, 
                             deformed_control_points=deformed_control_points,
                             func=func,
                             radius=radius,
                             extra_parameter=extra_parameter)
        CADDeformation.__init__(self,
                                u_knots_to_add=u_knots_to_add,
                                v_knots_to_add=v_knots_to_add,
                                t_knots_to_add=t_knots_to_add,
                                tolerance=tolerance)
