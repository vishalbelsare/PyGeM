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
import os
import numpy as np

from scipy.spatial.distance import cdist

from .rbf_factory import RBFFactory


class RBF(object):
    """
    Class that handles the Radial Basis Functions interpolation on the mesh
    points.

    :cvar numpy.matrix weights: the matrix formed by the weights corresponding
        to the a-priori selected N control points, associated to the basis
        functions and c and Q terms that describe the polynomial of order one
        p(x) = c + Qx.  The shape is (n_control_points+1+3)-by-3. It is computed
        internally.
    :cvar string basis: name of the basis functions to use in the
        transformation. The functions implemented so far are: gaussian spline,
        multi quadratic biharmonic spline, inv multi quadratic biharmonic
        spline, thin plate spline, beckert wendland c2 basis, polyharmonic
        splines. For a comprehensive list with details see the class
        :class:`~pygem.radialbasis.RBF`. The default value is 'gaussian_spline'.
    :cvar float radius: the scaling parameter r that affects the shape of the
        basis functions.  For details see the class
        :class:`~pygem.radialbasis.RBF`. The default value is 0.5.
    :Example:

        >>> from pygem import RBF
        >>> import numpy as np
        >>> rbf = RBF('gaussian_spline')
        >>> xv = np.linspace(0, 1, 20)
        >>> yv = np.linspace(0, 1, 20)
        >>> zv = np.linspace(0, 1, 20)
        >>> z, y, x = np.meshgrid(zv, yv, xv)
        >>> mesh = np.array([x.ravel(), y.ravel(), z.ravel()])
        >>> deformed_mesh = rbf(mesh)
    """

    def __init__(self, 
                 original_control_points=None, 
                 deformed_control_points=None,
                 func='gaussian_spline',
                 radius=0.5,
                 extra_parameter=None):

        if callable(func):
            self.basis = func
        elif isinstance(func, str):
            self.basis = RBFFactory(func)
        else:
            raise TypeError('`func` is not valid.')

        self.radius = radius

        if original_control_points is None:
            self.original_control_points = np.array(
                [[0., 0., 0.], [0., 0., 1.], [0., 1., 0.], [1., 0., 0.],
                 [0., 1., 1.], [1., 0., 1.], [1., 1., 0.], [1., 1., 1.]])
        else:
            self.original_control_points = original_control_points

        if deformed_control_points is None:
            self.deformed_control_points = np.array(
                [[0., 0., 0.], [0., 0., 1.], [0., 1., 0.], [1., 0., 0.],
                 [0., 1., 1.], [1., 0., 1.], [1., 1., 0.], [1., 1., 1.]])
        else:
            self.deformed_control_points = deformed_control_points

        self.extra = extra_parameter if extra_parameter else dict()

        self.weights = self._get_weights(
            self.original_control_points,
            self.deformed_control_points)

    @property
    def n_control_points(self):
        """
        Total number of control points.

        :rtype: int
        """
        return self.original_control_points.shape[0]

    def _get_weights(self, X, Y):
        """
        This private method, given the original control points and the deformed
        ones, returns the matrix with the weights and the polynomial terms, that
        is :math:`W`, :math:`c^T` and :math:`Q^T`. The shape is
        (n_control_points+1+3)-by-3.

        :param numpy.ndarray X: it is an n_control_points-by-3 array with the
            coordinates of the original interpolation control points before the
            deformation.
        :param numpy.ndarray Y: it is an n_control_points-by-3 array with the
            coordinates of the interpolation control points after the
            deformation.

        :return: weights: the matrix with the weights and the polynomial terms.
        :rtype: numpy.matrix
        """
        npts, dim = X.shape
        H = np.zeros((npts + 3 + 1, npts + 3 + 1))
        H[:npts, :npts] = self.basis(cdist(X, X), self.radius, **self.extra)
        H[npts, :npts] = 1.0
        H[:npts, npts] = 1.0
        H[:npts, -3:] = X
        H[-3:, :npts] = X.T

        rhs = np.zeros((npts + 3 + 1, dim))
        rhs[:npts, :] = Y
        weights = np.linalg.solve(H, rhs)
        return weights

    def read_parameters(self, filename='parameters_rbf.prm'):
        """
        Reads in the parameters file and fill the self structure.

        :param string filename: parameters file to be read in. Default value is
            parameters_rbf.prm.
        """
        if not isinstance(filename, str):
            raise TypeError('filename must be a string')

        # Checks if the parameters file exists. If not it writes the default
        # class into filename.  It consists in the vetices of a cube of side one
        # with a vertex in (0, 0, 0) and opposite one in (1, 1, 1).
        if not os.path.isfile(filename):
            self.write_parameters(filename)
            return

        config = configparser.RawConfigParser()
        config.read(filename)

        self.basis = config.get('Radial Basis Functions', 'basis function')
        self.radius = config.getfloat('Radial Basis Functions', 'radius')

        ctrl_points = config.get('Control points', 'original control points')
        lines = ctrl_points.split('\n')
        self.original_control_points = np.zeros((len(lines), 3))
        for line, i in zip(lines, list(range(0, self.n_control_points))):
            values = line.split()
            self.original_control_points[i] = np.array(
                [float(values[0]),
                 float(values[1]),
                 float(values[2])])

        mod_points = config.get('Control points', 'deformed control points')
        lines = mod_points.split('\n')

        if len(lines) != self.n_control_points:
            raise TypeError("The number of control points must be equal both in"
                            "the 'original control points' and in the 'deformed"
                            "control points' section of the parameters file"
                            "({0!s})".format(filename))

        self.deformed_control_points = np.zeros((self.n_control_points, 3))
        for line, i in zip(lines, list(range(0, self.n_control_points))):
            values = line.split()
            self.deformed_control_points[i] = np.array(
                [float(values[0]),
                 float(values[1]),
                 float(values[2])])

    def write_parameters(self, filename='parameters_rbf.prm'):
        """
        This method writes a parameters file (.prm) called `filename` and fills
        it with all the parameters class members. Default value is
        parameters_rbf.prm.

        :param string filename: parameters file to be written out.
        """
        if not isinstance(filename, str):
            raise TypeError("filename must be a string")

        output_string = ""
        output_string += '\n[Radial Basis Functions]\n'
        output_string += '# This section describes the radial basis functions'
        output_string += ' shape.\n'

        output_string += '\n# basis funtion is the name of the basis functions'
        output_string += ' to use in the transformation. The functions\n'
        output_string += '# implemented so far are: gaussian_spline,'
        output_string += ' multi_quadratic_biharmonic_spline,\n'
        output_string += '# inv_multi_quadratic_biharmonic_spline,'
        output_string += ' thin_plate_spline, beckert_wendland_c2_basis,'
        output_string += ' polyharmonic_spline.\n'
        output_string += '# For a comprehensive list with details see the'
        output_string += ' class RBF.\n'
        output_string += 'basis function: {}\n'.format(str(self.basis))

        output_string += '\n# radius is the scaling parameter r that affects'
        output_string += ' the shape of the basis functions. See the'
        output_string += ' documentation\n'
        output_string += '# of the class RBF for details.\n'
        output_string += 'radius: {}\n'.format(str(self.radius))

        output_string += '\n\n[Control points]\n'
        output_string += '# This section describes the RBF control points.\n'

        output_string += '\n# original control points collects the coordinates'
        output_string += ' of the interpolation control points before the'
        output_string += ' deformation.\n'

        output_string += 'original control points:'
        offset = 1
        for i in range(0, self.n_control_points):
            output_string += offset * ' ' + str(
                self.original_control_points[i][0]) + '   ' + str(
                    self.original_control_points[i][1]) + '   ' + str(
                        self.original_control_points[i][2]) + '\n'
            offset = 25

        output_string += '\n# deformed control points collects the coordinates'
        output_string += ' of the interpolation control points after the'
        output_string += ' deformation.\n'

        output_string += 'deformed control points:'
        offset = 1
        for i in range(0, self.n_control_points):
            output_string += offset * ' ' + str(
                self.deformed_control_points[i][0]) + '   ' + str(
                    self.deformed_control_points[i][1]) + '   ' + str(
                        self.deformed_control_points[i][2]) + '\n'
            offset = 25

        with open(filename, 'w') as f:
            f.write(output_string)

    def __str__(self):
        """
        This method prints all the RBF parameters on the screen. Its purpose is
        for debugging.
        """
        string = ''
        string += 'basis function = {}\n'.format(self.basis)
        string += 'radius = {}\n'.format(self.radius)
        string += '\noriginal control points =\n'
        string += '{}\n'.format(self.original_control_points)
        string += '\ndeformed control points =\n'
        string += '{}\n'.format(self.deformed_control_points)
        return string

    def plot_points(self, filename=None):
        """
        Method to plot the control points. It is possible to save the resulting
        figure.

        :param str filename: if None the figure is shown, otherwise it is saved
            on the specified `filename`. Default is None.
        """
        fig = plt.figure(1)
        axes = fig.add_subplot(111, projection='3d')
        orig = axes.scatter(
            self.original_control_points[:, 0],
            self.original_control_points[:, 1],
            self.original_control_points[:, 2],
            c='blue',
            marker='o')
        defor = axes.scatter(
            self.deformed_control_points[:, 0],
            self.deformed_control_points[:, 1],
            self.deformed_control_points[:, 2],
            c='red',
            marker='x')

        axes.set_xlabel('X axis')
        axes.set_ylabel('Y axis')
        axes.set_zlabel('Z axis')

        plt.legend(
            (orig, defor), ('Original', 'Deformed'),
            scatterpoints=1,
            loc='lower left',
            ncol=2,
            fontsize=10)

        # Show the plot to the screen
        if filename is None:
            plt.show()
        else:
            fig.savefig(filename)

    def __call__(self, src_pts):
        """
        This method performs the deformation of the mesh points. After the
        execution it sets `self.modified_mesh_points`.
        """
        H = np.zeros((n_mesh_points, self.n_control_points+3+1))
        H[:, :self.n_control_points] = self.basis(
            cdist(src_pts, self.original_control_points),
            self.radius, **self.extra)
        H[:, n_control_points] = 1.0
        H[:, -3:] = self.original_mesh_points
        self.modified_mesh_points = np.asarray(np.dot(H, self.weights))
