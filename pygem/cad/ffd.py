"""
Utilities for performing Free Form Deformation (FFD)

:Theoretical Insight:

    Free Form Deformation is a technique for the efficient, smooth and accurate
    geometrical parametrization. It has been proposed the first time in
    *Sederberg, Thomas W., and Scott R. Parry. "Free-form deformation of solid
    geometric models." ACM SIGGRAPH computer graphics 20.4 (1986): 151-160*. It
    consists in three different step:

    - Mapping the physical domain to the reference one with map
      :math:`\\boldsymbol{\\psi}`.  In the code it is named *transformation*.

    - Moving some control points to deform the lattice with :math:`\\hat{T}`.
      The movement of the control points is basically the weight (or
      displacement) :math:`\\boldsymbol{\\mu}` we set in the *parameters file*.

    - Mapping back to the physical domain with map
      :math:`\\boldsymbol{\\psi}^{-1}`.  In the code it is named
      *inverse_transformation*.

    FFD map (:math:`T`) is the composition of the three maps, that is

    .. math:: T(\\cdot, \\boldsymbol{\\mu}) = (\\Psi^{-1} \\circ \\hat{T} \\circ
            \\Psi) (\\cdot, \\boldsymbol{\\mu})

    In this way, every point inside the FFD box changes it position according to

    .. math:: \\boldsymbol{P} = \\boldsymbol{\\psi}^{-1} \\left( \\sum_{l=0}^L
            \\sum_{m=0}^M \\sum_{n=0}^N
            \\mathsf{b}_{lmn}(\\boldsymbol{\\psi}(\\boldsymbol{P}_0))
            \\boldsymbol{\\mu}_{lmn} \\right)

    where :math:`\\mathsf{b}_{lmn}` are Bernstein polynomials.  We improve the
    traditional version by allowing a rotation of the FFD lattice in order to
    give more flexibility to the tool.

    You can try to add more shapes to the lattice to allow more and more
    involved transformations.

"""

import numpy as np
from OCC.Core.TopoDS import (TopoDS_Shape, topods_Wire,     \
                             TopoDS_Compound, topods_Face,  \
                             topods_Edge, TopoDS_Face,      \
                             TopoDS_Wire)
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_WIRE)
from OCC.Core.TopTools import TopTools_ListOfShape
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeFace,    \
                                     BRepBuilderAPI_MakeWire,    \
                                     BRepBuilderAPI_MakeEdge,    \
                                     BRepBuilderAPI_NurbsConvert)
from OCC.Core.BRep import BRep_Tool, BRep_Tool_Curve
from OCC.Core.Geom import Geom_BSplineCurve, Geom_BSplineSurface
from OCC.Core.GeomConvert import (geomconvert_SurfaceToBSplineSurface, \
                                  geomconvert_CurveToBSplineCurve,     \
                                  GeomConvert_CompCurveToBSplineCurve)
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepTools import breptools_OuterWire

from pygem import FFD as OriginalFFD
from pygem.cad.igeshandler import IgesHandler

class FFD(OriginalFFD):
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

        >>> import pygem.freeform as ffd
        >>> import pygem.params as ffdp
        >>> import numpy as np
        >>> ffd = FFD()
        >>> ffd.read_parameters(
        >>>        'tests/test_datasets/parameters_test_ffd_iges.prm')
        >>> input_cad_file_name = "input.iges"
        >>> modified_cad_file_name = "output.iges"
        >>> ffd(input_cad_file_name,modified_cad_file_name)
    """
    def __init__(self,
                 n_control_points=None,
                 u_knots_to_add=30,
                 v_knots_to_add=30,
                 t_knots_to_add=30,
                 tolerance=1e-4):
        super().__init__(n_control_points=None)
        self.u_knots_to_add = u_knots_to_add
        self.v_knots_to_add = v_knots_to_add
        self.t_knots_to_add = t_knots_to_add
        self.tolerance = tolerance

    def _bspline_surface_from_face(self, face):
        """
        Takes a TopoDS_Face and transforms it into a Bspline_Surface.
        :param TopoDS_Face face: the face to be converted.
        :return: the output Bspline surface
        :rtype:  Geom_BSplineSurface
        """
        if not isinstance(face, TopoDS_Face):
            raise TypeError("face must be a TopoDS_Face")
        # TopoDS_Face converted to Nurbs
        nurbs_face = topods_Face(BRepBuilderAPI_NurbsConvert(face).Shape())
        # GeomSurface obtained from Nurbs face
        surface = BRep_Tool.Surface(nurbs_face)
        # surface is now further converted to a bspline surface
        bspline_surface = geomconvert_SurfaceToBSplineSurface(surface)
        return bspline_surface


    def _bspline_curve_from_wire(self, wire):
        """
        Takes a TopoDS_Wire and transforms it into a Bspline_Curve.
        :param TopoDS_Wire wire: the wire to be converted.
        :return: The output Bspline curve
        :rtype:  Geom_BSplineCurve
        """
        if not isinstance(wire, TopoDS_Wire):
            raise TypeError("face must be a TopoDS_Wire")
        # joining all the wire edges in a single curve here
        # composite curve builder (can only join Bspline curves)
        composite_curve_builder = GeomConvert_CompCurveToBSplineCurve()
        # iterator to edges in the TopoDS_Wire
        edge_explorer = TopExp_Explorer(wire, TopAbs_EDGE)
        while edge_explorer.More():
            # getting the edge from the iterator
            edge = topods_Edge(edge_explorer.Current())
            # edge can be joined only if it is not degenerated (zero
            # length)
            if not BRep_Tool.Degenerated(edge):
                # the edge must be converted to Nurbs edge
                # the Nurbs edge converter class
                nurbs_converter = BRepBuilderAPI_NurbsConvert(edge)
                nurbs_converter.Perform(edge)
                # the Nurbs edge
                nurbs_edge = topods_Edge(nurbs_converter.Shape())
                # here we extract the underlying curve from the Nurbs
                # edge
                nurbs_curve, _, _ = BRep_Tool_Curve(nurbs_edge)
                # we convert the Nurbs curve to Bspline curve
                bspline_curve = \
                    geomconvert_CurveToBSplineCurve(nurbs_curve)
                # we can now add the Bspline curve to
                # the composite wire curve
                composite_curve_builder.Add(bspline_curve,
                                            self.tolerance)
            edge_explorer.Next()

        # GeomCurve obtained by the builder after edges are joined
        comp_curve = composite_curve_builder.BSplineCurve()
        return comp_curve

    def _enrich_curve_knots(self, bsp_curve):
        """
        Takes a Geom_BSplineCurve and adds self.t_knots_to_add
        poles to it.
        :param Geom_BSplineCurve bsp_curve: Bspline curve to be enriched.
        """
        if not isinstance(bsp_curve, Geom_BSplineCurve):
            raise TypeError("bsp_curve must be a Geom_BSplineCurve")
        # number of knots is enriched here, if required, to
        # enhance precision
        # start parameter of composite curve
        first_param = bsp_curve.FirstParameter()
        # end parameter of composite curve
        last_param = bsp_curve.LastParameter()
        for i in range(self.t_knots_to_add):
            bsp_curve.InsertKnot(first_param+ \
                           i*(last_param-first_param)/self.t_knots_to_add, 1, \
                           self.tolerance)


    def _enrich_surface_knots(self, bsp_surface):
        """
        Takes a Geom_Bspline_Surface and adds self.u_knots_to_add
        and self.v_knots_to_add knots to it in u and v direction respectively.
        :param Geom_BSplineSurface bsp_surface: Bspline curve to be enriched.
        """
        if not isinstance(bsp_surface, Geom_BSplineSurface):
            raise TypeError("bsp_surface must be a Geom_BSplineSurface")
        # we will add the prescribed amount of nodes
        # both along u and v parametric directions

        # bounds (in surface parametric space) of the surface
        bounds = bsp_surface.Bounds()

        for i in range(self.u_knots_to_add):
            bsp_surface.InsertUKnot(bounds[0]+ \
                i*(bounds[1]-bounds[0])/self.u_knots_to_add, 1, self.tolerance)
        for i in range(self.v_knots_to_add):
            bsp_surface.InsertVKnot(bounds[2]+ \
                i*(bounds[3]-bounds[2])/self.v_knots_to_add, 1, self.tolerance)

    def _deform_bspline_curve(self, bsp_curve):
        """
        Takes a Geom_BSplineCurve and deforms it through FFD.
        :param Geom_BSplineCurve bsp_curve: Bspline curve to be deformed.
        """
        if not isinstance(bsp_curve, Geom_BSplineCurve):
            raise TypeError("bsp_curve must be a Geom_BSplineCurve")
        # we first extract the poles of the curve
        # poles number
        n_poles = bsp_curve.NbPoles()
        # array containing the poles coordinates
        poles_coordinates = np.zeros(shape=(n_poles, 3))
        # cycle over the poles to get their coordinates

        for pole_id in range(n_poles):
            # gp_Pnt corresponding to the pole
            pole = bsp_curve.Pole(pole_id + 1)
            # coordinates are added to array
            poles_coordinates[pole_id, :] = [pole.X(), pole.Y(), pole.Z()]


        # the new poles positions are computed through FFD
        new_pts = super().__call__(poles_coordinates)

        # the Bspline curve is now looped again to
        # set the poles positions  to new_points
        i = 0
        for pole in range(n_poles):
            # gp_Point corresponding to the new pole coordinates
            control_point = gp_Pnt(new_pts[i, 0],
                                   new_pts[i, 1],
                                   new_pts[i, 2])
            bsp_curve.SetPole(pole + 1, control_point)
            i += 1

    def _deform_bspline_surface(self, bsp_surface):
        """
        Takes a Geom_BSplineSurface and deforms it through FFD.
        :param Geom_BSplineSurface bsp_surface: Bspline curve to be deformed.
        """
        if not isinstance(bsp_surface, Geom_BSplineSurface):
            raise TypeError("bsp_surface must be a Geom_BSplineSurface")
        # we first extract the poles of the curve
        # number of poles in u direction
        n_poles_u = bsp_surface.NbUPoles()
        # number of poles in v direction
        n_poles_v = bsp_surface.NbVPoles()
        # array which will contain the coordinates of the poles
        poles_coordinates = np.zeros(shape=(n_poles_u * n_poles_v, 3))
        # cycle over the poles to get their coordinates
        i = 0
        for pole_u_direction in range(n_poles_u):
            for pole_v_direction in range(n_poles_v):
                # gp_Pnt containing the current pole
                pole = bsp_surface.Pole(pole_u_direction + 1,
                                        pole_v_direction + 1)
                poles_coordinates[i, :] = [pole.X(), pole.Y(), pole.Z()]
                i += 1

        # the new poles positions are computed through FFD
        new_pts = super().__call__(poles_coordinates)

        # the surface is now looped again to
        # set the new poles positions
        i = 0
        for pole_u_direction in range(n_poles_u):
            for pole_v_direction in range(n_poles_v):
                new_pole = gp_Pnt(new_pts[i, 0],
                                  new_pts[i, 1],
                                  new_pts[i, 2])
                bsp_surface.SetPole(pole_u_direction + 1,
                                    pole_v_direction + 1,
                                    new_pole)
                i += 1


    def __call__(self, obj, dst=None):
        """
        This method performs the deformation on the CAD file.
        """

        # Manage input
        if isinstance(obj, str): # if a input filename is passed
            iges_handler = IgesHandler()
            shape = iges_handler.load_shape_from_file(obj)
        elif isinstance(obj, TopoDS_Shape):
            shape = obj
        # Maybe do we need to handle also Compound?
        else:
            raise TypeError


        #create compound to store modified faces
        compound_builder = BRep_Builder()
        compound = TopoDS_Compound()
        compound_builder.MakeCompound(compound)


        # cycle on the faces to get the control points

        # iterator to faces (TopoDS_Shape) contained in the shape
        faces_explorer = TopExp_Explorer(shape, TopAbs_FACE)

        while faces_explorer.More():
            # performing some conversions to get the right
            # format (BSplineSurface)
            # TopoDS_Face obtained from iterator
            face = topods_Face(faces_explorer.Current())
            # performing some conversions to get the right
            # format (BSplineSurface)
            bspline_surface = self._bspline_surface_from_face(face)

            # add the required amount of poles in u and v directions
            self._enrich_surface_knots(bspline_surface)

            # deform the Bspline surface through FFD
            self._deform_bspline_surface(bspline_surface)

            # through moving the control points, we now changed the SURFACE
            # underlying FACE we are processing. we now need to obtain the
            # curves (actually, the WIRES) that define the bounds of the
            # surface and TRIM the surface with them, to obtain the new FACE


            #we now start really looping on the wires
            #we will create a single curve joining all the edges in the wire
            # the curve must be a bspline curve so we need to make conversions
            # through all the way

            # list that will contain the (single) outer wire of the face
            outer_wires = []
            # list that will contain all the inner wires (holes) of the face
            inner_wires = []
            # iterator to loop over TopoDS_Wire in the original (undeformed)
            # face
            wire_explorer = TopExp_Explorer(face, TopAbs_WIRE)
            while wire_explorer.More():
                # wire obtained from the iterator
                wire = topods_Wire(wire_explorer.Current())

                # getting a bpline curve joining all the edges of the wire
                composite_curve = self._bspline_curve_from_wire(wire)

                # adding all the required knots to the Bspline curve
                self._enrich_curve_knots(composite_curve)

                # deforming the Bspline curve through FFD
                self._deform_bspline_curve(composite_curve)

                # the GeomCurve corresponding to the whole edge has now
                # been deformed. Now we must make it become an proper
                # wire

                # list of shapes (needed by the wire generator)
                shapes_list = TopTools_ListOfShape()

                # edge (to be converted to wire) obtained from the modified
                # Bspline curve
                modified_composite_edge = \
                    BRepBuilderAPI_MakeEdge(composite_curve).Edge()
                # modified edge is added to shapes_list
                shapes_list.Append(modified_composite_edge)

                # wire builder
                wire_maker = BRepBuilderAPI_MakeWire()
                wire_maker.Add(shapes_list)
                # deformed wire is finally obtained
                modified_wire = wire_maker.Wire()

                # now, the wire can be outer or inner. we store the outer
                # and (possible) inner ones in different lists
                # this is because we first need to trim the surface
                # using the outer wire, and then we can trim it
                # with the wires corresponding to all the holes.
                # the wire order is important, in the trimming process
                if wire == breptools_OuterWire(face):
                    outer_wires.append(modified_wire)
                else:
                    inner_wires.append(modified_wire)
                wire_explorer.Next()


            # so once we finished looping on all the wires to modify them,
            # we first use the only outer one to trim the surface
            # face builder object
            face_maker = BRepBuilderAPI_MakeFace(bspline_surface,
                                                 outer_wires[0])

            # and then add all other inner wires for the holes
            for inner_wire in inner_wires:
                face_maker.Add(inner_wire)


            # finally, we get our trimmed face with all its holes
            trimmed_modified_face = face_maker.Face()

            # trimmed_modified_face is added to the modified faces compound
            compound_builder.Add(compound, trimmed_modified_face)

            # and move to the next face
            faces_explorer.Next()



        ## END SURFACES #################################################


        if isinstance(dst, str): # if a input filename is passed
            # save the shape exactly to the filename, aka `dst`
            iges_handler = IgesHandler()
            iges_handler.write_shape_to_file(compound, dst)
        else:
            return compound
