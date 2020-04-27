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
      The movement of the control points is basically the weight (or displacement)
      :math:`\\boldsymbol{\\mu}` we set in the *parameters file*.

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
try:
    import configparser as configparser
except ImportError:
    import ConfigParser as configparser
import os
import numpy as np
from scipy import special
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound

from pygem import FFD as OriginalFFD

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
        >>> ffd.read_parameters('tests/test_datasets/parameters_test_ffd_sphere.prm')
        # TODO
        >>> new = free_form.modified_mesh_points
    """

    def __call__(self, obj):
        """
        This method performs the deformation on the CAD file.
        """

        # Manage input
        if isinstance(obj, str): # if a input filename is passed
            shape = # extract topo_shape from filename
        elif isinstance(obj, TopoDS_Shape):
            shape = obj
        # Maybe do we need to handle also Compound?
        else:
            raise TypeError

        print("Modifying faces")

        compound_builder = BRep_Builder()
        compound = OCC.TopoDS.TopoDS_Compound()
        compound_builder.MakeCompound(compound)

        # cycle on the faces to get the control points
        # init some quantities
        faceCount = 0
        face_list = []
        control_point_position = [0]
        faces_explorer = TopExp_Explorer(color, TopAbs_FACE)
        mesh_points = np.zeros(shape=(0, 3))

        while faces_explorer.More():
            points_orig = []
            points_def = []
            # performing some conversions to get the right format (BSplineSurface)
            print("Processing face ",faceCount) 
            face = OCC.TopoDS.topods_Face(faces_explorer.Current())
            nurbs_converter = BRepBuilderAPI_NurbsConvert(face)
            nurbs_converter.Perform(face)
            nurbs_face = nurbs_converter.Shape()
            face_aux = OCC.TopoDS.topods_Face(nurbs_face)
            int_tools = IntTools_FClass2d(face_aux, 1e-3)
            brep_face = BRep_Tool.Surface(OCC.TopoDS.topods_Face(nurbs_face))
            old_brep_face = BRep_Tool.Surface(OCC.TopoDS.topods_Face(nurbs_face))
            
            bounds = 0.0
            bounds = brep_face.GetObject().Bounds()
            #print("Bounds: ", bounds)
            bspline_face = geomconvert_SurfaceToBSplineSurface(brep_face)

            # we will then add an amount of nodes that will grant us our prescribed resolution both along u and v
            uKnotsToAdd = 20;
            vKnotsToAdd = 20;
            print("Added U knots: ", uKnotsToAdd )
            print("Added V knots: ", vKnotsToAdd )
            for i in range(uKnotsToAdd):
                bspline_face.GetObject().InsertUKnot(bounds[0]+i*(bounds[1]-bounds[0])/uKnotsToAdd,1,1e-7)
            for i in range(vKnotsToAdd):
                bspline_face.GetObject().InsertVKnot(bounds[2]+i*(bounds[3]-bounds[2])/vKnotsToAdd,1,1e-7)

            # openCascade object
            occ_face = bspline_face.GetObject()

            bounds = 0.0
            bounds = occ_face.Bounds()
            u_min = bounds[0]
            u_max = bounds[1]
            v_min = bounds[2]
            v_max = bounds[3]
            center = occ_face.Value((u_min+u_max)/2.0,(v_min+v_max)/2.0)
            print("*Center: ",center.X(),center.Y(),center.Z())

            # extract the Control Points of each face
            n_poles_u = occ_face.NbUPoles()
            n_poles_v = occ_face.NbVPoles()
            control_polygon_coordinates = np.zeros(\
            shape=(n_poles_u * n_poles_v, 3))
            #print("Number of poles: ", n_poles_u * n_poles_v)
            # cycle over the poles to get their coordinates
            i = 0
            for pole_u_direction in range(n_poles_u):
                for pole_v_direction in range(n_poles_v):
                    control_point_coordinates = occ_face.Pole(\
                    pole_u_direction + 1, pole_v_direction + 1)
                    control_polygon_coordinates[i, :] = [control_point_coordinates.X(),\
                    control_point_coordinates.Y(),\
                    control_point_coordinates.Z()]
                    
                    #control_point = gp_Pnt(control_point_coordinates.X(),\
                    #control_point_coordinates.Y(),control_point_coordinates.Z())
                    #display.DisplayShape(BRepBuilderAPI_MakeVertex(control_point).Vertex(),update=True)
                    #print("Original: ",control_point.X()," ",control_point.Y()," ",control_point.Z())
                    i+=1
            
                    
            #orig_control_polygon_coordinates = deepcopy(control_polygon_coordinates)
            ## SURFACES PHASE #####################################################
            src_pts = control_polygon_coordinates
            new_pts = super().__call__(src_pts) # dont touch this line
            i = 0
            for pole_u_direction in range(n_poles_u):
                for pole_v_direction in range(n_poles_v):
                    control_point = gp_Pnt(control_polygon_coordinates[i,0],
                                           control_polygon_coordinates[i,1],
                                           control_polygon_coordinates[i,2])
                    occ_face.SetPole(pole_u_direction + 1, pole_v_direction + 1, control_point)
                    i += 1
            # through moving the control points, we now changed the SURFACE of the FACE we are processing
            # we now need to obtain the curves (actually, the WIRES) that define the bounds of the surface and TRIM the surface
            # with them, to obtain the new face
            
            faceFilling = BRepFill_Filling()
            # we start creating a face with the modified surface. we will cut this new face with all the wires
            # that the original face had
            brep = BRepBuilderAPI_MakeFace(occ_face.GetHandle(), tol).Face()
            #IGESControl_Controller_Init()
            #writer = IGESControl_Writer()
            #writer.AddShape(brep)
            #nomefile = "untrimmed_face"+str(faceCount)+".iges"
            #writer.Write(nomefile)
            
            # we here start looping on the wires of the original face
            # in this forst loop we do nothing but count the wires in this face
            # few faces have more than one wire: if they have, it is because they have holes, and
            # if this is the case one wire is the outer, and the others are the inner ones.
            # the loop will also tell us which wire is the outer one
            wire_count = 0
            wire_explorer = TopExp_Explorer(face_aux, TopAbs_WIRE)
            while wire_explorer.More():
                wire = OCC.TopoDS.topods_Wire(wire_explorer.Current())
                if (wire == breptools_OuterWire(face_aux)):
                    print("Wire", wire_count+1, "is outer wire")
                wire_count += 1
                wire_explorer.Next()
            print("This face has ",wire_count," wires")
            
            #we now start really looping on the wires
            wire_count = 0
            outer_wires = []
            inner_wires = []
            brep_face = BRep_Tool.Surface(brep)
            wire_explorer = TopExp_Explorer(face_aux, TopAbs_WIRE)
            while wire_explorer.More():
                wire = OCC.TopoDS.topods_Wire(wire_explorer.Current())
                wire_count += 1
                h_bspline_edge = GeomConvert_CompCurveToBSplineCurve()
                edge_explorer = TopExp_Explorer(wire, TopAbs_EDGE)
                edgesCount=0
                while edge_explorer.More():
                    # performing some conversions to get the right format (BSplineSurface)
                    #print("Edge in curve: ", edgesCount)
                    edge = OCC.TopoDS.topods_Edge(edge_explorer.Current())
                    if (BRep_Tool.Degenerated(edge) == False):
                        bspline_converter = BRepBuilderAPI_NurbsConvert(edge)
                        bspline_converter.Perform(edge)
                        bspline_tshape_edge = bspline_converter.Shape()
                        h_geom_edge, a, b = BRep_Tool_Curve(OCC.TopoDS.topods_Edge(bspline_tshape_edge))
                        this_bspline_edge = geomconvert_CurveToBSplineCurve(h_geom_edge)
                        bspline_geom_edge = this_bspline_edge.GetObject()
                        h_bspline_edge.Add(this_bspline_edge,1e-4)
                    edgesCount += 1
                    #print("Curve ", curve_count, "added edge ", edgesCount)
                    edge_explorer.Next()
                
                bspline_geom_hedge = h_bspline_edge.BSplineCurve()
                bspline_geom_edge = bspline_geom_hedge.GetObject()
                unified_edge = BRepBuilderAPI_MakeEdge(bspline_geom_hedge).Edge()
                #aa = bspline_geom_edge.FirstParameter()
                #bb = bspline_geom_edge.LastParameter()
                #print("Unif. Edge First Point:", bspline_geom_edge.Value(aa).X(),
                #bspline_geom_edge.Value(aa).Y(),
                #bspline_geom_edge.Value(aa).Z())
                #print("Unif. Edge Last Point:", bspline_geom_edge.Value(bb).X(),
                #bspline_geom_edge.Value(bb).Y(),
                #bspline_geom_edge.Value(bb).Z())
                #unified_edge = BRepBuilderAPI_MakeEdge(bspline_geom_hedge).Edge()
                #print("SAVING UNIFIED CURVE: ", curve_count)
                #IGESControl_Controller_Init()
                #writer = IGESControl_Writer()
                #writer.AddShape(unified_edge)
                #filename = nomedir+"unfied_curve_"+str(curve_count)+".iges"
                ##print(filename)
                #writer.Write(filename)


                firstParam = bspline_geom_edge.FirstParameter()
                lastParam = bspline_geom_edge.LastParameter()
                #print("First Parameter: ", firstParam)
                #print("Last Parameter: ", lastParam)
                knotsToAdd = 200
                for i in range(knotsToAdd):
                    bspline_geom_edge.InsertKnot(firstParam+i*(lastParam-firstParam)/knotsToAdd,1,1e-7)
                shapesList = TopTools_ListOfShape()
                # openCascade object
                occ_edge = bspline_geom_edge

                # extract the Control Points of each face
                n_poles = occ_edge.NbPoles()
                control_polygon_coordinates = np.zeros(\
                shape=(n_poles , 3))
                #print("Number of poles: ", n_poles)
                # cycle over the poles to get their coordinates. The idea here is to move poles
                # coordinates to deform the curves
                i = 0
                for pole in range(n_poles):
                    control_point_coordinates = occ_edge.Pole(pole + 1)
                    control_polygon_coordinates[i, :] = [control_point_coordinates.X(),\
                    control_point_coordinates.Y(),\
                    control_point_coordinates.Z()]
                    #control_point = gp_Pnt(control_point_coordinates.X(),\
                    #control_point_coordinates.Y(),control_point_coordinates.Z())
                    #orig_control_point = control_point
                    #display.DisplayShape(BRepBuilderAPI_MakeVertex(control_point).Vertex(),update=True)
                    #print("Original: ",control_point.X()," ",control_point.Y()," ",control_point.Z())

                ## CURVES PHASE #######################################################
                src_pts = control_point_coordinates

                new_pts = super().__call__(src_pts) # dont touch this line

                # save here the `new_pts` into the shape
                ## END CURVES #########################################################

                i = 0
                for pole in range(n_poles):
                    #print("Function test:", mod_test_control_point.X(),mod_test_control_point.Y(),mod_test_control_point.Z())
                    #print("**")
                    control_point = gp_Pnt(control_polygon_coordinates[i,0],
                                           control_polygon_coordinates[i,1],
                                           control_polygon_coordinates[i,2])
                    occ_edge.SetPole(pole + 1, control_point)
                    i += 1
                        
                #occ_new_hedge = geomprojlib.Project(occ_edge.GetHandle(),occ_face.GetHandle())
                #new_curve = occ_new_hedge.GetObject()
                #print("??", new_curve.IsClosed())
                modified_edge= BRepBuilderAPI_MakeEdge(occ_edge.GetHandle()).Edge()
                shapesList.Append(modified_edge)
                #if new_curve.IsClosed()==False:
               #fixer_edge = BRepBuilderAPI_MakeEdge(new_curve.Value(new_curve.LastParameter()),new_curve.Value(new_curve.FirstParameter())).Edge()
                    #shapesList.Append(fixer_edge)
                #display.DisplayShape(modified_edge,update=True,color="BLUE1")
                
                wire_maker = BRepBuilderAPI_MakeWire()
                wire_maker.Add(shapesList)
                result_wire = wire_maker.Wire()


                # now, the wire can be oute or inner. we store the outer and (possible) inner ones in different lists
                # this is because we first need to trim the surface using the outer wire, and then we can trim it
                # with the wires corresponding to al the holse. if this is not done, the procedure will not work
                if (wire == breptools_OuterWire(face_aux)):
                    outer_wires.append(result_wire)
                else:
                    inner_wires.append(result_wire)
                wire_count += 1
                wire_explorer.Next()

            #faceFilling.LoadInitSurface(brep)
            #faceFilling.Build()
            # checking if things worked out
            #print("Is face Filling Done? ",faceFilling.IsDone())

            # so once we finished looping on all the wires to modify them, we use the only outer one to trim the surface
            face_maker = BRepBuilderAPI_MakeFace(occ_face.GetHandle(),outer_wires[0])

            # and then add all other inner wires for the holes
            for inner_wire in inner_wires:
                face_maker.Add(inner_wire)
                            
            # finally, we get our trimmed face with all its holes
            brep_surf = face_maker.Face()
            #brep_surf = faceFilling.Face()
            #IGESControl_Controller_Init()
            #writer = IGESControl_Writer()
            #writer.AddShape(brep_surf)
            #nomefile = nomedir+"trimmed_face"+str(faceCount)+".iges"
            #writer.Write(nomefile)
            
            
            # and move to the next face
            faceCount += 1
            faces_explorer.Next()
            


        ## END SURFACES #######################################################






        if isinstance(obj, str): # if a input filename is passed
            # save the shape exactly to the filename, aka `obj`
        elif isinstance(obj, TopoDS_Shape):
            return shape
