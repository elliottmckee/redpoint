import gmsh
import sys

'''
See notebook for picture definitions

'''


if __name__ == '__main__':


    xyz = (0.0, 0.0, 0.0)
    L = 1.0

    output_file_base = 'cube'


    gmsh.initialize()
    gmsh.model.add("model1")

    # Add Points
    # gmsh.model.geo.addPoint(0, 0, 0, meshsize=1, tag=1)
    gmsh.model.geo.addPoint(xyz[0],     xyz[1],     xyz[2], 0.5, tag=1)
    gmsh.model.geo.addPoint(xyz[0]+L,   xyz[1],     xyz[2], 0.5, tag=2)
    gmsh.model.geo.addPoint(xyz[0]+L,   xyz[1]+L,   xyz[2], 0.5, tag=3)
    gmsh.model.geo.addPoint(xyz[0],     xyz[1]+L,   xyz[2], 0.5, tag=4)
    gmsh.model.geo.addPoint(xyz[0],     xyz[1],     xyz[2]+L, 0.5, tag=5)
    gmsh.model.geo.addPoint(xyz[0]+L,   xyz[1],     xyz[2]+L, 0.5, tag=6)
    gmsh.model.geo.addPoint(xyz[0]+L,   xyz[1]+L,   xyz[2]+L, 0.5, tag=7)
    gmsh.model.geo.addPoint(xyz[0],     xyz[1]+L,   xyz[2]+L, 0.5, tag=8)

    # Add Lines
    gmsh.model.geo.addLine(1, 2, tag=11) # Front Square
    gmsh.model.geo.addLine(2, 3, tag=12)
    gmsh.model.geo.addLine(3, 4, tag=13)
    gmsh.model.geo.addLine(4, 1, tag=14)
    gmsh.model.geo.addLine(5, 6, tag=15) # Back Square
    gmsh.model.geo.addLine(6, 7, tag=16)
    gmsh.model.geo.addLine(7, 8, tag=17)
    gmsh.model.geo.addLine(8, 5, tag=18)
    gmsh.model.geo.addLine(1, 5, tag=19) # Corner Connectors
    gmsh.model.geo.addLine(2, 6, tag=20)
    gmsh.model.geo.addLine(3, 7, tag=21)
    gmsh.model.geo.addLine(4, 8, tag=22)


    # Add Curve Loops
    gmsh.model.geo.addCurveLoop([11, 12, 13, 14], tag=201) # Front
    gmsh.model.geo.addCurveLoop([15, 16, 17, 18], tag=202) # Back NOT COUNTERCLOCKWISE
    gmsh.model.geo.addCurveLoop([11, 20,-15,-19], tag=203) # Bottom
    gmsh.model.geo.addCurveLoop([20, 16,-21,-12], tag=204) # Right
    gmsh.model.geo.addCurveLoop([-13,21, 17,-22], tag=205) # Top
    gmsh.model.geo.addCurveLoop([-19,-14,22, 18], tag=206) # Left

    # Add Plane Surface
    gmsh.model.geo.addPlaneSurface([201], tag=21)
    gmsh.model.geo.addPlaneSurface([202], tag=22)
    gmsh.model.geo.addPlaneSurface([203], tag=23)
    gmsh.model.geo.addPlaneSurface([204], tag=24)
    gmsh.model.geo.addPlaneSurface([205], tag=25)
    gmsh.model.geo.addPlaneSurface([206], tag=26)

    # Add volume
    gmsh.model.geo.addSurfaceLoop([21,22,23,24,25,26], tag=301)
    gmsh.model.geo.addVolume([301], tag=31)

    # Synchronize geometry with GMSH
    gmsh.model.geo.synchronize()


    # Add physical groups (NOTE: If physical groups are defined, GMSH only outputs mesh elements that belong to one such group)
    # See gmsh.option.setNumber("Mesh.SaveAll", 1)
    # Add Point (0D) Tag Groups
    gmsh.model.addPhysicalGroup(0, [1], tag=1, name='PT_FRNT_SW')
    gmsh.model.addPhysicalGroup(0, [2], tag=2, name='PT_FRNT_SE')
    gmsh.model.addPhysicalGroup(0, [3], tag=3, name='PT_FRNT_NE')
    gmsh.model.addPhysicalGroup(0, [4], tag=4, name='PT_FRNT_NW')
    gmsh.model.addPhysicalGroup(0, [5], tag=5, name='PT_BCK_SW')
    gmsh.model.addPhysicalGroup(0, [6], tag=6, name='PT_BCK_SE')
    gmsh.model.addPhysicalGroup(0, [7], tag=7, name='PT_BCK_NE')
    gmsh.model.addPhysicalGroup(0, [8], tag=8, name='PT_BCK_NW')
 

    # Add Edge (1D) Tag Groups
    gmsh.model.addPhysicalGroup(1, [11], tag=11, name='LN_FRNT_S')
    gmsh.model.addPhysicalGroup(1, [12], tag=12, name='LN_FRNT_E')
    gmsh.model.addPhysicalGroup(1, [13], tag=13, name='LN_FRNT_N')
    gmsh.model.addPhysicalGroup(1, [14], tag=14, name='LN_FRNT_W')
    gmsh.model.addPhysicalGroup(1, [15], tag=15, name='LN_BCK_S')
    gmsh.model.addPhysicalGroup(1, [16], tag=16, name='LN_BCK_E')
    gmsh.model.addPhysicalGroup(1, [17], tag=17, name='LN_BCK_N')
    gmsh.model.addPhysicalGroup(1, [18], tag=18, name='LN_BCK_W')
    gmsh.model.addPhysicalGroup(1, [19], tag=19, name='LN_CRNR_SW')
    gmsh.model.addPhysicalGroup(1, [20], tag=20, name='LN_CRNR_SE')
    gmsh.model.addPhysicalGroup(1, [21], tag=21, name='LN_CRNR_NE')
    gmsh.model.addPhysicalGroup(1, [22], tag=22, name='LN_CRNR_NW')

    # Add Surface (2D) Tag Groups
    gmsh.model.addPhysicalGroup(2, [21], tag=21, name='SRFC_FRNT')
    gmsh.model.addPhysicalGroup(2, [22], tag=22, name='SRFC_BCK')
    gmsh.model.addPhysicalGroup(2, [23], tag=23, name='SRFC_S')
    gmsh.model.addPhysicalGroup(2, [24], tag=24, name='SRFC_E')
    gmsh.model.addPhysicalGroup(2, [25], tag=25, name='SRFC_N')
    gmsh.model.addPhysicalGroup(2, [26], tag=26, name='SRFC_W')

    # Add Volume (3D) Tag Groups
    gmsh.model.addPhysicalGroup(3, [31], tag=31, name='VOL')

    # https://scicomp.stackexchange.com/questions/29385/creating-a-proper-quad-mesh-in-gmsh-for-an-i-shaped-geometry
    # gmsh.option.setNumber("Mesh.Algorithm", 9)  # 3 corresponds to quadrilateral elements


    # Generate 3D Mesh
    gmsh.model.mesh.generate(3)


    # Save out
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_file_base + '.msh')
    gmsh.finalize()


        # To visualize the model we can run the graphical user interface with
        # `gmsh.fltk.run()'. Here we run it only if "-nopopup" is not provided in the
        # command line arguments:
        # if '-nopopup' not in sys.argv:
        #     gmsh.fltk.run()



        # Note that starting with Gmsh 3.0, models can be built using other geometry
        # kernels than the default "built-in" kernel. To use the OpenCASCADE CAD kernel
        # instead of the built-in kernel, you should use the functions with the
        # `gmsh.model.occ' prefix.
        # Different CAD kernels have different features. With OpenCASCADE, instead of
        # defining the surface by successively defining 4 points, 4 curves and 1 curve
        # loop, one can define the rectangular surface directly with
        # gmsh.model.occ.addRectangle(.2, 0, 0, .1, .3)

        # After synchronization with the Gmsh model with
        # gmsh.model.occ.synchronize()
        # the underlying curves and points could be accessed with
        # gmsh.model.getBoundary().
        # See e.g. `t16.py', `t18.py', `t19.py' or `t20.py' for complete examples based
        # on OpenCASCADE, and `examples/api' for more.

        # This should be called when you are done using the Gmsh Python API:
        
