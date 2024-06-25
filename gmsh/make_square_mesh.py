import gmsh
import sys



if __name__ == '__main__':


    xyz = (0.0, 0.0, 0.0)
    L = 1.0

    output_file_base = 'square'


    gmsh.initialize()
    gmsh.model.add("model1")

    # Add Points
    # gmsh.model.geo.addPoint(0, 0, 0, meshsize=1, tag=1)
    gmsh.model.geo.addPoint(xyz[0],     xyz[1],     xyz[2], 0.5, tag=1)
    gmsh.model.geo.addPoint(xyz[0]+L,   xyz[1],     xyz[2], 0.5, tag=2)
    gmsh.model.geo.addPoint(xyz[0]+L,   xyz[1]+L,   xyz[2], 0.5, tag=3)
    gmsh.model.geo.addPoint(xyz[0],     xyz[1]+L,   xyz[2], 0.5, tag=4)

    # Add Lines
    gmsh.model.geo.addLine(1, 2, tag=11)
    gmsh.model.geo.addLine(2, 3, tag=12)
    gmsh.model.geo.addLine(3, 4, tag=13)
    gmsh.model.geo.addLine(4, 1, tag=14)

    # Add Curve Loops
    gmsh.model.geo.addCurveLoop([11, 12, 13, 14], tag=99)

    # Add Plane Surface
    gmsh.model.geo.addPlaneSurface([99], tag=21)

    # Synchronize geometry with GMSH
    gmsh.model.geo.synchronize()


    # Add physical groups (NOTE: If physical groups are defined, GMSH only outputs mesh elements that belong to one such group)
    # See gmsh.option.setNumber("Mesh.SaveAll", 1)
    # Add Point (0D) Tag Groups
    gmsh.model.addPhysicalGroup(0, [1], tag=1, name='PT_SW')
    gmsh.model.addPhysicalGroup(0, [2], tag=2, name='PT_SE')
    gmsh.model.addPhysicalGroup(0, [3], tag=3, name='PT_NE')
    gmsh.model.addPhysicalGroup(0, [4], tag=4, name='PT_NW')
 
    # Add Edge (1D) Tag Groups
    gmsh.model.addPhysicalGroup(1, [11], tag=11, name='LN_S')
    gmsh.model.addPhysicalGroup(1, [12], tag=12, name='LN_E')
    gmsh.model.addPhysicalGroup(1, [13], tag=13, name='LN_N')
    gmsh.model.addPhysicalGroup(1, [14], tag=14, name='LN_W')

    # Add Surface (2D) Tag Groups
    gmsh.model.addPhysicalGroup(2, [21], tag=21, name='SRFC')

    # https://scicomp.stackexchange.com/questions/29385/creating-a-proper-quad-mesh-in-gmsh-for-an-i-shaped-geometry
    # gmsh.option.setNumber("Mesh.Algorithm", 9)  # 3 corresponds to quadrilateral elements


    # Generate 2D Mesh
    gmsh.model.mesh.generate(2)



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
        
