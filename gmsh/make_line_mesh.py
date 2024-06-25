import gmsh
import sys



if __name__ == '__main__':

    # Values at the start of refinement
    L = 1.0
    h_0 = 1.0

    output_file_base = 'line_L1p0_r'

    for r_val in  [0,1,2,5]:

        h = h_0 * 0.5**r_val

        gmsh.initialize()
        gmsh.model.add("model1")

        # Add Points
        # gmsh.model.geo.addPoint(0, 0, 0, meshsize=1, tag=1)
        gmsh.model.geo.addPoint(0, 0, 0, meshSize=h, tag=1)
        gmsh.model.geo.addPoint(1, 0, 0, meshSize=h, tag=2)

        # Add Line
        gmsh.model.geo.addLine(1, 2, tag=11) 

        # Synchronize geometry with GMSH
        gmsh.model.geo.synchronize()

        # Add physical groups (NOTE: If physical groups are defined, GMSH only outputs mesh elements that belong to one such group)
        # See gmsh.option.setNumber("Mesh.SaveAll", 1)
        
        # Add Point (0D) Tag Groups
        gmsh.model.addPhysicalGroup(0, [1], tag=1, name='PT_L')
        gmsh.model.addPhysicalGroup(0, [2], tag=2, name='PT_R')

        # Add Edge (1D) Tag Groups
        gmsh.model.addPhysicalGroup(1, [11], tag=11, name='LN')

        # Generate Mesh
        gmsh.model.mesh.generate(1)

        # ... and save it to disk
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(output_file_base + str(r_val)+'.msh')




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
        gmsh.finalize()
