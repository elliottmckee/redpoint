'''

Just trying to make sure that I can interact with tagged regions efficiently


TAKEAWAYS:
    - Boundary integrators use coefficients to perform integatrion
        - They are keyed by using the attribute idx in a array
        - qdot = [0 0 10 0] will have attribute (3? 0-indexing sucks) have a load on it, others not


NOTES:
    - I think mfem doesn't actually maintain attribute data on the vertices...?
    - Might need to make dumbn wrappers to get DOFs we care about
    - Make sure you pay attention to AddDomainIntegrator vs. AddBoundaryIntegrator

'''

import os
from os.path import expanduser, join
import numpy as np

import mfem.ser as mfem


def run(order=1,  meshfile=''):

    # Hard coding dumb shit
    K_COND = 1.0
    QDOT = 10.0

    # tags for essential BC, and load applied BC tags
    load_attr_list  = [11]
    ess_attr_list   = [12, 13, 14]
    

    # Read in the mesh: 
    # Mesh(const char *mesh_file, int generate_edges, int refine, bool fix_orientation = True)
    mesh = mfem.Mesh(meshfile, 1, 1)
    mesh.UniformRefinement()


    # Define FE Space
    fec = mfem.H1_FECollection(order,  mesh.Dimension())
    fespace = mfem.FiniteElementSpace(mesh, fec)
    print('Number of finite element unknowns: ' + str(fespace.GetTrueVSize()))


    # Pull attribute sets from mesh
    bdr_attrs = mesh.bdr_attributes
    bdr_attrs_array = mesh.GetBdrAttributeArray()
    bdr_dofs = mesh.GetBdrArray(13)
    vol_attrs = mesh.GetAttributeArray()
    # bdr_vertex     = GetEdgeVertices(i)
    print(*bdr_attrs)
    print(bdr_attrs_array)
    print(bdr_dofs)
    print(vol_attrs)
    # print(bdr_vertex)
    # print('bdr_attrs: ', *mesh.bdr_attributes)


    # MAKE ARRAY MASKS FOR LOADS, ESSENTIAL DOFS

    # Get list of integers, up to the length of the maximum bdr_attribute being assigned
    bdr_attr_list = [*range(1, mesh.bdr_attributes.Max()+1)]
    print('bdr_attr_list: ', bdr_attr_list)

    # Make this into the dumb array that mfem expects
    ess_bdr_bool = mfem.intArray([1 if attr in ess_attr_list else 0 for attr in bdr_attr_list])
    print('ess_bdr_bool: ', *ess_bdr_bool)

    # Get TDOFS for these boundaries
    ess_tdof_list = mfem.intArray()
    fespace.GetEssentialTrueDofs(ess_bdr_bool, ess_tdof_list)
    # print(*ess_tdof_list)


    # Get load function
    # See ex2.py
    load_bdr_bool = [QDOT if attr in load_attr_list else 0 for attr in bdr_attr_list]
    print('load_bdr_bool: ', *load_bdr_bool)

    load_bdr_coeff = mfem.PWConstCoefficient(mfem.Vector(load_bdr_bool))



    # Define the solution x as a finite element grid function in fespace
    x = mfem.GridFunction(fespace)
    x.Assign(0.0)

    # 6. Set up the linear form b(.) corresponding to the right-hand side.
    # qdot_coeff = mfem.ConstantCoefficient(QDOT)
    b = mfem.LinearForm(fespace)
    b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(load_bdr_coeff))
    b.Assemble()

    # 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
    k_coeff = mfem.ConstantCoefficient(K_COND)
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.DiffusionIntegrator(k_coeff))
    a.Assemble()

    # 8. Form the linear system A X = B. This includes eliminating boundary
    #    conditions, applying AMR constraints, and other transformations.
    A = mfem.SparseMatrix()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B)
    # print("Size of linear system: " + str(A.Height()))

    # 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
    M = mfem.GSSmoother(A)
    mfem.PCG(A, M, B, X, 1, 200, 1e-12, 0.0)

    # 10. Recover the solution x as a grid function and save to file. The output
    #     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
    a.RecoverFEMSolution(X, b, x)
    

    sol_sock = mfem.socketstream("localhost", 19916)
    sol_sock.precision(8)
    sol_sock.send_solution(mesh, x)





if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description='Ex1 (Laplace Problem)')
    parser.add_argument('-m', '--mesh',
                        default='square.msh',
                        action='store', type=str,
                        help='Mesh file to use.')
    parser.add_argument('-o', '--order',
                        action='store', default=1, type=int,
                        help="Finite element order (polynomial degree) or -1 for isoparametric space.")

    args = parser.parse_args()
    parser.print_options(args)

    order = args.order
    meshfile = expanduser( join(os.path.dirname(__file__), 'gmsh', args.mesh) )

    run(order=order, meshfile=meshfile)


