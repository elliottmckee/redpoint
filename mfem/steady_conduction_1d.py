'''
First cut at trying to get steady 1-D conduction to work
'''

import os
# from os.path import expanduser, join
import numpy as np

import mfem.ser as mfem


class InitialTemperature(mfem.PyCoefficient):
    def EvalValue(self, x):
        if x <= 1.0e-6:
            return 400.0
        return 300.0



if __name__ == '__main__':
    
    # Dumb Inputs
    mesh_file = 'data/line_L1p0_r2.msh'
    K_COND = 1.0
    QDOT = 10.0


    # Read the mesh from the given mesh file
    mesh = mfem.Mesh(mesh_file, 1, 1)


    # 3. Define a finite element space on the mesh. Here we use H1 continuous Lagrange finite elements of the given order.
    fec = mfem.H1_FECollection(1,  mesh.Dimension())
    fespace = mfem.FiniteElementSpace(mesh, fec)
    print('Number of finite element unknowns: ' +
          str(fespace.GetTrueVSize()))
    
    # 4. Extract the list of all the boundary DOFs. These will be marked as Dirichlet in order to enforce zero boundary conditions.
    boundary_dofs = mfem.intArray([0]) # Just setting RHS to be dirichlet rn
    # fespace.GetBoundaryTrueDofs(boundary_dofs)
    # print(boundary_dofs.ToList())

    # 5. Define the solution x as a finite element grid function in fespace. Set the initial guess to zero, which also sets the boundary conditions.
    x = mfem.GridFunction(fespace)
    x0 = InitialTemperature()
    # x.Assign(300.0)
    x.ProjectCoefficient(x0)


    # 6. Set up the linear form b(.) corresponding to the right-hand side. (See ex2.py)
    q = mfem.ConstantCoefficient(QDOT)
    b = mfem.LinearForm(fespace)
    b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(q))
    b.Assemble()

    # 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
    k = mfem.ConstantCoefficient(K_COND)
    a = mfem.BilinearForm(fespace)
    a.AddDomainIntegrator(mfem.DiffusionIntegrator(k))
    a.Assemble()

    # 8. Form the linear system A X = B. This includes eliminating boundary conditions, applying AMR constraints, and other transformations.
    A = mfem.SparseMatrix()
    B = mfem.Vector()
    X = mfem.Vector()
    a.FormLinearSystem(boundary_dofs, x, b, A, X, B)
    print("Size of linear system: " + str(A.Height()))

    # 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
    M = mfem.GSSmoother(A)
    mfem.PCG(A, M, B, X, 1, 300, 1e-12, 0.0)

    # 10. Recover the solution x as a grid function and save to file. The output
    #     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
    a.RecoverFEMSolution(X, b, x)
    # x.Save('sol.gf')
    # mesh.Save('mesh.mesh')

    for w in x:
        print(w)


    # IDK why this doesn't work but it works in example 2  

    nodes = mesh.GetNodes()
    print(nodes)






    


    