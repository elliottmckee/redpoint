'''
   MFEM example 16

   How to run:
      python <arguments>

   Example of arguments:
      ex16.py -m inline-tri.mesh
      ex16.py -m disc-nurbs.mesh -tf 2
      ex16.py -s 1 -a 0.0 -k 1.0
      ex16.py -s 2 -a 1.0 -k 0.0
      ex16.py -s 3 -a 0.5 -k 0.5 -o 4
      ex16.py -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
      ex16.py -m fichera-q2.mesh
      ex16.py -m escher.mesh
      ex16.py -m beam-tet.mesh -tf 10 -dt 0.1
      ex16.py -m amr-quad.mesh -o 4 -r 0
      ex16.py -m amr-hex.mesh -o 2 -r 0

'''
import sys
import matplotlib.pyplot as plt


from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np
from mfem import path

import mfem.ser as mfem
from mfem.ser import intArray


class ConductionOperator(mfem.PyTimeDependentOperator):
    def __init__(self, fespace, alpha, u):
        mfem.PyTimeDependentOperator.__init__(self, fespace.GetTrueVSize(), 0.0)
        rel_tol = 1e-8
        self.alpha = alpha
        self.T = None
        self.K = None
        self.M = None
        self.fespace = fespace

        u_alpha_gf = mfem.GridFunction(self.fespace)
        u_alpha_gf.SetFromTrueDofs(u)


        self.ess_tdof_list = intArray()
        self.Mmat = mfem.SparseMatrix()
        self.Kmat = mfem.SparseMatrix()
        
        # self.T_solver = mfem.CGSolver()
        # self.T_prec = mfem.DSmoother()
        self.z = mfem.Vector(self.Height())


        self.M = mfem.BilinearForm(fespace)
        self.M.AddDomainIntegrator(mfem.MassIntegrator())
        self.M.Assemble()
        self.M.FormSystemMatrix(self.ess_tdof_list, self.Mmat)

        # Pulling this in from setparms since want this to be constant
        self.K = mfem.BilinearForm(fespace)
        self.K.AddDomainIntegrator(mfem.DiffusionIntegrator(mfem.ConstantCoefficient(alpha)))
        self.K.Assemble()
        self.K.FormSystemMatrix(self.ess_tdof_list, self.Kmat)

        # Linear form (boundary qdots)
        self.b = mfem.LinearForm(fespace)
        self.b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(mfem.ConstantCoefficient(1.0)))
        self.b.Assemble()

        self.M_prec = mfem.DSmoother()
        self.M_solver = mfem.CGSolver()
        self.M_solver.iterative_mode = False
        self.M_solver.SetRelTol(rel_tol)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(30)
        self.M_solver.SetPrintLevel(0)
        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.SetOperator(self.Mmat)

        # self.T_solver.iterative_mode = False
        # self.T_solver.SetRelTol(rel_tol)
        # self.T_solver.SetAbsTol(0.0)
        # self.T_solver.SetMaxIter(100)
        # self.T_solver.SetPrintLevel(0)
        # self.T_solver.SetPreconditioner(self.T_prec)

        # self.SetParameters(u)

    def Mult(self, u, du_dt):
        # Compute:
        #  du_dt = M^{-1}*-K(u) for du_dt
        self.Kmat.Mult(u, self.z)
        self.z.Neg()   # z = -z
        self.z += self.b
        self.M_solver.Mult(self.z, du_dt)

    def ImplicitSolve(self, dt, u, du_dt):
        # Solve the equation:
        #    du_dt = M^{-1}*[-K(u + dt*du_dt)]
        #    for du_dt
        if self.T is None:
            self.T = mfem.Add(1.0, self.Mmat, dt, self.Kmat)
            current_dt = dt
            self.T_solver.SetOperator(self.T)
        self.Kmat.Mult(u, self.z)
        self.z.Neg()
        self.T_solver.Mult(self.z, du_dt)

    # def SetParameters(self, u):
    #     u_alpha_gf = mfem.GridFunction(self.fespace)
    #     u_alpha_gf.SetFromTrueDofs(u)
    #     for i in range(u_alpha_gf.Size()):
    #         u_alpha_gf[i] = self.kappa + self.alpha * u_alpha_gf[i]

    #     self.K = mfem.BilinearForm(self.fespace)
    #     u_coeff = mfem.GridFunctionCoefficient(u_alpha_gf)
    #     self.K.AddDomainIntegrator(mfem.DiffusionIntegrator(u_coeff))
    #     self.K.Assemble()
    #     self.K.FormSystemMatrix(self.ess_tdof_list, self.Kmat)
    #     self.T = None


class InitialTemperature(mfem.PyCoefficient):
    def EvalValue(self, x):
        xx = np.array(x)
        norm2 = np.sqrt(np.sum(xx**2))
        if norm2 < 0.5:
            return 2.0
        return 1.0

# class InitialTemperature(mfem.PyCoefficient):
#     def EvalValue(self, x):
#         if x == 0.0:
#             return 400.0
#         return 300.0


# Input Options
order = 1
dt = 0.01
t_final = 1
alpha = 0.1
# kappa = args.kappa
visualization = False
vis_steps = 10
ode_solver = mfem.ForwardEulerSolver()
# ode_solver = mfem.BackwardEulerSolver()
ode_solver = mfem.RK4Solver()
meshfile = 'data/line_L1p0_r2.msh'


# Read the mesh, get coordinates of vertices for plotting
mesh = mfem.Mesh(meshfile, 1, 1)
dim = mesh.Dimension()
mesh_coords = np.array(mesh.GetVertexArray())


# Define the vector finite element space representing the current and the initial temperature, u_ref.
fe_coll = mfem.H1_FECollection(order, dim)
fespace = mfem.FiniteElementSpace(mesh, fe_coll)
fe_size = fespace.GetTrueVSize()
print("Number of temperature unknowns (pre-BC elim): " + str(fe_size))


# Set the initial conditions for u. 
u_gf = mfem.GridFunction(fespace)
u_0 = InitialTemperature()
u_gf.ProjectCoefficient(u_0) # Project IC function onto grid
u = mfem.Vector()
u_gf.GetTrueDofs(u) #Save just true (non-hanging) DOFs into u
u_gf.SetFromTrueDofs(u)


# Initialize the conduction operator and the visualization.
oper = ConductionOperator(fespace, alpha , u)


# GLVIS THINGS
# mesh.Print('ex16.mesh', 8)
# u_gf.Save('ex16-init.gf', 8)
# if visualization:
#     sout = mfem.socketstream("localhost", 19916)
#     sout.precision(8)
#     sout << "solution\n" << mesh << u_gf
#     sout << "pause\n"
#     sout.flush()
#     print("GLVis visualization paused.")
#     print(" Press space (in the GLVis window) to resume it.")


# Perform time-integration 
ode_solver.Init(oper)
t = 0.0
ti = 1
last_step = False

data_out = np.reshape(np.array(u_gf.GetDataArray()), (1, -1))

while not last_step:
    if (t + dt >= t_final - dt/2):
        last_step = True

    t, dt = ode_solver.Step(u, t, dt)

    if (last_step or (ti % vis_steps) == 0):
        # if True:
        print("step " + str(ti) + ", t = " + "{:g}".format(t))
        u_gf.SetFromTrueDofs(u) # Pull back onto original grid (with hanging vertsif exist?)

        data_out = np.concatenate((data_out,  np.reshape(np.array(u_gf.GetDataArray()), (1, -1))), axis=0)
        # if (visualization):
        #     sout << "solution\n" << mesh << u_gf
        #     sout.flush()

    ti = ti + 1
    # oper.SetParameters(u)

# RECOVER FEM SOLUTION?


plt.figure()
plt.plot(mesh_coords, np.array(data_out).T,  marker='o', linestyle='')
plt.show()


# 9. Save the final solution. This output can be viewed later using GLVis:
#    "glvis -m ex16.mesh -g ex16-final.gf".
# u_gf.Save('ex16-final.gf', 8)