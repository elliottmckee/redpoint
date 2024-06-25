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
from scipy.sparse import csr_matrix


from mfem.common.arg_parser import ArgParser
from os.path import expanduser, join, dirname
import numpy as np
from mfem import path

import mfem.ser as mfem
from mfem.ser import intArray


class ConductionOperator(mfem.PyTimeDependentOperator):
    def __init__(self, M, K, B, u, ess_tdof_list):
        mfem.PyTimeDependentOperator.__init__(self, fespace.GetTrueVSize()+1, 0.0)
        
        
        self.M = M
        # print(self.M.ToDenseMatrix().GetDataArray())
        self.K = K
        # print(self.K.GetDataArray())
        self.B = B
        # print(self.B.GetDataArray())

        self.ess_tdof_list = ess_tdof_list

        self.z = mfem.Vector(self.Height())

        self.M_prec = mfem.DSmoother()
        self.M_solver = mfem.CGSolver()
        self.M_solver.iterative_mode = False
        rel_tol = 1e-8
        self.M_solver.SetRelTol(rel_tol)
        self.M_solver.SetAbsTol(0.0)
        self.M_solver.SetMaxIter(30)
        self.M_solver.SetPrintLevel(0)
        self.M_solver.SetPreconditioner(self.M_prec)
        self.M_solver.SetOperator(self.M)

        # self.alpha = alpha
        # self.T = None
        # self.K = None
        # self.M = None
        # self.fespace = fespace

        # u_alpha_gf = mfem.GridFunction(self.fespace)
        # u_alpha_gf.SetFromTrueDofs(u)

        # self.T_solver = mfem.CGSolver()
        # self.T_prec = mfem.DSmoother()


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
        # print(self.K.GetDataArray())
        # print(u.GetDataArray())
        self.K.Mult(u, self.z)
        # print(*self.z)
        # quit()

        
        self.z.Neg()   # z = -z
        self.z += self.B

        
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


# class InitialTemperature(mfem.PyCoefficient):
#     def EvalValue(self, x):
#         xx = np.array(x)
#         norm2 = np.sqrt(np.sum(xx**2))
#         if norm2 < 0.5:
#             return 2.0
#         return 1.0

class InitialTemperature(mfem.PyCoefficient):
    def EvalValue(self, x):
        if x == 0.0:
            return 300.0
        return 300.0


# Input Options
order = 1
dt = 0.001
t_final = 1
alpha = 0.1
QDOT = 0.0
# kappa = args.kappa
visualization = False
vis_steps = 10
ode_solver = mfem.ForwardEulerSolver()
# ode_solver = mfem.BackwardEulerSolver()
# ode_solver = mfem.RK4Solver()
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
u_gf.ProjectCoefficient(u_0) # Project IC function onto grid function
u = mfem.Vector()
u_gf.GetTrueDofs(u) #Save just true (non-hanging) DOFs into u


# Dumb hack to check augmenting works
u_aug = mfem.Vector(np.append(u.GetDataArray(), 305))


# Initialize initial boundaries
dirichlet_mkr = [0,0]
neumann_mkr = [int(not elem) for elem in dirichlet_mkr]

ess_bdr = mfem.intArray(dirichlet_mkr)
neu_bdr = mfem.intArray(neumann_mkr)
# ess_bdr = mfem.intArray()
# if mesh.bdr_attributes.Size():
#     ess_bdr.SetSize(mesh.bdr_attributes.Max())
#     ess_bdr.Assign(0) # Initialize entire array to zeros


# Get boundary dofs
ess_tdof_list = mfem.intArray()
fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list)
# print(*ess_tdof_list)


# Assemble Mass Matrix (BILLINEAR)
Mmat = mfem.SparseMatrix()
M = mfem.BilinearForm(fespace)
M.AddDomainIntegrator(mfem.MassIntegrator())
M.Assemble()
M.FormSystemMatrix(ess_tdof_list, Mmat)
# print(Mmat.ToDenseMatrix().GetDataArray())


### Dumb hacking at the matrices to get a augmented matrix
#Expand Matrix
Mmat_aug_np = np.zeros((fe_size+1, fe_size+1))
Mmat_aug_np[0:fe_size, 0:fe_size] = Mmat.ToDenseMatrix().GetDataArray()

#Add Values corresponding to final node
Mmat_aug_np[fe_size,fe_size] = 0.08333333
# Mmat_aug_np[1, fe_size] = 0.04166667
# Mmat_aug_np[fe_size, 1] = 0.04166667

#Convert to Mfem densematrix
Mmat_aug = mfem.DenseMatrix(fe_size+1, fe_size+1)
Mmat_aug.Assign(Mmat_aug_np)
# print(Mmat_aug.GetDataArray())
### 



# Convert the MFEM Matrix to a SparseMatrix
csr_matrix_format = csr_matrix(Mmat_aug_np)

# Extract the indptr, indices, and data arrays
indptr = np.array(csr_matrix_format.indptr, dtype=np.int32)
indices = np.array(csr_matrix_format.indices, dtype=np.int32)
data = np.array(csr_matrix_format.data, dtype=float)

# Make Sparse Matrix
Mmat_aug_sparse = mfem.SparseMatrix([indptr, indices, data, fe_size+1, fe_size+1])
# print(Mmat_aug_sparse.ToDenseMatrix().GetDataArray())




# Assemble Stiffness Matrix (BILLINEAR)
Kmat = mfem.SparseMatrix()
alpha_coeff = mfem.ConstantCoefficient(alpha)
K = mfem.BilinearForm(fespace)
K.AddDomainIntegrator(mfem.DiffusionIntegrator(alpha_coeff))
K.Assemble()
K.FormSystemMatrix(ess_tdof_list, Kmat)

Kmat_dense = Kmat.ToDenseMatrix()
# print(Kmat.ToDenseMatrix().GetDataArray())

### Dumb hacking at the matrices to get a augmented matrix
#Expand Matrix
Kmat_aug_np = np.zeros((fe_size+1, fe_size+1))
Kmat_aug_np[0:fe_size, 0:fe_size] = Kmat.ToDenseMatrix().GetDataArray()


#Add Values corresponding to final node
Kmat_aug_np[fe_size,fe_size] = 0.4
Kmat_aug_np[1,1] += 0.4
Kmat_aug_np[1, fe_size] = -0.4
Kmat_aug_np[fe_size, 1] = -0.4
# print(Kmat_aug_np)
#Convert to Mfem densematrix
Kmat_aug = mfem.DenseMatrix(fe_size+1, fe_size+1)
Kmat_aug.Assign(Kmat_aug_np)
### 






# Assemble RHS Loadings (LINEAR)
B = mfem.LinearForm(fespace)
QDOT_coeff = mfem.ConstantCoefficient(QDOT)
B.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(QDOT_coeff), neu_bdr)
# b.AddDomainIntegrator()
B.Assemble()
# print('B: ', B.GetDataArray())

### Dumb hacking at the matrices to get a augmented matrix

B_aug = mfem.Vector([0,0,0,0,0,0])






# Do we need to so some sort of essential BC removal here?

# FormLinearSystem does do some fancy assembly things- see documentation. 
# Not sure if can leverage any of them. But maybe look into in future

# K_TEST = mfem.SparseMatrix()
# B_TEST = mfem.Vector()
# X_TEST = mfem.Vector()
# K.FormLinearSystem(mfem.intArray([1, 0]), u_gf, B, K_TEST, X_TEST, B_TEST)
# quit()



# Initialize the conduction operator and the visualization.
# print('M_aug: ', Mmat_aug_sparse.ToDenseMatrix().GetDataArray())
# print('K_aug: ', Kmat_aug.GetDataArray())
# print('B_aug:', u_aug.GetDataArray())
# print('ess_tdof_list: ', *ess_tdof_list)

oper = ConductionOperator(Mmat_aug_sparse, Kmat_aug, B_aug, u_aug, ess_tdof_list)


u_gf.SetFromTrueDofs(mfem.Vector([*u_aug[0:5]]))


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
ti = 0
last_step = False

data_out = np.reshape(np.array(u_gf.GetDataArray()), (1, -1))

while not last_step:
    ti = ti + 1
    if (t + dt >= t_final - dt/2):
        last_step = True

    t, dt = ode_solver.Step(u_aug, t, dt)



    print(*u_aug)

    u_aug[5] = 305


    # u[0] = 300.0

    if (last_step or (ti % vis_steps) == 0):
        # if True:
        print("step " + str(ti) + ", t = " + "{:g}".format(t))

        u_gf.SetFromTrueDofs(mfem.Vector([*u_aug[0:5]])) # Pull back onto original grid (with hanging vertsif exist?)

        data_out = np.concatenate((data_out,  np.reshape(np.array(u_gf.GetDataArray()), (1, -1))), axis=0)
        # if (visualization):
        #     sout << "solution\n" << mesh << u_gf
        #     sout.flush()

    
    # oper.SetParameters(u)

# RECOVER FEM SOLUTION?
plt.figure()
plt.plot(mesh_coords, np.array(data_out).T,  marker='o', linestyle='')
plt.show()


# 9. Save the final solution. This output can be viewed later using GLVis:
#    "glvis -m ex16.mesh -g ex16-final.gf".
# u_gf.Save('ex16-final.gf', 8)