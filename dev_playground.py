'''
Wonder how far we'll get here
'''

import numpy as np
import scipy
import networkx as nx
import matplotlib.pyplot as plt
import gmsh
import math
import mfem.ser as mfem


from src import nodes, network, solver, ureg, Q_, fe
from src.materials import STD_DATABASE, Material



def assemble_system(components:list):
    
    # Get size of total system
    n_dofs_comp = [comp.get_n_dofs() for comp in components]
    n_dofs_tot = sum(n_dofs_comp)

    # Maintain map of which DOFS correpond to who (list of lists)
    # TODO: See Node Mapping in Solver.py
    # TODO: Check this!
    dofMap = [list(range(start, start+count)) for start, count in zip([0] + list(range(len(n_dofs_comp) + 1)), n_dofs_comp)]
    # print(dofMap)

    # Assemble global vectors
    X_glob_list = []
    for comp in components: 
        X_glob_list = X_glob_list+[*comp.X.GetDataArray()]

    B_glob_list = []
    for comp in components: 
        B_glob_list = B_glob_list+[*comp.B.GetDataArray()]

    X_glob = mfem.Vector(X_glob_list)
    B_glob = mfem.Vector(B_glob_list)

    # Assemble Global Matrices
    # Get list of matrices
    A_list = [[*comp.A.ToDenseMatrix().GetDataArray()] for comp in components]
    # Assemble block diagonally
    A_assembled_numpy = scipy.linalg.block_diag(*A_list)

    # Add connections
    # TODO


    # Convert back to sparse matrix
    csr_matrix_format = scipy.sparse.csr_matrix(A_assembled_numpy)
    indptr = np.array(csr_matrix_format.indptr, dtype=np.int32)
    indices = np.array(csr_matrix_format.indices, dtype=np.int32)
    data = np.array(csr_matrix_format.data, dtype=float)

    # Make Sparse Matrix
    A_glob_sparse = mfem.SparseMatrix([indptr, indices, data, n_dofs_tot, n_dofs_tot])


    return A_glob_sparse, B_glob, X_glob, dofMap

    





### MATPROP IMPLEMENTATION DEV
# ALU6061 = Material(STD_DATABASE['ALU6061_CONST'])
ALU6061 = Material(STD_DATABASE['ALU6061_VARIABLE'])

# print('Alu rho: ', ALU6061.get_density())
# print('Alu CP: ', ALU6061.get_specific_heat(Q_(100,'kelvin')))
# print('Alu k: ', ALU6061.get_thermal_conductivity(Q_(100,'kelvin')))
# print('Alu emis: ', ALU6061.get_emissivity())
# quit()

# Checking interpolation looks right
# x_query = Q_(np.linspace(0, 250, 50), 'kelvin')
# plt.figure()
# plt.plot(x_query, ALU6061.get_thermal_conductivity(x_query))
# plt.show()
# quit()
### 


MFEM_Square = fe.MFEM_Component('gmsh/square.msh', ALU6061)
# print('bdr_attrs', MFEM_Square.bdr_attrs)
MFEM_Square.mfem_set_bound_ess_tdofs([11, 12, 13])
MFEM_Square.mfem_add_const_bound_integrator_to_linform(10.0, [14])
MFEM_Square.mfem_assemble()
MFEM_Square.mfem_form_linear_system()

# M = mfem.GSSmoother(MFEM_Square.A)
# mfem.PCG(MFEM_Square.A, M, MFEM_Square.B, MFEM_Square.X, 1, 200, 1e-12, 0.0)
# MFEM_Square.bilinform.RecoverFEMSolution(MFEM_Square.X, MFEM_Square.linform, MFEM_Square.x)
# sol_sock = mfem.socketstream("localhost", 19916)
# sol_sock.precision(8)
# sol_sock.send_solution(MFEM_Square.mesh, MFEM_Square.x)


MFEM_Square_2 = fe.MFEM_Component('gmsh/square_2_0_0.msh', ALU6061)
# print('bdr_attrs', MFEM_Square.bdr_attrs)
MFEM_Square_2.mfem_set_bound_ess_tdofs([12])
# MFEM_Square_2.mfem_add_const_bound_integrator_to_linform(10.0, [14])
MFEM_Square_2.mfem_assemble()
MFEM_Square_2.mfem_form_linear_system()


# print(MFEM_Square.mesh.GetVertexArray())


A_G, B_G, X_G, DOFMAP  = assemble_system([MFEM_Square, MFEM_Square_2])


M = mfem.GSSmoother(A_G)
mfem.PCG(A_G, M, B_G, X_G, 1, 200, 1e-12, 0.0)



#Unpack solution
X_G_list = [*X_G.GetDataArray()]

X1 = [X_G_list[i] for i in DOFMAP[0]]
X2 = [X_G_list[i] for i in DOFMAP[1]]


MFEM_Square.bilinform.RecoverFEMSolution(mfem.Vector(X1), MFEM_Square.linform, MFEM_Square.x)
MFEM_Square_2.bilinform.RecoverFEMSolution(mfem.Vector(X2), MFEM_Square_2.linform, MFEM_Square_2.x)



sol_sock = mfem.socketstream("localhost", 19916)
sol_sock.precision(8)
sol_sock.send_solution(MFEM_Square.mesh, MFEM_Square.x)








# Thermal Network Solution
# ThermalNetwork = network.ThermalNetwork();

# ThermalNetwork.build_network_from_2D_mesh_centered_nodes('data/square.msh', visualize=False)

# ThermalSolver = solver.TransientSolver(ThermalNetwork)
# ThermalSolver.initialize()
# ThermalSolver.X_init[-1] = 400;

# X, t_vec = ThermalSolver.solve(t_bounds=(0,10), dt=.1)

# ThermalNetwork.plot_result(X, t_vec)









    


    

























