
import numpy as np
import inspect

import mfem.ser as mfem
from . import ureg, Q_
from .materials import STD_DATABASE



'''
TODO:
    - MAKE ROBUST TO CHANGES IN UNDERLYING DEFAULT QUANTITY REGISTRY DEFINITIONS 
    - Check that point-wise thermal conductivity works
    - Support SPARSE Matrices, not DENSE
'''



class MFEM_Component():
    '''
    Class representing an MFEM component

    Builds the default, static parts of the FE representation, but allows for external objects (heatloads, connections, etc.) 
    to go in and add themselves to things

    TODO:
        - 
    '''

    def __init__(self, meshfile, Material, T_init:Q_=Q_(300.0,'kelvin'), order:int=1) -> None:
        
        self.Material = Material
        self.T_init = T_init
        
        # Mesh(const char *mesh_file, int generate_edges, int refine, bool fix_orientation = True)
        self.mesh = mfem.Mesh(meshfile, 1, 1)
        # self.mesh.UniformRefinement(), self.mesh.UniformRefinement()
        self.dim = self.mesh.Dimension()
        

        # Define FE Space
        self.fe_coll = mfem.H1_FECollection(order, self.mesh.Dimension())
        self.fespace = mfem.FiniteElementSpace(self.mesh, self.fe_coll)

        self.n_raw_dofs = self.fespace.GetTrueVSize()

        # Attributes (both in actual attr list, and ready-to-be-masked format)
        self.bdr_attrs = [*self.mesh.bdr_attributes]
        self.bdr_attrs_range = [*range(1, max(self.bdr_attrs)+1)]
        self.attrs = [*self.mesh.attributes]
        self.attrs_range = [*range(1, max(self.attrs)+1)]


        self.mfem_initialize_solution_vector()
        self.mfem_initialize_linear_form()
        self.mfem_initialize_billinear_form()



    def mfem_initialize_solution_vector(self) -> None:
        self.x = mfem.GridFunction(self.fespace)
        self.x.Assign(self.T_init.to_base_units().magnitude)

        self.x_units = self.T_init.to_base_units().units

    
    def mfem_initialize_linear_form(self) -> None:
        self.linform = mfem.LinearForm(self.fespace)
        # b.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(load_bdr_coeff))
        # b.Assemble()
        pass
    
    def mfem_initialize_billinear_form(self) -> None:
        self.bilinform = mfem.BilinearForm(self.fespace)

        # TODO: Check that this sets the element-wise properties correctly
        # TODO: THIS IS BROKEN
        # k_vect = mfem.Vector(self.Material.get_thermal_conductivity( Q_([*self.x], self.x_units) ).magnitude)
        # k_coeff = mfem.PWConstCoefficient(k_vect)

        # TODO: This works tho
        k_coeff = mfem.ConstantCoefficient(self.Material.get_thermal_conductivity().magnitude)

        # Initialize default integrator
        self.bilinform.AddDomainIntegrator(mfem.DiffusionIntegrator(k_coeff))


    def mfem_set_bound_ess_tdofs(self, attr_tags) -> None:

        self.ess_bdr_list = mfem.intArray([1 if attr in attr_tags else 0 for attr in self.bdr_attrs_range])
        self.ess_tdof_arr = mfem.intArray()
        self.fespace.GetEssentialTrueDofs(self.ess_bdr_list, self.ess_tdof_arr)


    def mfem_add_const_bound_integrator_to_linform(self, const_value, attr_tags ) -> None:
        
        load_bdr_list = [const_value if attr in attr_tags else 0 for attr in self.bdr_attrs_range]
        load_bdr_coeff = mfem.PWConstCoefficient(mfem.Vector(load_bdr_list))

        self.linform.AddBoundaryIntegrator(mfem.BoundaryLFIntegrator(load_bdr_coeff))


    def mfem_assemble(self) -> None:
        self.linform.Assemble()
        self.bilinform.Assemble()


    def mfem_form_linear_system(self) -> None:
        self.A = mfem.SparseMatrix()
        self.B = mfem.Vector()
        self.X = mfem.Vector()
        self.bilinform.FormLinearSystem(self.ess_tdof_arr, self.x, self.linform, self.A, self.X, self.B)

    def get_n_dofs(self) -> int:
        return self.n_raw_dofs

    
        
        
        



        