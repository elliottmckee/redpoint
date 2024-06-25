'''
For definitions of nodes and otherwise
'''

import numpy as np
import scipy
import networkx as nx
from typing import Literal


_MATERIALS = Literal['DEFAULT, ''Al_7050', 'SS316']


class ThermalNode:
    '''
    The fundamental data structure in this architecture.

    Holds definition information
    - mCp, capacitance
    - x,y,z location

    And state information
    - T, temperature

    TODO:
        - Properties that vary with temperature
    '''


    def __init__(self, x: float=0.0, y: float=0.0, z: float=0.0, T_init: float=300, **kwargs):
        
        self.x = x
        self.y = y
        self.z = z
        self.T_init = T_init
        
        if 'mCp' in kwargs:
            self.mCp = kwargs['mCp']
        elif all(key in kwargs for key in ['volume', 'density', 'cp']):
            # Calculate mCp from volume, density, and cp if all are provided
            self.mCp = kwargs['volume'] * kwargs['density'] * kwargs['cp']
        else:
            # Raise an error for invalid arguments
            raise ValueError("Invalid arguments. Provide either ['mCp'] or ['volume', 'density', and 'cp']")
        




class MaterialNode(ThermalNode):
    '''
    One step higher than a ThermalNode. Has to be assigned a material, by which it derives its properties
    '''
    # Keeping this material type hint checking as a reminder to myself of this technique for later
    def __init__(self, material: _MATERIALS='DEFAULT', volume: float=0.0, x: float=0.0, y: float=0.0, z: float=0.0, T_init: float=0.0):
        self.material = material
        self.volume = volume
        self.x = x
        self.y = y
        self.z = z
        self.T_init = T_init

        self.set_capacitance(self)

    def set_capacitance(self):
        pass
        















