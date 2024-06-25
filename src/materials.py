import numpy as np
from scipy.interpolate import interp1d

from . import ureg, Q_


'''
Materials Database and Helper Object

Why:
- I think that the eventual materials database is likely going to be some form of YAML, JSON, XML, etc. that I can easily read into a DICT
- Helper functions are going to make my life easier whenever I need to pull material properties at a specific temperature

TODO:
    - Make constant matprop implementation be able to handle vectors, since variable specification can
    - Make custom units registry with F, C not defined as dumb E&M units
    

'''


STD_DATABASE = { 
    
    # Aluminum 6061 Constant Matprops
    'ALU6061_CONST': {
        'rho':  Q_(2700.0, 'kg/m^3'), # Density
        'cp':   Q_(896.0, 'J/kg/kelvin'),  # Specific Heat
        'k':    Q_(167.0, 'W/m/kelvin'),  # Thermal Conductivity
        'emissivity': Q_(0.8, 'dimensionless')    # Black Body Emissivity Coefficient
    },

    # Aluminum 6061 Variable Matprops
    'ALU6061_VARIABLE': {
        # Source: https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061
        # Eyeballing rn
        'rho':  Q_(2700.0, 'kg/m^3'), # Density

        'cp': {'temperature': Q_(np.array([50.0, 100.0, 150.0, 200.0]), 'kelvin'), # Specific Heat 
                 'value': Q_(np.array([180.0, 500.0, 700.0, 810.0]), 'J/kg/kelvin')},

        'k': {'temperature': Q_(np.array([50.0, 100.0, 150.0, 200.0]), 'kelvin'), # Thermal Conductivity
                 'value': Q_(np.array([60.0, 98.0, 120.0, 138.0]), 'W/m/kelvin')},
        
        'emissivity': Q_(0.8, 'dimensionless')    #[] Black Body Emissivity Coefficient
    },

    # Aluminum 6061 CURVEFIT MATPROPS
    # TODO: Just leaving this here as a hint that the below implementation can handle callable functions, like curve fits
    'ALU6061_CURVEFIT': {
        # # Source: https://www.nist.gov/mml/acmd/aluminum-6061-t6-uns-aa96061
        # # Eyeballing rn
        # 'rho':  Q_(2700.0, 'kg/m^3'), #[kg/m^3] Density

        # 'cp': {'temperature': Q_(np.array([50.0, 100.0, 150.0, 200.0]), 'K'),
        #          'value': Q_(np.array([180.0, 500.0, 700.0, 810.0]), 'J/kg/K')},

        # 'k': {'temperature': Q_(np.array([50.0, 100.0, 150.0, 200.0]), 'K'),
        #          'value': Q_(np.array([60.0, 98.0, 120.0, 138.0]), 'W/m/K')},
        
        # 'emissivity': Q_(0.8, 'dimensionless' )    #[] Black Body Emissivity Coefficient
    },


}


class Material():
    '''
    Quick material class to add helper functions to interact with varied matprop specifications
    '''
    def __init__(self, properties):
        """
        Initialize the Material object with thermal properties.

        Parameters:
        - properties (dict): A dictionary containing 'rho' (density),
                            'cp' (specific heat), and 'k' (thermal conductivity).
                            'cp' and 'k' can be specified as functions of temperature
                            or static/constant values (or callables?)
        """
        #Assume rho, emissivity  constant
        self.rho = properties['rho']
        self.emissivity = properties['emissivity']

        # Check if cp, k are provided as a function or static/constant value
        self.cp_function = self.create_interpolation_function(properties['cp'])
        self.k_function = self.create_interpolation_function(properties['k'])

    def get_density(self):
        return self.rho    

    def get_emissivity(self):
        return self.emissivity
    
    def get_specific_heat(self, temperature=Q_(300,'kelvin')):
        return self.cp_function(temperature)

    def get_thermal_conductivity(self, temperature=Q_(300,'kelvin')):
        return self.k_function(temperature)

    def create_interpolation_function(self, property_data):
        """
        Create function to return a given material property, based on input data format provided

        Parameters:
        - property_data (Union[Pint Quantity, Dict, Callable]): Property values or function.

        Returns:
        - Callable: Interpolation function.
        """
        if callable(property_data):
            # If property_data is already a callable function, return it
            return property_data
        elif isinstance(property_data, Q_):
            # If property_data is a static/constant value, create a function that simply returns value
            return lambda T: property_data
        elif isinstance(property_data, dict):
            # If property_data is a dictionary, assume it contains temperatures, and values corresponding to those temperatures
            return lambda x_query: np.interp(x_query, property_data.get('temperature'), property_data.get('value'))
        else:
            raise ValueError("Invalid format for property_data. Must be a static Pint Quantity, dictionary with 'temperature' and 'value' keys, or function")

    
    


    




