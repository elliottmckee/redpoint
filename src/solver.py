''' Solver things '''

import numpy as np


class NodeMapping:
    def __init__(self, node_ids):
        self.node_id_to_index = {node_id: index for index, node_id in enumerate(node_ids)}
        self.index_to_node_id = {index: node_id for index, node_id in enumerate(node_ids)}

    def get_index(self, node_ids):
        if isinstance(node_ids, int):
            return self.node_id_to_index[node_ids]
        elif isinstance(node_ids, list):
            return [self.node_id_to_index[node_id] for node_id in node_ids]
        else:
            raise ValueError("Input must be an int or a list of ints.")

    def get_node_id(self, indices):
        if isinstance(indices, int):
            return self.index_to_node_id[indices]
        elif isinstance(indices, list):
            return [self.index_to_node_id[index] for index in indices]
        else:
            raise ValueError("Input must be an int or a list of ints.")



class TransientSolver():

    def __init__(self, Network) -> None:
        self.Network = Network
        self.NodeMap = NodeMapping(self.Network.nodes())
   

    def initialize(self) -> None:

        # Initialize State Vector
        self.X_init = np.array([]);
        for nodeID, data in self.Network.nodes(data=True):
            self.X_init = np.append(self.X_init, data['node'].T_init)
        self.X_init = np.reshape(self.X_init, (-1,1))
        # print(self.X_init)


        # Initialize Dynamics Matrix
        self.A = np.zeros((self.Network.number_of_nodes(), self.Network.number_of_nodes()))

        # For each edge
        for edge in self.Network.edges():

            # Get arbitrary DOFS into 0-indexed array
            n0 = self.NodeMap.get_index(edge[0])
            n1 = self.NodeMap.get_index(edge[1])

            # Aliasing
            R = self.Network.edges[edge]['R']
            mCp0 = self.Network.nodes[edge[0]]['node'].mCp
            mCp1 = self.Network.nodes[edge[1]]['node'].mCp

            #Eqn 1
            self.A[n0, n0] = self.A[n0, n0] - 1 / R / mCp0
            self.A[n0, n1] = self.A[n0, n1] + 1 / R / mCp0

            # Eqn 2
            self.A[n1, n0] = self.A[n1, n0] + 1 / R /mCp1
            self.A[n1, n1] = self.A[n1, n1] -1 / R /mCp1

        # print(self.A)
        # quit()

    def solve(self, t_bounds: tuple, dt: float):

        # Initialize
        t_vec = np.arange(t_bounds[0], t_bounds[1], dt)
        
        X = np.zeros((self.Network.number_of_nodes(), len(t_vec)))
        X[:,0] = np.squeeze(self.X_init)

        if 'T_fixed' in self.Network.BC_Set:
            fixed_node_DOFs, fixed_values = self.Network.BC_Set['T_fixed']

            # CONVERT DOFS to IDXS
            fixed_node_idxs = self.NodeMap.get_index(fixed_node_DOFs)

            # Set initial condition
            X[fixed_node_idxs,0] = fixed_values

        print(X[:,0]) 


        for i, t in enumerate(t_vec[:-1]):
            print('Time: ', t, '\nResults: ', X[:,i])

            # Get rates of change
            deltaX = np.dot(self.A, X[:,i]) * dt

            # Zero out fixed BC's
            deltaX[fixed_node_idxs] = 0

            # Update
            X[:,i+1] = X[:,i] + deltaX

        return X, t_vec


        





 














