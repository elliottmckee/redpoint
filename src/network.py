
import networkx as nx
import gmsh
import numpy as np
import matplotlib.pyplot as plt

from . import nodes



class ThermalNetwork(nx.Graph):
    '''
    Class for the specifcation of a Thermal Network Submodel

    Defined by a graph with thermal nodes as the nodes, and conductances as the edges

    TODO:
        -Will probably want to overwrite some of the add functionality and whatnot
    '''

    def __init__(self) -> None:
        # Just make a simple, empty networkX graph
        super(ThermalNetwork, self).__init__()

        self.BC_Set = {}

        pass



    def draw_network(self):

        #Unpack node coords
        node_coords = {node: [self.nodes[node]['node'].x, self.nodes[node]['node'].y] for node in self.nodes()}

        #Draw
        nx.draw(self, pos = node_coords, with_labels=True, font_weight='bold', node_color='skyblue', edge_color='gray', node_size=800)
        edge_labels = nx.get_edge_attributes(self, 'R')
        nx.draw_networkx_edge_labels(self, node_coords, edge_labels=edge_labels, font_color='red')
        plt.show()


    def plot_result(self, results, times):
        plt.figure()
        plt.plot(times, results.T)
        plt.show()



    def gmsh_1d_to_network(self, mesh_path: str, visualize: bool=False):
        '''
        Converts a 1D mesh from GMSH to thermal network

        Could make detection of boundaries smarter w/ smart use of gmsh.model.mesh.getNodes(dim1, includeBoundary=True)
        to include, and not include the boundary, and then take the reverse-intersection to get boundary nodes. Also likely an easy GMSH function somewhere 

        Assumes linear elements (order=1)
        Assumes uniform mesh spacing
        Assuming all thermal properties are 1.0
        '''

        #GMSH init
        gmsh.initialize()
        gmsh.open(mesh_path)

        #Alias to make GMSH more readable
        dim0 = 0
        dim1 = 1
        
        # Get Global nodelists, convert to dict
        nodeTags_global, nodeCoords_global, _ = gmsh.model.mesh.getNodes(dim1, includeBoundary=True)
        globCoords_dict = {node_id: coord for node_id, coord in zip(nodeTags_global, np.reshape(nodeCoords_global, (-1,3)))}

        # Get BC nodes
        tags_BC1= gmsh.model.getEntitiesForPhysicalGroup(dim0, 1)
        nodeTags_BC1, nodeCoords_BC1, _ = gmsh.model.mesh.getNodes(dim0, int(tags_BC1))
        # print('nodeTags_BC1 ', nodeTags_BC1)
        tags_BC2= gmsh.model.getEntitiesForPhysicalGroup(dim0, 2)
        nodeTags_BC2, nodeCoords_BC2, _ = gmsh.model.mesh.getNodes(dim0, int(tags_BC2))
        # print('nodeTags_BC2 ', nodeTags_BC2)

        # Get domain/interior nodes
        tags_dom = gmsh.model.getEntitiesForPhysicalGroup(dim1, 1)
        nodeTags_int, nodeCoords_int, _ = gmsh.model.mesh.getNodes(dim1, int(tags_dom))
        # print('nodeTags_int ', nodeTags_int)
   
        # Get the mesh elements for the entity (dim, tag):
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim1, 1)
        # Convert to list of tuples
        elem_tuples = [tuple(row) for row in np.reshape(elemNodeTags,(-1,2)).tolist()]
        # print(elem_tuples)



        # Add edges to define graph
        self.add_edges_from(elem_tuples)

        # Get mesh spacing
        h_el = np.round(nodeCoords_int[3] - nodeCoords_int[0], decimals=7)

        # Initialize Nodal properties
        for n_idx in nodeTags_global:
            # Make boundary nodes have half of the capacitance
            if n_idx in np.concatenate((nodeTags_BC1, nodeTags_BC2)):
                vol = h_el / 2;
            else:
                vol = h_el

            c_xyz = globCoords_dict[n_idx]
            self.nodes[n_idx]['node'] = nodes.ThermalNode(x=c_xyz[0], y=c_xyz[1], z=c_xyz[2], T_init=300, volume=vol, density=1.0, cp=1.0)

        # Initialize edge properties
        for edge in self.edges:
            self.edges[edge]['R'] = h_el / 1.0


        # Visualize the graph
        if visualize: self.draw_network();

























    def build_network_from_2D_mesh_centered_nodes(self, mesh_path: str, visualize: bool=False):
        '''
        Builds a Thermal Network object from a mesh

        Inputs:
            -pyGMSH mesh object?
        '''

        gmsh.initialize()
        gmsh.open(mesh_path)
        # gmsh.fltk.initialize()
        # gmsh.fltk.run()

        ### Get adjacency list ### 
        # Since GMSH doesn't do this shit, we have to do it ourselves
        # See: https://onelab.info/pipermail/gmsh/2019/012963.html

        # get all 2D elements
        elementType = gmsh.model.mesh.getElementType("quadrangle", 1)
        elemlist, nodelist = gmsh.model.mesh.getElementsByType(elementType)
        nodelist = np.reshape(nodelist, (-1,4))
        unique_nodelist = np.unique(nodelist)
        # print('elemlist:\n', elemlist)
        # print('nodelist:\n', nodelist)
        
        # Get nodal coordinates
        node_coords = np.zeros((len(unique_nodelist),3))
        for n_i in range(unique_nodelist.shape[0]):

            node_coords[n_i, :], _, _, _ = gmsh.model.mesh.getNode(unique_nodelist[n_i])
        # print(node_coords)


        # get edge nodelist
        bounds = gmsh.model.mesh.getElementEdgeNodes(3, 1)
        # print(f'edges: {bounds}')

        n_bounds_el = 4 # 4 for square
        n_bounds = len(elemlist)*n_bounds_el
        glob_boundlist = np.zeros((n_bounds,2))

        # Reshape to get a n_bounds x 2 list of edges/boundaries
        glob_boundlist = np.sort(np.reshape(bounds, (n_bounds, 2)), axis=1)
        # print('glob_boundlist', glob_boundlist)

        # Dumbass loop to get unique edges
        edgelist = set()
        node_map = np.floor( np.arange(n_bounds)/4 ).astype(int)

        for n_edge in range(0, n_bounds):

            # Remove self from boundary list
            temp_boundlist = np.copy(glob_boundlist)
            temp_boundlist[n_edge,:] = -1

            # Find where else current edge occurs
            idx_match = np.where( np.all(glob_boundlist[n_edge,:] == temp_boundlist, axis=1))
            
            # If that edge occurs elsewhere in the list
            if idx_match[0].size > 0:
                edge = tuple(sorted([elemlist[node_map[n_edge]], elemlist[node_map[idx_match[0][0]]]]))
                edgelist.add(edge)

        # print('edgelist', edgelist)


        # Add Nodes
        for el_idx in range(len(elemlist)):

            # Get element center coords
            n_idxs = nodelist[el_idx, :]  - 1
            curr_coords = node_coords[n_idxs, :]
            [x_c, y_c, z_c] = np.mean(node_coords[n_idxs, :], 0)

            # Get vector of xyz coords
            x_v, y_v, z_v =  curr_coords[:, 0],  curr_coords[:, 1],  curr_coords[:, 2]
            
            # Get Area
            area = abs(sum(x_v[i] * (y_v[i + 1] - y_v[i - 1]) for i in range(1, len(x_v) - 1)))
            # print(f"Area of the element: {area}")

            # Add Node
            self.add_node(elemlist[el_idx], node=nodes.ThermalNode(mCp=1.0*area, x=x_c, y=y_c, z=z_c, T_init=300.0))
            self.node_coords.update({elemlist[el_idx]: [x_c, y_c]})

        # Add Edges
        self.add_edges_from(edgelist, R = 10)

        # Visualize the graph
        if visualize: self.draw_network();

        gmsh.finalize()













