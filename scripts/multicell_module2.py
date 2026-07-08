from typing import Union
import pandas as pd
import numpy as np
import warnings
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import hummuspy
import re
azim = np.random.randint(-180, 180)

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)


def split_layer_name(s, separator='_'):
    # Use re.escape to handle special characters in the separator
    pattern = rf"(.*){re.escape(separator)}(\w+)$"
    match = re.match(pattern, s)
    if match:
        return match.group(1), match.group(2)
    else:
        return s, ''  # In case there's no separator, return the original string and empty


def illustrate_multicell(
    lamb,
    figsize=(12, 10),
    azim = 90,
    elev = 20,
    display_layer_axis=True,
    display_self_proba=True,
    display_layer_names=False,
    alpha_layers = 0.5,
    cell_communication_layer_name="cell_communication"):

    # Adjust the viewing angle to be more horizontal
    fig = plt.figure(figsize=(12, 10), dpi=200)
    ax = fig.add_subplot(111, projection='3d')

    # Set a lower elevation angle to create a more horizontal view
    ax.view_init(elev=elev, azim=azim)  # Set elevation to 20 for a more side view

    layers = lamb.index.values
    lamb = lamb.loc[layers, layers]
    cell_groups = set([split_layer_name(layer)[0] for layer in layers if layer != cell_communication_layer_name])
    cell_groups.add(cell_communication_layer_name)
    
    # Define the x-coordinate for each cell type's column
    celltypes = set([split_layer_name(layer)[0] for layer in layers if layer != cell_communication_layer_name])
    x_positions = {celltype: i for i, celltype in enumerate(celltypes)}
    x_max, x_min = max(x_positions.values()), min(x_positions.values())
    x_positions[cell_communication_layer_name] = (x_min + x_max) / 2

    layer_types = list(set([split_layer_name(layer)[-1] for layer in layers if layer != cell_communication_layer_name]))
    layer_heights = {
        layer_type: i for i, layer_type in enumerate(layer_types)
    }
    layer_heights[cell_communication_layer_name] = len(layer_heights)
    layer_types.append(cell_communication_layer_name)

    z_layers = {
        layer: layer_heights[split_layer_name(layer)[-1]]
        for layer in layers if layer != cell_communication_layer_name}
    z_layers[cell_communication_layer_name] = layer_heights[cell_communication_layer_name]

    if len(cell_groups) < 10:
        celltype_colors = {
            celltype: color for celltype, color in zip(cell_groups, plt.cm.tab10.colors)}
    elif len(cell_groups) < 20:
        celltype_colors = {
            celltype: color for celltype, color in zip(cell_groups, plt.cm.tab20.colors)}
    else:
        celltype_colors = {
            celltype: plt.cm.tab20.colors[i] for i, celltype in enumerate(cell_groups)}

    # Draw the "cell_communication" layer as a large plane spanning across columns
    Z = z_layers[cell_communication_layer_name]

    u = np.linspace(x_min - 0.4, x_max + 0.4, 2)
    v = np.linspace(-0.3, 0.3, 2)
    U, V = np.meshgrid(u, v)
    W = Z * np.ones_like(U)
    ax.plot_surface(
        U, V, W,
        color=celltype_colors[cell_communication_layer_name],
        label=cell_communication_layer_name,
        alpha=alpha_layers)

    if display_layer_names:
        ax.text(
            (x_max+x_min)/2+0.3, 0.3, Z,  # Position of the layer name
            cell_communication_layer_name,  # Layer name
            "x", color="black", ha="center", va="top", fontsize=10)

    # Draw planes for each layer in their respective columns
    for layer in layers:
        if layer == cell_communication_layer_name:
            continue  # Skip since cell_communication is already drawn

        celltype = split_layer_name(layer)[0]
        x = x_positions.get(celltype, 0)
        z = z_layers[layer]
        color = celltype_colors.get(celltype, "gray")

        # Draw each layer's plane at specified x and z positions
        u = np.linspace(x - 0.4, x + 0.4, 2)
        v = np.linspace(-0.3, 0.3, 2)
        U, V = np.meshgrid(u, v)
        W = z * np.ones_like(U)
        ax.plot_surface(U, V, W, color=color, alpha=alpha_layers)
        if display_layer_names:
            ax.text(
                x, 0.3, z,  # Position of the layer name
                layer,  # Layer name
                "x", color="black", ha="center", va="top", fontsize=10)

    # Draw edges based on probabilities with positions and heights specified
    for i, layer_from in enumerate(lamb.index):
        for j, layer_to in enumerate(lamb.columns):
            prob = lamb.iloc[i, j]
            if i == j:
                if layer_from != cell_communication_layer_name:
                    x1, z1 = x_positions.get(
                        split_layer_name(layer_from)[0], 0), z_layers[layer_from]
                else:
                    x1, z1 = x_positions.get(layer_from, 0), z_layers[layer_from]
                if display_self_proba:
                    ax.text(
                        x1, 0, z1,  # Position of the probability of self-transition
                        round(prob, 2),  # Displayed probability of self-transition
                        "x", color="black", ha="center", va="center", fontsize=10)
            elif prob != 0:
                x1, z1 = x_positions.get(
                    split_layer_name(layer_from)[0], 0), z_layers[layer_from]
                x2, z2 = x_positions.get(
                    split_layer_name(layer_to)[0], 0), z_layers[layer_to]
                if layer_from == cell_communication_layer_name:
                    x1 = x2
                elif layer_to == cell_communication_layer_name:
                    x2 = x1

                if x1 > x2:
                    y1 = 0.2
                    y2 = 0.2
                elif x1 < x2:
                    y1 = -0.2
                    y2 = -0.2

                if z1 > z2:
                    y1, y2 = -0.0, -0.0,
                    x1-=0.2
                    x2-=0.2
                elif z1 < z2:
                    y1, y2 = 0.0, 0.0
                    x1+=0.2
                    x2+=0.2

                ax.add_artist(
                    Arrow3D(
                        [x1, x2], [y1, y2], [z1, z2],
                        mutation_scale=15, lw=prob * 5,
                        arrowstyle="-|>", color="black", alpha=0.8))

    # Set plot limits and labels
    ax.set_xlim(x_min-0.4, x_max+0.4)
    ax.set_ylim(-0.5, 0.3)
    ax.set_zlim(-0.0, len(layer_heights)-1)

    ax.set_xlabel("Cell Types")
    ax.set_ylabel("")
    ax.set_zlabel("Layer Height")

    ax.set_title("Layer Transition Probabilities in Multilayer Graph (Adjusted View)")

    # Generate x-ticks and labels for cell types
    x_dict = {
        celltype: x for celltype, x in x_positions.items()
        if celltype != cell_communication_layer_name}
    ax.set_xticks(list(x_dict.values()))
    ax.set_xticklabels(list(x_dict.keys()))

    # Generate z-ticks and labels for layer types
    z_dict = {layer_type: z for layer_type, z in layer_heights.items()}
    ax.set_zticks(list(z_dict.values()))
    ax.set_zticklabels(list(z_dict.keys()))

    # Hide y-axis
    ax.set_yticks([])  # Hide x-axis ticks
    ax.set_yticklabels([])  # Hide x-axis tick labels

    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = (tmp_planes[3], tmp_planes[2], 
                        tmp_planes[5], tmp_planes[4], 
                        tmp_planes[0], tmp_planes[1])

    if not display_layer_axis:
        ax.set_zticks([])  # Hide z-axis ticks
        ax.set_zticklabels([])  # Hide z-axis tick labels
        ax.set_zlabel("")
        ax.set_xticks([])  # Hide x-axis ticks
        ax.set_xticklabels([])  # Hide x-axis tick labels
        ax.set_xlabel("")

    # Hide the grid and axis lines
    ax.grid(False)
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)

    # Hide axis lines while keeping ticks
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide x-axis line
    ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide y-axis line
    ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide z-axis line

    ax.set_box_aspect((np.sqrt(len(lamb)),1,1.5))
plt.show()


class Celltype:
    def __init__(
        self,
        receptor_graph:Union[str, pd.DataFrame],
        grn_graph:Union[str, pd.DataFrame],
        bipartite:Union[str, pd.DataFrame],
        celltype_name:str,
        lamb=None,
        eta=None,
        receptor_graph_directed:bool=True,
        receptor_graph_weighted:bool=False,
        grn_graph_directed:bool=False,
        grn_graph_weighted:bool=True,
        bipartite_graph_directed:bool=True,
        bipartite_graph_weighted:bool=False,
        copy_graphs=True,
        seeds:Union[list, dict]=[]
    ):
        """
        Create a Celltype object to explore the multilayer graph of a celltype.
        It can be used to predict the targets/regulators of a celltype based on
        the expression of its genes and the ligands it produces.

        It can be used in a Multicell object to explore the interactions between
        different celltypes.

        Parameters
        ----------
        receptor_graph : Union[str, pd.DataFrame]
            The graph of the receptors of the celltype.
            It can be a string with the path to a file containing the graph or a pandas DataFrame.
        grn_graph : Union[str, pd.DataFrame]
            The graph of the genes of the celltype.
            It can be a string with the path to a file containing the graph or a pandas DataFrame.
        bipartite : Union[str, pd.DataFrame]
            The bipartite graph between the receptors and the genes of the celltype.
            It can be a string with the path to a file containing the graph or a pandas DataFrame.
        celltype_name : str
            The name of the celltype.
        lamb : pd.DataFrame, optional
            The transition matrix between the layers of the celltype.
            If None, it will be set to a default one allowing transition between receptor,
            and receptor and grn layer exploration 
        eta : pd.Series, optional
            The probability of restarting at each layer.
            If None, it will be set to a vector with the same probability for each layer.
        receptor_graph_directed : bool, optional
            If True, the receptor graph is directed.
            The default is True.
        receptor_graph_weighted : bool, optional
            If True, the receptor graph is weighted.
            The default is False.
        grn_graph_directed : bool, optional
            If True, the grn graph is directed.
            The default is False.
        grn_graph_weighted : bool, optional 
            If True, the grn graph is weighted.
            The default is True.
        bipartite_graph_directed : bool, optional
            If True, the bipartite graph is directed.
            The default is True.
        bipartite_graph_weighted : bool, optional
            If True, the bipartite graph is weighted.
            The default is False.
        copy_graphs : bool, optional
            If True, the graphs are copied.
            The default is True.
        seeds : Union[list, dict], optional
            The seeds to explore the multilayer graph.
            If a list, the seeds are the nodes to explore.
            If a dictionary, the seeds are the nodes with their weights.
            The default is [].
        
        Returns
        -------
        Celltype
        """
        
        self.seeds = seeds
        self.celltype_name = celltype_name
        
        # Copy networks, True if no problem of memory
        if copy_graphs:
            if isinstance(receptor_graph, pd.DataFrame):
                receptor_graph = receptor_graph.copy()
            if isinstance(grn_graph, pd.DataFrame):
                grn_graph = grn_graph.copy()
            if isinstance(bipartite, pd.DataFrame):
                bipartite = bipartite.copy()

        # Type of graphs
        receptor_graph_type = str(int(receptor_graph_directed)) + str(int(receptor_graph_weighted))
        grn_graph_type = str(int(grn_graph_directed)) + str(int(grn_graph_weighted))
        bipartite_graph_type = str(int(bipartite_graph_directed)) + str(int(bipartite_graph_weighted))
        
        # Format multiplex dictionary
        self.multiplexes = {
            celltype_name + "_receptor": {
                "names":["receptor"],
                "graph_type":[receptor_graph_type],
                "layers":[receptor_graph]
            },
            celltype_name + "_grn": {
                "names":["grn"],
                "graph_type":[grn_graph_type], 
                "layers":[grn_graph]
            }
        }

        # Format bipartite dictionary
        bipartite = bipartite.rename({
            'receptor': 'col2',
            'grn': 'col1',
            'score': 'weight'}, axis=1)

        self.bipartites = {
            celltype_name+"_grn-" + celltype_name + "_receptor": {
                "source": celltype_name + "_receptor",
                "target": celltype_name + "_grn",
                "graph_type": bipartite_graph_type,
                "edge_list_df": bipartite
            }
        }
        
        # Prepare lamb
        if lamb is None:    
            lamb = pd.DataFrame(
                np.zeros((len(self.multiplexes), len(self.multiplexes))),
                index=list(self.multiplexes.keys()),
                columns=list(self.multiplexes.keys())
            )
            lamb.loc[
                self.celltype_name + "_grn",
                self.celltype_name + "_grn"] = 1
            lamb.loc[
                self.celltype_name + "_receptor", self.celltype_name + "_grn"] = 1
        self.lamb = lamb
        
        # Fake eta
        if eta is None:
            self.eta = pd.Series(np.ones(2)/2, index=self.multiplexes.keys())
        else:
            self.eta = eta

    def Multixrank(
        self,
        self_loops=True,
        restart_proba=0.7,
        verbose=True):

        if self.seeds==[]:
            warnings.warn("""
            No seeds provided to explore the multilayer, all the scoree will be null (np.nan).
            You can pass it as an argument earlier or set up the .seeds attribute.
            """)
        # if seeds are provided as a dictionary
        if isinstance(self.seeds, dict):
            # Every value should be a numeric value (int or float)
            if not all(
                [isinstance(v, int) or isinstance(v, float)
                for v in self.seeds.values()]):
                raise ValueError("Seeds should be a dictionary with numeric values (weight per seed).")

            # Create the multixrank object
            multilayer = hummuspy.create_multilayer.Multixrank(
                multiplex=self.multiplexes,
                bipartite=self.bipartites,
                lamb=self.lamb.values.T,
                seeds=[],
                self_loops=self_loops,
                restart_proba=restart_proba,
                eta = self.eta.tolist(),
                verbose=verbose)

            # Create a probability vector from the seeds' names and weights
            node_list = [node for node_list in  multilayer.multiplexall_node_list2d for node in node_list]
            prox_vector  = np.zeros(len(node_list))

            # Position is matching node order in the multilayer.
            for seed in self.seeds:
                idx = np.where(np.array(node_list)==seed)[0]
                if len(idx)!=0:
                    prox_vector[idx[0]] = self.seeds[seed]

            # Values are normalized.
            multilayer.pr = prox_vector/prox_vector.sum()

        else:
            #should be a list or an array
            if not isinstance(self.seeds, (list, np.ndarray)):
                raise ValueError("Seeds should be a list or an array.")

            # Create the multixrank object
            multilayer = hummuspy.create_multilayer.Multixrank(
                multiplex=self.multiplexes,
                bipartite=self.bipartites,
                lamb=self.lamb.values.T,
                seeds=self.seeds,
                self_loops=self_loops,
                restart_proba=restart_proba,
                eta = self.eta.tolist(),
                verbose=verbose)

        return multilayer


    def rename_celltype(
        self,
        celltype_name
    ):

        self.multiplexes = {k.replace(self.celltype_name, celltype_name, 1): v
                            for k, v in self.multiplexes.items()
                            }

        self.bipartites = {
            celltype_name+"_grn-" + celltype_name + "_receptor": self.bipartites[
                self.celltype_name+"_grn-" + self.celltype_name + "_receptor"]\
        }
        for bipartite in self.bipartites:
            self.bipartites[bipartite]["source"] = self.bipartites[bipartite]["source"].replace(
                self.celltype_name, celltype_name)
            self.bipartites[bipartite]["target"] = self.bipartites[bipartite]["target"].replace(
                self.celltype_name, celltype_name)            
        
        lamb = self.lamb
        lamb.index = lamb.index.str.replace(self.celltype_name+"_", celltype_name+"_")
        lamb.columns = lamb.columns.str.replace(self.celltype_name+"_", celltype_name+"_")
        self.lamb = lamb

        eta = self.eta
        eta.index = eta.index.str.replace(self.celltype_name+"_", celltype_name+"_")
        self.eta = eta

        self.celltype_name = celltype_name

            
class Multicell(Celltype):
    def __init__(
        self,
        celltypes:Union[list[Celltype], list[dict], dict[dict], dict[Celltype]],
        cell_communication_graph:pd.DataFrame,
        lamb=None,
        eta=None,
        cell_communication_graph_directed:bool=True,
        cell_communication_graph_weighted:bool=True,
        # bipartite parameters can be -1, 0, 1 here
        bipartite_grn_cell_communication_directed:Union[bool, int]=True,
        bipartite_grn_cell_communication_weighted:bool=False,
        bipartite_cell_communication_receptor_directed:bool=True,
        bipartite_cell_communication_receptor_weighted:bool=False,
        copy_graphs=True,
        seeds = [],
        verbose=True
    ):

        if not isinstance(bipartite_grn_cell_communication_directed, (bool, int)):
            raise ValueError("bipartite_grn_cell_communication_directed must be True, False, 0, 1, or -1.")
        if not isinstance(bipartite_cell_communication_receptor_directed, (bool, int)):
            raise ValueError("bipartite_cell_communication_receptor_directed must be True, False, 0, 1, or -1.")
        if not isinstance(celltypes, (list, dict)):
            raise ValueError("celltypes must be a list or a dictionary.")
        if not isinstance(bipartite_grn_cell_communication_weighted, (bool, int)):
            raise ValueError("bipartite_grn_cell_communication_weighted must be True, False, 0, or 1.")
        if not isinstance(bipartite_cell_communication_receptor_weighted, (bool, int)):
            raise ValueError("bipartite_cell_communication_receptor_weighted must be True, False, 0, or 1.")
        if not isinstance(cell_communication_graph, pd.DataFrame):
            raise ValueError("cell_communication_graph must be a pandas DataFrame.")
        if not isinstance(cell_communication_graph_directed, (bool, int)):
            raise ValueError("cell_communication_graph_directed must be a boolean.")
        if cell_communication_graph_weighted not in [True, False, 1, 0]:
            raise ValueError("cell_communication_graph_weighted must be a boolean.")

        # Check if celltypes is a dictionary and convert it to a list of Celltype objects
        # If celltypes is a dictionary, the keys will be the new names of the celltypes
        self.celltypes_names = None
        if isinstance(celltypes, dict):
            self.celltypes_names = list(celltypes.keys())
            celltypes = list(celltypes.values())
            if verbose:
                warnings.warn("The celltypes dictionary was converted to a list of Celltype objects.\n" +
                              "The keys of the dictionary will be the celltype names.")
        else:
            self.celltypes_names = [
                celltype.celltype_name if isinstance(celltype, Celltype)
                else celltype["celltype_name"] for celltype in celltypes
                ]
        # Loop over celltypes and create Celltype objects
        for i in range(len(celltypes)):
            if isinstance(celltypes[i], dict):
                celltype_dict = celltypes[i]
                if "lamb" not in celltypes[i].keys():
                    celltype_dict["lamb"]=None
                celltypes[i] = Celltype(**celltype_dict)
            elif isinstance(celltypes[i], Celltype):
                if self.celltypes_names!=celltypes[i].celltype_name:
                    celltypes[i].rename_celltype(self.celltypes_names[i])
            elif not isinstance(celltypes[i], Celltype):
                raise ValueError("celltypes must be a list of Celltype objects or dictionaries.")

        # Rename celltypes if celltype_new_names is not None    
        if self.celltypes_names is not None:
            for i, celltype in enumerate(celltypes):
                celltype.celltype_name = self.celltypes_names[i]

        if copy_graphs: #test if when false not all copied already
            celltypes = list(celltypes)
            cell_communication_graph = cell_communication_graph.copy()

        # Store celltypes in a dictionary
        celltypes = {celltype.celltype_name: celltype for celltype in celltypes}

        # Save general parameters
        self.celltypes_names = list(celltypes.keys())
        self.seeds = seeds

        # Create bipartites between cell communication and celltype layer
        self.bipartites = {}
        for celltype in self.celltypes_names:
            mask_source_cell = cell_communication_graph["celltype_source"] == celltype
            mask_target_cell = cell_communication_graph["celltype_target"] == celltype
            
            # Create bipartites for each celltype
            # Link genes of each celltype to ligands of this celltype in cell communication
            bipartite_ligand = pd.DataFrame({
                "col2": cell_communication_graph[mask_source_cell]["source"].unique() + '-' + celltype,
                "col1": cell_communication_graph[mask_source_cell]["source"].unique(),
            })
            self.bipartites[celltype + "_to_ligands"] = {
                "source": "cell_communication",
                "target": celltype + "_grn",
                "graph_type": str(int(bool(bipartite_grn_cell_communication_directed)))\
                    + str(int(bipartite_grn_cell_communication_weighted)),
                "edge_list_df": bipartite_ligand
            }
            # Link receptors of each celltype to genes of this celltype in cell communication
            bipartite_receptor = pd.DataFrame({
                "col1": cell_communication_graph[mask_target_cell][
                    "target"].unique()+"_receptor",
                "col2": cell_communication_graph[mask_target_cell][
                    "target"].unique() + '-' + celltype
            })
            self.bipartites[celltype + "_to_receptors"] = {
                "source": "cell_communication",
                "target": celltype + "_receptor",
                "graph_type": str(int(bool(bipartite_cell_communication_receptor_directed)))\
                    + str(int(bipartite_cell_communication_receptor_weighted)),
                "edge_list_df": bipartite_receptor
            }

        # Create cell communication layer
        cell_communication_graph_type = str(int(cell_communication_graph_directed)) + str(int(cell_communication_graph_weighted))
        cell_communication_graph["source"] = cell_communication_graph["source"] + '-' + cell_communication_graph["celltype_source"]
        cell_communication_graph["target"] = cell_communication_graph["target"] + '-' + cell_communication_graph["celltype_target"]
        cell_communication_graph.rename(columns={"lr_means": "weight"}, inplace=True)

        self.multiplexes = {
            "cell_communication": {
                "names":["cell_communication"],
                "graph_type":[cell_communication_graph_type],
                "layers":[cell_communication_graph]
            }
        }

        # Update with layers and bipartites of each celltype
        for celltype in self.celltypes_names:
            self.multiplexes.update(celltypes[celltype].multiplexes)
            self.bipartites.update(celltypes[celltype].bipartites)
        
        # Update nodes of each celtype to add a celltype specific suffix
        for celltype in self.celltypes_names:
            for multiplex in self.multiplexes:
                if celltype in multiplex:
                    self.multiplexes[multiplex]["layers"][0]["source"] = \
                        self.multiplexes[multiplex]["layers"][0]["source"] + \
                            "::" + celltype
                    self.multiplexes[multiplex]["layers"][0]["target"] = \
                        self.multiplexes[multiplex]["layers"][0]["target"] + \
                            "::" + celltype
            for bipartite in self.bipartites:
                if celltype in self.bipartites[bipartite]["source"]:
                    self.bipartites[bipartite]["edge_list_df"]["col2"] = \
                        self.bipartites[bipartite]["edge_list_df"]["col2"] + \
                            "::" + celltype
                if celltype in self.bipartites[bipartite]["target"]:
                    self.bipartites[bipartite]["edge_list_df"]["col1"] = \
                        self.bipartites[bipartite]["edge_list_df"]["col1"] + \
                            "::" + celltype

        # Inverse source and target to match multixrank error
        for multiplex in self.multiplexes:
            self.multiplexes[multiplex]["layers"][0].rename(
                columns={"source": "target", "target": "source"}, inplace=True)

        # Prepare lamb matrix if not provided
        if lamb is None:
            lamb = pd.DataFrame(
                np.zeros((len(self.multiplexes), len(self.multiplexes))),
                index=list(self.multiplexes.keys()),
                columns=list(self.multiplexes.keys())
            )

            grn_layer_mask = lamb.index.str.endswith("_grn")
            receptor_layer_mask = lamb.index.str.endswith("_receptor")
            ccc_layer_mask = lamb.index.str.endswith("cell_communication")

            for celltype in self.celltypes_names:
                cell_type_mask = lamb.index.str.contains(celltype)

                lamb.loc[
                    grn_layer_mask*cell_type_mask,
                    receptor_layer_mask*cell_type_mask
                    + grn_layer_mask*cell_type_mask] = 1
                lamb.loc[
                    ccc_layer_mask,
                    ccc_layer_mask
                    + grn_layer_mask*cell_type_mask] += 1
                lamb.loc[
                    receptor_layer_mask*cell_type_mask,
                    ccc_layer_mask] = 1
            lamb = lamb.transpose().div(lamb.transpose().sum(1), 0)
        self.lamb = lamb

        # Fake eta
        if eta is None:
            self.eta = pd.Series(
                np.ones(len(self.multiplexes))/len(self.multiplexes),
                index=self.multiplexes.keys())
        else:
            self.eta = eta


    def illustrate_multicell(
        self,
        figsize=(12, 10),
        azim = 90,
        elev = 20,
        display_layer_axis=True,
        display_self_proba=True,
        display_layer_names=False,
        alpha_layers = 0.5,
        cell_communication_layer_name="cell_communication"
        ):

        lamb = self.lamb

        illustrate_multicell(
            lamb=lamb,
            figsize=figsize,
            azim=azim,
            elev=elev,
            display_layer_axis=display_layer_axis,
            display_self_proba=display_self_proba,
            display_layer_names=display_layer_names,
            alpha_layers=alpha_layers,
            cell_communication_layer_name=cell_communication_layer_name
        )


def set_lambda(
    multicell=None,
    celltypes=None,
    cell_communication_layer_name="cell_communication",
    direction: Union["upstream", "downstream"]="downstream",
    strategy: Union["intracell", "intercell"]="intercell"
):
    if multicell is None and celltypes is None:
        raise ValueError("Either multicell or celltype should not be None.")
    elif multicell is not None:
        if celltypes is not None:
            raise warnings.warn("Both multicell and celltypes are provided."+
                "multicell will be used.")
        else:
            multiplexes = list(multicell.multiplexes.keys())
            celltypes = multicell.celltypes_names
    else:
        for celltype in celltypes:
            multiplexes.append(f"{celltype}_receptor")
            multiplexes.append(f"{celltype}_grn")
        multiplexes.append(cell_communication_layer_name)

    lamb = pd.DataFrame(
        np.zeros((len(multiplexes), len(multiplexes))),
        index=multiplexes,
        columns=multiplexes
    )
    
    # define masks
    grn_layer_mask = lamb.index.str.endswith("_grn")
    receptor_layer_mask = lamb.index.str.endswith("_receptor")
    ccc_layer_mask = lamb.index.str.endswith(cell_communication_layer_name)

    for celltype in celltypes:
        cell_type_mask = lamb.index.str.contains(celltype)

        # Transition FROM receptor and grn TO grn
        lamb.loc[
            receptor_layer_mask*cell_type_mask
            + grn_layer_mask*cell_type_mask,
            grn_layer_mask*cell_type_mask
            ] = 1
        
    # Transition FROM cell communication TO cell communication and receptor
    lamb.loc[
        ccc_layer_mask,
        ccc_layer_mask + receptor_layer_mask,
        ] = 1

    if strategy=="intercell":
        # Transition FROM grn TO cell communication
        lamb.loc[
            grn_layer_mask,
            ccc_layer_mask
            ] = 1

    # transpose if direction is upstream
    if direction=="upstream":
        lamb = lamb.transpose()

    lamb = lamb.div(lamb.sum(1), 0)

    return lamb



# To find targets : 

# Run intracellular predictions
## --> Direct targets

# Run intercellular predictions from prodcuted ligands
# (Corresponding to intracellular genes in the CCC)
## --> Indirect targets

def multicell_targets(
    seeds,
    celltypes,
    grn,
    receptor_layer,
    receptor_grn,
    ccc,
    grn_graph_directed=False,
    grn_graph_weighted=True,
    bipartite_graph_directed=False,
    bipartite_graph_weighted=True,
    receptor_graph_directed=False,
    receptor_graph_weighted=False,
    cell_communication_graph_directed=False,
    cell_communication_graph_weighted=True,
    bipartite_grn_cell_communication_directed=False,
    bipartite_grn_cell_communication_weighted=False,
    bipartite_cell_communication_receptor_directed=False,
    bipartite_cell_communication_receptor_weighted=False,
    restart_proba=0.6,
    extend_seeds=False
):
    if extend_seeds:
        if type(seeds) is not list and type(seeds) is not dict:
            seeds = list(seeds)
        starting_nodes = [f"{seed}-{celltype}" for seed in seeds for celltype in celltypes]
    else:
        starting_nodes = seeds

    generic_multicell = Multicell(
        celltypes={
            celltype: Celltype(
                receptor_graph=receptor_layer,
                grn_graph=grn,
                bipartite=receptor_grn,
                celltype_name=celltype,
                receptor_graph_directed=receptor_graph_directed,
                receptor_graph_weighted=receptor_graph_weighted,
                grn_graph_directed=grn_graph_directed,
                grn_graph_weighted=grn_graph_weighted,
                bipartite_graph_directed=bipartite_graph_directed,
                bipartite_graph_weighted=bipartite_graph_weighted,
            )
            for celltype in celltypes
        },
        cell_communication_graph=ccc,
        cell_communication_graph_directed=cell_communication_graph_directed,
        cell_communication_graph_weighted=cell_communication_graph_weighted,
        bipartite_grn_cell_communication_directed=bipartite_grn_cell_communication_directed,
        bipartite_grn_cell_communication_weighted=bipartite_grn_cell_communication_weighted,
        bipartite_cell_communication_receptor_directed=bipartite_cell_communication_receptor_directed,
        bipartite_cell_communication_receptor_weighted=bipartite_cell_communication_receptor_weighted,
        seeds=starting_nodes,
    )

    # Intracellular direct targets
    generic_multicell.lamb = set_lambda(
        generic_multicell, direction="downstream", strategy="intracell"
    )

    multilayer = generic_multicell.Multixrank(restart_proba=restart_proba)
    intracell = multilayer.random_walk_rank()
    intracell = intracell.sort_values(ascending=False, by="node")
    intracell = intracell[intracell["layer"] == "grn"].set_index("node")
    # Extracellular indirect regulations
    extracell = {}
    generic_multicell.lamb = set_lambda(
        generic_multicell, direction="downstream", strategy="intercell"
    )

    for celltype in celltypes:
        cell_seeds = intracell[intracell["multiplex"] == f"{celltype}_grn"]["score"]
        cell_seeds.index = cell_seeds.index.str.replace("::", "-")
        cell_seeds = cell_seeds[
            cell_seeds.index.isin(
                generic_multicell.multiplexes["cell_communication"]["layers"][0]["target"].unique()
            )
        ]
        generic_multicell.seeds = cell_seeds.to_dict()

        multilayer = generic_multicell.Multixrank(restart_proba=restart_proba)

        extracell[celltype] = (
            multilayer.random_walk_rank()
            .sort_values("node", ascending=False)
            .query("layer == 'grn'")
            .set_index("node")[["score", "multiplex"]]
        )
        extracell[celltype]["score"] = extracell[celltype]["score"] / \
            extracell[celltype]["score"].sum() * cell_seeds.sum()

    # Store cell type contributions
    cell_contributions = {}

    for celltype_target in celltypes:
        cell_contributions[celltype_target] = intracell[
            intracell["multiplex"] == f"{celltype_target}_grn"
        ]["score"].to_frame()
        cell_contributions[celltype_target].columns = [celltype_target+"_direct"]

        for celltype_source in celltypes:
            cell_contributions[celltype_target][celltype_source] = extracell[celltype_source][
                extracell[celltype_source]["multiplex"] == f"{celltype_target}_grn"
            ]["score"]

    return cell_contributions