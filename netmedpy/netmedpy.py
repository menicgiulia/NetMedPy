# -*- coding: utf-8 -*-
"""
Functionality
----------------

This module provides the main NetMedPy functions for Network Medicine and
network topology analysis. Functions allow for different notions of distance
between nodes and null models.

The main functions in this module include:

- `extract_lcc`: Computes the Largest Connected Component (LCC) for a specified subset of nodes within a graph.
- `lcc_significance`: Calculates the statistical significance of the LCC, as defined by a subset of nodes in a graph, according to a specified null model.
- `all_pair_distances`: Calculates distances between every pair of nodes in a graph.
- `save_distances`: Saves a precomputed `DistanceMatrix` to a compressed NumPy `.npz` file.
- `load_distances`: Loads a precomputed `DistanceMatrix` from a compressed NumPy `.npz` file or from a legacy pickle file.
- `get_amspl`: Calculates the Average Minimum Shortest Path Length between nodes.
- `proximity`: Computes the proximity between two sets of nodes in a graph.
- `separation`: Calculates the separation between two sets of nodes in a network.
- `separation_z_score`: Determines the z-score of the separation between two node sets based on randomized samples.
- `screening`: Screens for proximity/separation between sets of source and target nodes.

Distance metrics.
------------------

When calculating the distance matrix, five information flow metrics are available to the user:

- `shortest_path`: Distance between nodes is based on the length of the path with the least number of edges or lowest total weight that connects two nodes.
- `random_walk`: Distance between nodes is based on the probability of reaching one node from another via a random walk with restart.
- `biased_random_walk`: Same as `random_walk` but compensating for the bias induced by the degree of the target node.
- `communicability`: Distance between nodes is based on the concept of communicability, defined as the ability of nodes to communicate through all available paths in a network, considering both indirect and direct connections.
- `custom`: Distance between nodes is computed with a user-provided function.

Null models.
-------------

The NetMedPy functions involving statistical analysis allow the user to select among the following null models:

- `degree_match`: Selects random samples replicating the original node-set degree distribution.
- `log_binning`: Categorizes node degrees into logarithmically sized bins and samples from matching bins.
- `strength_binning`: Analogous to `log_binning`, using node strength instead of degree.
- `uniform`: Randomly selects nodes from the entire network, disregarding degree or strength.
- `custom`: Allows users to specify custom node proxies through `node_bucket`.

These functions use both exact and approximate methods for node-set randomization.
Precomputed distance matrices are leveraged for efficient computation, and some
operations use multiprocessing with shared memory for scalability.

Required packages:
-------------------
    - networkx
    - numpy
    - pickle
    - multiprocessing
    - random
    - scipy
    - pandas
"""

import networkx as nx
import numpy as np
import pickle
import random
import warnings
import os
import pandas as pd

import multiprocessing as mp
from multiprocessing import cpu_count, get_context
from multiprocessing.shared_memory import SharedMemory
from concurrent.futures import ProcessPoolExecutor

class DistanceMatrix:
    """A class to manage a square matrix representing the distances between nodes.
    
    Attributes
    ----------
    nodes : list
        Identifiers for all nodes within the matrix.
    node_to_idx : dict
        Maps node identifiers to their respective indices in the matrix.
    matrix : np.ndarray
        A 2D numpy array where each element [i][j] is the distance from node i to node j,
        initially set to infinity.
    """
    
    @classmethod
    def from_components(cls, nodes, node_to_idx, matrix):
        obj = cls()
        
        # Validation
        assert len(nodes) == matrix.shape[0] == matrix.shape[1]
        assert all(node_to_idx[n] == i for i, n in enumerate(nodes))
        
        obj.nodes = nodes
        obj.node_to_idx = node_to_idx
        obj.matrix = matrix
        
        return obj

    def _from_name_list(self, nodes):
        """Initializes the matrix using a list of node identifiers, setting all distances to infinity.
        
        Parameters
        ----------
        nodes : list
            List of identifiers for the nodes.
        
        Notes
        -----
        This method is typically used when node connections are unknown or not yet defined,
        allowing for subsequent updates with actual distances.
        """
        self.nodes = nodes
        self.node_to_idx = {node: idx for idx, node in enumerate(self.nodes)}
        self.matrix = np.full((len(self.nodes), len(self.nodes)), np.inf)

    def _from_dictionary(self, distance_dict):
        """Initializes the distance matrix from a nested dictionary that specifies the distances
        between node pairs.
        
        Parameters
        ----------
        distance_dict : dict
            A dictionary where the keys are node identifiers and the values are dictionaries
            that map connected node identifiers to distances.
        
        Notes
        -----
        This method populates the matrix with known distances, ideal for networks with predefined
        connectivity.
        """
        self.nodes = list(distance_dict.keys())
        self.node_to_idx = {node: idx for idx, node in enumerate(self.nodes)}
        self.matrix = np.full((len(self.nodes), len(self.nodes)), np.inf)

        for node_i, neighbors in distance_dict.items():
            for node_j, distance in neighbors.items():
                self.matrix[self.node_to_idx[node_i], self.node_to_idx[node_j]] = distance

    def get(self, A, B):
        """Fetches the distance between two specified nodes.
        
        Parameters
        ----------
        A : str
            Identifier for the source node.
        B : str
            Identifier for the target node.
        
        Returns
        -------
        float
            The distance between the source and target nodes.
        """
        idx_A = self.node_to_idx[A]
        idx_B = self.node_to_idx[B]
        return self.matrix[idx_A, idx_B]

    def put(self, A, B, d):
        """Sets the distance between two specified nodes.
        
        Parameters
        ----------
        A : str
            Identifier for the source node.
        B : str
            Identifier for the target node.
        d : float
            Distance to set between the source and target nodes.
        """
        idx_A = self.node_to_idx[A]
        idx_B = self.node_to_idx[B]
        self.matrix[idx_A, idx_B] = d

    def get_matrix(self):
        """Returns the full distance matrix.
        
        Returns
        -------
        np.ndarray
            The complete matrix, showing distances between all node pairs.
        """
        return self.matrix

_SHM_MATRIX = None
_SHM_ARRAY = None
_GRAPH = None
_NODE_TO_IDX = None
_CUSTOM_DISTANCE = None
_CUSTOM_KWARGS = None

_SCREEN_NETWORK = None
_SCREEN_SOURCES = None
_SCREEN_TARGETS = None
_SCREEN_NODE_BUCKET = None
_SCREEN_DISTANCE_MATRIX = None

_SCREEN_SHM = None
_SCREEN_MATRIX = None

def _init_distance_worker(shm_name, shape, dtype_str, graph, node_to_idx,
                          custom_distance=None, custom_kwargs=None):
    global _SHM_MATRIX, _SHM_ARRAY, _GRAPH, _NODE_TO_IDX, _CUSTOM_DISTANCE, _CUSTOM_KWARGS

    _SHM_MATRIX = SharedMemory(name=shm_name)
    dtype = np.dtype(dtype_str)
    _SHM_ARRAY = np.ndarray(shape, dtype=dtype, buffer=_SHM_MATRIX.buf)

    _GRAPH = graph
    _NODE_TO_IDX = node_to_idx
    _CUSTOM_DISTANCE = custom_distance
    _CUSTOM_KWARGS = {} if custom_kwargs is None else custom_kwargs

def _split_nodes_with_row_ranges(nodes, n_chunks):
    nodes = list(nodes)
    k, m = divmod(len(nodes), n_chunks)

    out = []
    start = 0
    for i in range(n_chunks):
        size = k + (1 if i < m else 0)
        end = start + size
        if start < end:
            out.append((start, end, nodes[start:end]))
        start = end
    return out



def _split_into_chunks(lst,n):
    k, m = divmod(len(lst), n)
    return (lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def _get_amspl_all_operations(net,A,B):
    net_nodes = set(net.nodes)
    sA = set(A)
    sB = set(B)

    #Determine the nodes in A and B that are also in the network
    valid_a = sA & net_nodes
    valid_b = sB & net_nodes

    min_distances = np.zeros(len(valid_a))
    idx = 0
    for n in valid_a:
        min_dist = float('inf')
        for m in valid_b:
            db = nx.shortest_path_length(net,n,m)
            if(db is not None):
                min_dist = min(min_dist,db)
        min_distances[idx]=min_dist
        idx+=1
    avg = min_distances.mean()
    return avg


def _get_amspl_dmatrix(A,B,D):
    min_distances = np.zeros(len(A))
    idx = 0
    for n in A:
        min_dist = float('inf')
        for m in B:
            db = D.get(n,m)
            if(db is not None):
                min_dist = min(min_dist,db)
        min_distances[idx]=min_dist
        idx+=1
    avg = min_distances.mean()
    return avg


def get_amspl(net, A, B, D):
    """
    Compute the average minimum shortest path length (AMSPL) from a set of
    source nodes `A` to a set of target nodes `B` using a precomputed
    `DistanceMatrix`.

    For each valid node in `A`, the function finds its minimum distance to any
    valid node in `B`, then averages these minima.

    Parameters
    ----------
    net : networkx.Graph
        Input graph. Only nodes present in the graph are considered.

    A : iterable
        Source node set.

    B : iterable
        Target node set.

    D : DistanceMatrix
        Precomputed distance matrix.

    Returns
    -------
    float
        Average minimum distance from nodes in `A` to nodes in `B`.
    """
    net_nodes = set(net.nodes)
    sA = set(A)
    sB = set(B)

    #Determine the nodes in A and B that are also in the network
    valid_a = sA & net_nodes
    valid_b = sB & net_nodes

    return _get_amspl_dmatrix(valid_a,valid_b,D)


def _single_shortest_path_chunk(task):
    row_start, row_end, source_nodes = task

    graph = _GRAPH
    node_to_idx = _NODE_TO_IDX
    out = _SHM_ARRAY

    weighted = bool(nx.get_edge_attributes(graph, "weight"))

    for local_i, s in enumerate(source_nodes):
        row_idx = row_start + local_i
        row = out[row_idx]
        row.fill(np.inf)

        if weighted:
            d = nx.shortest_path_length(graph, s, weight="weight")
        else:
            d = nx.shortest_path_length(graph, s)

        for target, distance in d.items():
            col_idx = node_to_idx[target]
            row[col_idx] = distance

    return (row_start, row_end)

def _custom_distance_chunk(task):
    row_start, row_end, source_nodes = task

    graph = _GRAPH
    node_to_idx = _NODE_TO_IDX
    out = _SHM_ARRAY
    distance_fn = _CUSTOM_DISTANCE
    kwargs = _CUSTOM_KWARGS

    for local_i, s in enumerate(source_nodes):
        row_idx = row_start + local_i
        row = out[row_idx]
        row.fill(np.inf)

        res_dict = distance_fn(s, graph, **kwargs)
        for target, d in res_dict.items():
            col_idx = node_to_idx[target]
            row[col_idx] = d

    return (row_start, row_end)

def _run_shared_matrix_jobs(graph, node_to_idx, worker_fn, num_cpus, n_tasks,
                            dtype=np.float32, custom_distance=None, custom_kwargs=None):
    n = len(graph)
    shape = (n, n)
    dtype = np.dtype(dtype)

    shm = SharedMemory(create=True, size=int(np.prod(shape)) * dtype.itemsize)
    result = np.ndarray(shape, dtype=dtype, buffer=shm.buf)
    result.fill(np.inf)

    tasks = _split_nodes_with_row_ranges(list(graph.nodes), n_tasks)

    ctx = get_context("spawn")

    try:
        with ProcessPoolExecutor(
            max_workers=num_cpus,
            mp_context=ctx,
            initializer=_init_distance_worker,
            initargs=(
                shm.name,
                shape,
                dtype.str,
                graph,
                node_to_idx,
                custom_distance,
                custom_kwargs,
            ),
        ) as ex:
            list(ex.map(worker_fn, tasks, chunksize=1))

        final = result.copy()
    finally:
        shm.close()
        shm.unlink()

    return final


def _spl_distance(graph, node_to_idx, num_cpus, n_tasks, dtype=np.float32):
    return _run_shared_matrix_jobs(
        graph=graph,
        node_to_idx=node_to_idx,
        worker_fn=_single_shortest_path_chunk,
        num_cpus=num_cpus,
        n_tasks=n_tasks,
        dtype=dtype,
    )



def _RWR_Matrix(G, c):
    # Map from node names to indices and vice versa
    node_to_index = {node: i for i, node in enumerate(G.nodes())}
    index_to_node = {i: node for node, i in node_to_index.items()}

    A = nx.to_numpy_array(G)

    # Normalize A to get transition probabilities
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1  # Avoid division by zero for isolated nodes
    T = A / col_sums

    # Calculate the modified transition matrix for RWR
    n = G.number_of_nodes()
    I = np.eye(n)
    M = c * np.linalg.inv(I - (1 - c) * T)
    M = M.T

    return (M,node_to_index,index_to_node)

def _custom_all_distance(graph, node_to_idx, distance, num_cpus, n_tasks,
                         kwargs, dtype=np.float32):
    return _run_shared_matrix_jobs(
        graph=graph,
        node_to_idx=node_to_idx,
        worker_fn=_custom_distance_chunk,
        num_cpus=num_cpus,
        n_tasks=n_tasks,
        dtype=dtype,
        custom_distance=distance,
        custom_kwargs=kwargs,
    )


def all_pair_distances(graph, distance="shortest_path", custom_distance=None,
                       reset=0.2, n_processors=None, n_tasks=None,
                       matrix_dtype=np.float32, **kwargs):
    """
    Compute an all-pairs distance matrix for a connected graph.

    Depending on the selected `distance` mode, distances are computed from:
    shortest paths, random-walk-with-restart scores, biased random-walk scores,
    communicability, or a user-defined custom distance function.

    Parameters
    ----------
    graph : networkx.Graph
        Input graph. The graph must be connected.

    distance : {'shortest_path', 'random_walk', 'biased_random_walk',
                'communicability', 'custom'}, optional
        Distance definition to use. Default is `'shortest_path'`.

    custom_distance : callable, optional
        Required when `distance='custom'`. The function must accept
        `(source_node, graph, **kwargs)` and return a dictionary mapping
        target nodes to distances.

    reset : float, optional
        Restart probability used for random-walk-based distances.
        Must satisfy `0 <= reset <= 1`. Default is `0.2`.

    n_processors : int, optional
        Number of worker processes to use for parallel computations.
        If `None`, uses all available CPUs.

    n_tasks : int, optional
        Number of row chunks used to split the computation. Must be greater
        than or equal to the number of processors. If `None`, defaults to
        `n_processors`.

    matrix_dtype : numpy dtype, optional
        Data type used to store the resulting matrix. Default is `np.float32`.

    **kwargs
        Additional keyword arguments passed to `custom_distance` when
        `distance='custom'`.

    Returns
    -------
    DistanceMatrix
        A `DistanceMatrix` object containing node order, node-to-index mapping,
        and the computed matrix.

    Raises
    ------
    ValueError
        If the graph is disconnected, if an invalid distance mode is provided,
        if `reset` is outside `[0, 1]`, if `custom_distance` is missing in
        custom mode, or if `n_tasks < n_processors`.

    Notes
    -----
    For `'shortest_path'` and `'custom'`, the computation uses multiprocessing
    with shared memory to avoid duplicating the full matrix across workers.
    """
    if nx.number_connected_components(graph) > 1:
        raise ValueError("The network is not connected (it contains more than one connected component)")

    if distance not in ["shortest_path", "random_walk", "biased_random_walk", "communicability", "custom"]:
        raise ValueError("distance must be shortest_path|random_walk|biased_random_walk|communicability|custom")

    if distance in ["random_walk", "biased_random_walk"]:
        if reset < 0 or reset > 1:
            raise ValueError("Reset for random walks must comply 0 <= reset <= 1")

    num_cpus = cpu_count() if n_processors is None else n_processors
    n_tasks = num_cpus if n_tasks is None else n_tasks

    if n_tasks < num_cpus:
        raise ValueError("Number of tasks should be larger or equal the number of processors.")

    D = DistanceMatrix()
    D._from_name_list(list(graph.nodes()))
    D.matrix = None

    if distance == "shortest_path":
        D.matrix = _spl_distance(
            graph,
            D.node_to_idx,
            num_cpus=num_cpus,
            n_tasks=n_tasks,
            dtype=matrix_dtype,
        )

    elif distance == "random_walk":
        M, node_to_index, index_to_node = _RWR_Matrix(graph, reset)
        M = np.log(M)
        ma = np.max(M)
        mi = np.min(M)
        M = 1 - ((M - mi) / (ma - mi))

        D.nodes = list(graph.nodes())
        D.node_to_idx = node_to_index
        D.matrix = M.astype(matrix_dtype, copy=False)

    elif distance == "biased_random_walk":
        M, node_to_index, index_to_node = _RWR_Matrix(graph, reset)

        degs = np.array([graph.degree(index_to_node[i]) for i in range(len(graph))], dtype=np.float64)

        B = M / degs
        B = B.T
        B = B / B.sum(axis=0)
        B = B.T

        B = np.log(B)
        ma = np.max(B)
        mi = np.min(B)
        B = 1 - ((B - mi) / (ma - mi))

        D.nodes = list(graph.nodes())
        D.node_to_idx = node_to_index
        D.matrix = B.astype(matrix_dtype, copy=False)

    elif distance == "communicability":
        comm = nx.communicability_exp(graph)

        mat = np.full((len(D.nodes), len(D.nodes)), np.inf, dtype=np.float64)
        for a, dd in comm.items():
            i = D.node_to_idx[a]
            for b, c in dd.items():
                j = D.node_to_idx[b]
                mat[i, j] = c

        B = np.log(mat)
        ma = np.max(B)
        mi = np.min(B)
        B = 1 - ((B - mi) / (ma - mi))

        D.matrix = B.astype(matrix_dtype, copy=False)

    elif distance == "custom":
        if custom_distance is None:
            raise ValueError("custom_distance must be provided when distance='custom'")

        D.matrix = _custom_all_distance(
            graph,
            D.node_to_idx,
            distance=custom_distance,
            num_cpus=num_cpus,
            n_tasks=n_tasks,
            kwargs=kwargs,
            dtype=matrix_dtype,
        )

    return D



def save_distances(distances, filename):
    """
    Save a precomputed distance matrix to a compressed NumPy `.npz` file.

    This function stores the contents of a `DistanceMatrix` object in a compact
    archive that can later be reloaded with `load_distances`. The saved file
    contains the distance matrix, the node list, and the node-to-index mapping.

    Parameters
    ----------
    distances : DistanceMatrix
        Distance matrix object to save.

    filename : str
        Output file path. This is typically expected to use the `.npz` extension.
        If the file already exists, it will be overwritten.

    Notes
    -----
    The file is saved using `numpy.savez_compressed`, which usually produces
    much smaller files than pickle for dense numeric matrices.

    See Also
    --------
    load_distances
    """
    nodes = distances.nodes
    node_to_idx = distances.node_to_idx
    matrix = distances.matrix 

    n = max(node_to_idx.values()) + 1
    names_idx = [None] * n
        
    for name, idx in node_to_idx.items():
        names_idx[idx] = name

    np.savez_compressed(
        filename,
        matrix=matrix,
        nodes=nodes,
        node_to_idx = names_idx
    )


def load_distances(filename):
    """
    Load a precomputed distance matrix from disk.

    This function supports two formats:
    1. Legacy pickle files (`.pkl`, `.pickle`)
    2. Compressed NumPy archives (`.npz`), which are the current default format
       produced by `save_distances`

    Parameters
    ----------
    filename : str
        Path to the saved distance matrix file.

    Returns
    -------
    DistanceMatrix
        A reconstructed `DistanceMatrix` object.

    Notes
    -----
    Pickle loading is retained for backward compatibility. Loading pickled data
    from untrusted sources is unsafe and should be avoided.

    See Also
    --------
    save_distances
    """

    if filename.endswith((".pkl", ".pickle")):
        with open(filename, 'rb') as file:
            distances = pickle.load(file)
        return distances
    else:
        data = np.load(filename, allow_pickle=True)

        matrix = data["matrix"]
        nodes = data["nodes"].tolist()
        names_idx = data["node_to_idx"].tolist()

        node_to_idx = {name: idx for idx, name in enumerate(names_idx)}

        distances = DistanceMatrix.from_components(nodes, node_to_idx, matrix)

        return distances

def _degree_match_null_model(graph):
    degree_dict = {}

    for node in graph.nodes():
        degree = graph.degree(node)

        if degree not in degree_dict:
            degree_dict[degree] = []

        degree_dict[degree].append(node)

    return degree_dict


def _sample_node_proxy(G,S,bucket):
    sampled_nodes = set()

    for node in S:
        available_nodes = bucket[node]

        flag = True
        while(flag):
            chosen_node = random.choice(available_nodes)
            flag = chosen_node in sampled_nodes

        sampled_nodes.add(chosen_node)
    return sampled_nodes


def _sample_preserving_degrees(G,S,bucket):
    sampled_nodes = set()

    for node in S:
        degree = G.degree(node)

        available_nodes = bucket[degree]

        flag = True
        while(flag):
            chosen_node = random.choice(available_nodes)
            flag = chosen_node in sampled_nodes

        sampled_nodes.add(chosen_node)
    return sampled_nodes


def _get_strength_binning(network, bin_size):
    #Calculate the strength of each node
    strengths = {node: sum(weight for _, _, weight in network.edges(node, data='weight')) for node in network.nodes()}

    #Order nodes by their strength in decreasing order
    nodes_sorted_by_strength = sorted(network.nodes(), key=lambda node: strengths[node], reverse=True)

    #Create bins bins for the nodes
    N = len(network.nodes())
    bins = [nodes_sorted_by_strength[i:i+bin_size] for i in range(0, N, bin_size)]

    # 4. Associate each node with its corresponding bin in a dictionary
    node_to_bin = {}
    for bin_list in bins:
        for node in bin_list:
            node_to_bin[node] = bin_list

    return node_to_bin

def _get_degree_binning(net, bin_size=100):
    """Return a histogram of the degrees of the PPI.

    The histogram should have bins with size at least equal to
    `bin_size`. For each bin, the bin bounds l and u should be optimized
    such that a bin with bounds l and u - 1 does is of size smaller
    than `bin_size`. Original code: Rodrigo Dorantes

    Note that for a node with degree d to be in a bin with bounds (l, u],
    we have that l < d <= u.

    Parameters:
    ------------
    net: proximity.Network
        Usually the protein-protein interaction network. **NOTE** If `net` is a
        gt.Graph instance, it should have a `gt.VertexPropertyMap` with the
        node names caled "ids".
    bin_size: int

    Returns
    ---------
    nodes: list
        The nodes of each bin.
    lower: list
        The lower bound of each bin
    upper: list
        The upper bound of each bin.
    """
    #graph = net.Graph
    graph = net
    degree2nodes = {}

    deg = nx.degree(graph)
    for v, d in deg:
        degree2nodes.setdefault(d, list())
        degree2nodes[d].append(v)

    assert type(degree2nodes) is dict, "Not a dict!"

    degrees = sorted(degree2nodes)
    lower = []
    upper = []
    nodes = []
    counter = 0
    cumnodes = []
    l = 0
    for d in degrees:
        counter += len(degree2nodes[d])
        cumnodes.extend(degree2nodes[d])
        if counter < bin_size:
            continue
        u = d
        lower.append(l)
        upper.append(u)
        nodes.append(cumnodes)
        l = u
        cumnodes = []
        counter = 0
    if not nodes:
        raise Exception(f"There should be at least {bin_size} nodes in the graph!")
    if counter:
        upper[-1] = d
        nodes[-1].extend(cumnodes)

    upper[-1] +=1

    return lower, upper, nodes


def _dictionary_from_binning(lower,upper,values):
    degree_dict = {}
    for l, u, nodes in zip(lower, upper, values):
        for degree in range(l, u):
            degree_dict[degree] = nodes
    return degree_dict


def _proximity_dmatrix(net,T,S,D,null_model,node_bucket, n_iter,bin_size):
    if null_model == 'degree_match':
        bucket = _degree_match_null_model(net)
    elif null_model == 'log_binning':
        lower, upper, nodes = _get_degree_binning(net,bin_size=bin_size)
        bucket = _dictionary_from_binning(lower, upper, nodes)
    elif null_model == 'strength_binning':
        bucket = _get_strength_binning(net, bin_size=bin_size)
    elif null_model == 'custom':
        bucket = node_bucket
    elif null_model == 'uniform':
        bucket = set(net.nodes)
    else:
        raise ValueError("Null model should be: 'degree_match'|'log_binning'|'uniform'|'custom'")

    d_c = _get_amspl_dmatrix(T,S,D)
    distribution = []
    for _ in range(n_iter):
        if null_model in ["degree_match","log_binning"]:
            ran_T  = _sample_preserving_degrees(net, T,bucket)
            ran_S  = _sample_preserving_degrees(net, S, bucket)
        elif null_model in ["custom","strength_binning"]:
            ran_T = _sample_node_proxy(net, T, bucket)
            ran_S = _sample_node_proxy(net,S,bucket)
        else:
            ran_T = random.sample(list(bucket), len(T))
            ran_S = random.sample(list(bucket), len(S))
        ran_d_c = _get_amspl_dmatrix(ran_T, ran_S,D)
        distribution.append(ran_d_c)
    mu = np.mean(distribution)
    sigma = np.std(distribution)
    z = (d_c - mu) / sigma

    distribution = np.array(distribution)
    tail = distribution

    tail = tail[tail <= d_c]
    pval_single = len(tail)/len(distribution)

    tail = distribution
    dm = abs(d_c - mu)
    tail = np.abs(tail - mu)

    tail = tail[tail >= dm]
    pval_double = len(tail)/len(distribution)

    return {'d_mu':mu,'d_sigma':sigma,'z_score':z,'p_value_single_tail':pval_single,
            'p_value_double_tail':pval_double, 'raw_amspl':d_c,'dist':distribution}


def _proximity_symmetric(net,T,S,D,null_model, node_bucket,n_iter,bin_size):
    if null_model == 'degree_match':
        bucket = _degree_match_null_model(net)
    elif null_model == 'log_binning':
        lower, upper, nodes = _get_degree_binning(net,bin_size=bin_size)
        bucket = _dictionary_from_binning(lower, upper, nodes)
    elif null_model == 'strength_binning':
        bucket = _get_strength_binning(net, bin_size=bin_size)
    elif null_model == 'custom':
        bucket = node_bucket
    elif null_model == 'uniform':
        bucket = set(net.nodes)
    else:
        raise ValueError("Null model should be: 'degree_match'|'log_binning'|'uniform'|'custom'")

    dts = _get_amspl_dmatrix(T, S,D)
    dst = _get_amspl_dmatrix(S, T,D)
    d = (dts + dst)/2
    distribution = []
    for _ in range(n_iter):
        if null_model in ["degree_match","log_binning"]:
            ran_T  = _sample_preserving_degrees(net, T,bucket)
            ran_S  = _sample_preserving_degrees(net, S, bucket)
        elif null_model in ["custom","strength_binning"]:
            ran_T = _sample_node_proxy(net, T, bucket)
            ran_S = _sample_node_proxy(net,S,bucket)
        else:
            ran_T = random.sample(list(bucket), len(T))
            ran_S = random.sample(list(bucket), len(S))

        rants = _get_amspl_dmatrix(ran_T, ran_S,D)
        ranst = _get_amspl_dmatrix(ran_S, ran_T,D)
        distribution.append((rants+ranst)/2)
    mu = np.mean(distribution)
    sigma = np.std(distribution)
    z = (d - mu) / sigma

    distribution = np.array(distribution)
    tail = distribution

    tail = tail[tail <= d]
    pval_single = len(tail)/len(distribution)

    tail = distribution
    dm = abs(d - mu)
    tail = np.abs(tail - mu)

    tail = tail[tail >= dm]
    pval_double = len(tail)/len(distribution)
    return {'d_mu':mu,'d_sigma':sigma,'z_score':z,'p_value_single_tail':pval_single,
            'p_value_double_tail':pval_double, 'raw_amspl':d,'dist':distribution}



def proximity(net, T, S, D, null_model='degree_match', node_bucket=None,
              n_iter=1000, bin_size=100, symmetric=False):
    """
    Calculate the proximity between two node sets in a network.

    Proximity is defined from the average minimum shortest path length (AMSPL)
    between node sets `T` and `S`, and its statistical significance is assessed
    through repeated random sampling under a chosen null model.

    If `symmetric=False`, the function computes the directed AMSPL from `T` to `S`.
    If `symmetric=True`, it computes the average of AMSPL(T, S) and AMSPL(S, T).

    Parameters
    ----------
    net : networkx.Graph
        Input graph. The graph must be connected.

    T : iterable
        Source node set.

    S : iterable
        Target node set.

    D : DistanceMatrix
        Precomputed distance matrix.

    null_model : {'degree_match', 'log_binning', 'strength_binning',
                  'uniform', 'custom'}, optional
        Null model used to generate randomized node sets.

    node_bucket : dict, optional
        Required when `null_model='custom'`. Maps each original node to a list
        of admissible proxy nodes for resampling.

    n_iter : int, optional
        Number of random iterations. Default is `1000`.

    bin_size : int, optional
        Bin size used by `log_binning` and `strength_binning`. Default is `100`.

    symmetric : bool, optional
        Whether to compute the symmetric version of proximity. Default is `False`.

    Returns
    -------
    dict
        Dictionary with:
        - 'd_mu': mean randomized AMSPL
        - 'd_sigma': standard deviation of randomized AMSPL
        - 'z_score': z-score of the observed AMSPL
        - 'p_value_single_tail': lower-tail empirical p-value
        - 'p_value_double_tail': two-sided empirical p-value
        - 'raw_amspl': observed AMSPL
        - 'dist': array of randomized AMSPL values

    Raises
    ------
    ValueError
        If the network is disconnected, if `n_iter <= 0`, if `bin_size < 1`
        where applicable, if an invalid null model is provided, or if
        `node_bucket` is missing in custom mode.

    Warns
    -----
    UserWarning
        If some nodes in `T` or `S` are not present in the network.
    """

    net_nodes = set(net.nodes)
    sT = set(T)
    sS = set(S)

    #Determine the nodes in A and B that are also in the network
    valid_T = sT & net_nodes
    valid_S = sS & net_nodes


    if nx.number_connected_components(net) > 1:
        raise ValueError("The network is not connected (it contains more than one connected component)")
    if len(sT.difference(valid_T)) > 0:
        warnings.warn("T contains elements that are not present in the network.", UserWarning)
    if len(sS.difference(valid_S)) > 0:
        warnings.warn("S contains elements that are not present in the network.", UserWarning)
    if n_iter <= 0:
        raise ValueError("n_iter must be greater than 0")


    if not null_model in ['degree_match','log_binning','uniform','custom',"strength_binning"]:
        raise ValueError("Null model should be: 'degree_match'|'log_binning'|'uniform'|'custom'")

    if null_model == 'log_binning':
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")

    elif null_model == 'strength_binning':
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
    elif null_model == 'custom':
        if node_bucket == None:
            raise ValueError("In custom mode, node_bucket must be provided")


    if symmetric:
        return _proximity_symmetric(net, valid_T, valid_S, D,null_model,node_bucket,n_iter,bin_size)
    else:
        return _proximity_dmatrix(net, valid_T, valid_S, D,null_model,node_bucket,n_iter,bin_size)


def _sep_dxx(A,D):

    min_distances = np.zeros(len(A))
    idx = 0
    for n in A:
        min_dist = float('inf')
        for m in A:
            if n != m:
                # dis_from_a = D[n]
                # db = dis_from_a[m]
                db = D.get(n,m)
                if(db != None):
                    min_dist = min(min_dist,db)
        min_distances[idx]=min_dist
        idx+=1
    avg = min_distances.mean()
    return avg



def _sep_dab(A,B,D):
    min_distances = np.zeros(len(A))
    idx = 0
    for n in A:
        min_dist = float('inf')
        for m in B:
            # dis_from_a = D[n]
            # db = dis_from_a[m]
            db = D.get(n,m)
            if(db != None):
                min_dist = min(min_dist,db)
        min_distances[idx]=min_dist
        idx+=1
    avg = min_distances.mean()
    return avg



def _sep(A,B,D):
    saa = _sep_dxx(A, D)
    sbb = _sep_dxx(B,D)

    ssab = _sep_dab(A, B, D)*len(A)
    ssba = _sep_dab(B, A, D)*len(B)

    sab = (ssab + ssba)/(len(A) + len(B))

    s = sab + - ((saa + sbb)/2)
    return s



def separation(net,A,B,D):
    """
    Computes the separation between two sets of nodes, A and B, in the network `net` as defined by Menche et al., 2015.
    This measure indicates the relative distance or closeness of two sets of nodes (e.g., genes or proteins) in a network,
    which is useful for understanding relationships like disease-disease interactions or connections between different
    groups of nodes in a biological network.

    The function first checks if the network is connected. If it contains more than one connected component, a ValueError
    is raised. It also verifies the existence of all nodes in sets A and B within the network, issuing warnings for any
    nodes not found.

    Parameters:
    ------------
    net : networkx.Graph
        The input network or graph in which the separation between node sets A and B will be calculated.

    A : container (list, set, etc.)
        A subset of nodes in `net` representing the first group.

    B : container (list, set, etc.)
        A subset of nodes in `net` representing the second group.

    D : DistanceMatrix
        A precomputed distance matrix where D[i][j] provides the shortest distance between nodes i and j.
        This matrix should be generated using the `all_pair_distances` function or an equivalent method.

    Returns:
    ---------
    float
        The separation value between node sets A and B in the network `net`. A smaller value indicates
        that the two sets are closer in the network, while a larger value indicates that they are more
        distantly located.

    Raises:
    --------
    ValueError:
        - If the network is not connected (contains more than one connected component).

    Warnings:
    ----------
    UserWarning:
        - If any elements in `A` or `B` do not exist in the network.

    References:
    ------------
    .. [1] Menche, Jörg, et al.
           "Uncovering disease-disease relationships through the incomplete interactome."
           Science 347.6224 (2015).

    Notes:
    -------
    - Ensure the network is connected to obtain meaningful separation values. Disconnected components may skew the results.
    """

    net_nodes = set(net.nodes)
    sA = set(A)
    sB = set(B)

    #Determine the nodes in A and B that are also in the network
    valid_a = sA & net_nodes
    valid_b = sB & net_nodes

    if nx.number_connected_components(net) > 1:
        raise ValueError("The network is not connected (it contains more than one connected component)")

    if len(sA.difference(valid_a)) > 0:
        warnings.warn("A contains elements that are not present in the network.", UserWarning)
    if len(sB.difference(valid_b)) > 0:
        warnings.warn("B contains elements that are not present in the network.", UserWarning)

    return _sep(valid_a,valid_b,D)



def separation_z_score(net,A,B,D,null_model='degree_match', node_bucket = None, n_iter=1000,bin_size=100):
    """
    Calculates the z-score of the separation between two sets of nodes, A and B, in the network `net`,
    based on randomized node sets with degree-preserving properties.

    The function first checks if the network is connected. If it contains more than one connected component,
    it raises a ValueError. It also checks if all nodes in A and B exist in the network. If any nodes are not
    found, it issues a warning and proceeds with the existing nodes. The function calculates the actual
    separation between node sets A and B, then derives a reference distribution of separations through
    degree-preserving randomizations of the node sets. The z-score indicates how many standard deviations
    the actual separation is from the mean of the reference distribution.

    Parameters:
    ------------
    net : networkx.Graph
        The input network or graph in which the separation between node sets A and B will be assessed.

    A : container (list, set, etc.)
        A subset of nodes in `net` representing the first group.

    B : container (list, set, etc.)
        A subset of nodes in `net` representing the second group.

    D : DistanceMatrix
        A precomputed distance matrix where D[i][j] gives the shortest distance between nodes i and j.
        This matrix should be generated using the `all_pair_distances` function or an equivalent method.

    null_model : str, optional
        Method for degree-preserving randomization. Valid options are 'degree_match', 'log_binning', 'uniform',
        'strength_binning' and 'custom'. Default is 'degree_match'.

    node_bucket : dictionary, optional
        A collection of nodes to be used in 'custom' mode, mandatory when the null_model is set to 'custom'.
        This parameter should be a dictionary where each key represents a node ('node_k') from the network,
        and the corresponding value is a list of alternative nodes ('proxy_i').
        These alternatives are used by the null model for resampling:
        node_bucket[node_k] = [proxy_1, proxy_2, ..., proxy_m].
        Here, each 'proxy_i' serves as a potential substitute to be sampled in place of 'node_k'.


    n_iter : int, optional
        Number of random sampling iterations used to derive the reference distribution of separations.
        Default is 1000.

    bin_size : int, optional
        Bin size used when `null_model` is `'log_binning'` or
        `'strength_binning'`. Default is `100`.

    Returns:
    ---------
    dict
        A dictionary containing:
        - 'd_mu': Mean separation from the randomized samples.
        - 'd_sigma': Standard deviation of separations from the randomized samples.
        - 'z_score': Z-score of the actual separation against the randomized samples.
        - 'p_value_single_tail': One-tail P-value associated with the proximity z-score
        - 'p_value_double_tail': Two-tail P-value associated with the proximity z-score
        - 'raw_separation': Actual separation value between node sets A and B.
        - 'dist': List of separations from each randomization iteration.

    Raises:
    --------
    ValueError:
        - If the network is not connected (contains more than one connected component).
        - If 'n_iter' is less than or equal to 0.
        - If 'bin_size' is less than 1 when 'log_binning' or 'strength_binning' is used.
        - If 'null_model' is not one of ['degree_match', 'log_binning', 'strength_binning', 'uniform', 'custom'].
        - If 'node_bucket' is not provided when 'null_model' is 'custom'.

    Warnings:
    ----------
    UserWarning:
        - If any elements in `A` or `B` do not exist in the network.

    Notes:
    -------
    - The degree-preserving randomization ensures that the randomized node samples have a degree distribution similar
      to the original sets, ensuring a fair comparison.
    - The network should be connected to obtain meaningful separation values. Disconnected components may skew the results.
    """
    net_nodes = set(net.nodes)
    sA = set(A)
    sB = set(B)

    #Determine the nodes in A and B that are also in the network
    valid_a = sA & net_nodes
    valid_b = sB & net_nodes

    if nx.number_connected_components(net) > 1:
        raise ValueError("The network is not connected (it contains more than one connected component)")

    if len(sA.difference(valid_a)) > 0:
        warnings.warn("A contains elements that are not present in the network.", UserWarning)
    if len(sB.difference(valid_b)) > 0:
        warnings.warn("B contains elements that are not present in the network.", UserWarning)
    if n_iter <= 0:
        raise ValueError("n_iter must be greater than 0")


    if null_model == 'degree_match':
        bucket = _degree_match_null_model(net)
    elif null_model == 'log_binning':
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
        lower, upper, nodes = _get_degree_binning(net,bin_size=bin_size)
        bucket = _dictionary_from_binning(lower, upper, nodes)
    elif null_model == 'strength_binning':
        bucket = _get_strength_binning(net, bin_size=bin_size)
    elif null_model == 'uniform':
        bucket = set(net.nodes)
    elif null_model == 'custom':
        bucket = node_bucket
        if node_bucket == None:
            raise ValueError("In custom mode, node_bucket must be provided")
    else:
        raise ValueError("Null model should be: 'degree_match'|'log_binning'|'uniform'|'custom'")

    s = _sep(valid_a, valid_b, D)
    distribution = []
    for _ in range(n_iter):
        if null_model in ["degree_match","log_binning"]:
            ran_A  = _sample_preserving_degrees(net, valid_a,bucket)
            ran_B  = _sample_preserving_degrees(net, valid_b, bucket)
        elif null_model in ["custom","strength_binning"]:
            ran_A = _sample_node_proxy(net, valid_a, bucket)
            ran_B = _sample_node_proxy(net,valid_b,bucket)
        else:
            ran_A = random.sample(list(bucket), len(valid_a))
            ran_B = random.sample(list(bucket), len(valid_b))
        ran_sep = _sep(ran_A, ran_B,D)
        distribution.append(ran_sep)
    mu = np.mean(distribution)
    sigma = np.std(distribution)
    z = (s - mu) / sigma

    distribution = np.array(distribution)
    tail = distribution

    tail = tail[tail <= s]
    pval_single = len(tail)/len(distribution)

    tail = distribution
    dm = abs(s - mu)
    tail = np.abs(tail - mu)

    tail = tail[tail >= dm]
    pval_double = len(tail)/len(distribution)

    return {'d_mu':mu,'d_sigma':sigma,'z_score':z,'p_value_single_tail':pval_single,
            'p_value_double_tail':pval_double,'raw_separation':s,'dist':distribution}


def extract_lcc(A, net):
    """
    Extracts the Largest Connected Component (LCC) from a subgraph of `net` induced by the set of nodes `A`.

    This function first verifies that all nodes in `A` exist in `net`. If any nodes in `A` do not exist in `net`,
    it issues a warning and proceeds with the existing nodes. It then creates a subgraph of `net` that includes
    only the nodes in `A`. The largest connected component (LCC) within this subgraph is identified and returned
    as a subgraph of `net`.

    Parameters:
    -------------
    A : iterable
        An iterable (e.g., list, set) of nodes. The subgraph is induced by these nodes.

    net : networkx.Graph
        The original graph from which the subgraph and LCC are derived.

    Returns:
    ----------
    networkx.Graph
        A subgraph of `net` containing only the nodes in the largest connected component of the subgraph
        induced by `A`.

    Notes:
    --------
    - If `A` is empty or consists of nodes not in `net`, the function will issue a warning and return an empty graph.
    - The function utilizes NetworkX's methods for subgraph creation and connected component identification.

    Example:
    ----------
    >>> G = nx.Graph([(1, 2), (2, 3), (3, 4), (5, 6)])
    >>> A = [1, 2, 3, 7]
    >>> G_sub = extract_lcc(A, G)
    >>> list(G_sub.nodes)
    [1, 2, 3]
    """

    # Check for nodes in A that are not in net
    net_nodes = set(net.nodes())
    A_nodes = set(A)
    missing_nodes = A_nodes - net_nodes

    if missing_nodes:
        warnings.warn(f"The following nodes are not in the network and will be ignored: {missing_nodes}", UserWarning)
        A_nodes = A_nodes - missing_nodes

    if not A_nodes:
        return nx.Graph()  # Return an empty graph if no valid nodes are left

    G = net.subgraph(A_nodes)
    largest_cc = max(nx.connected_components(G), key=len)
    G_sub = G.subgraph(largest_cc)

    return G_sub


def _strength_binning_null_model(net, bin_size=100, weight='weight'):
    if bin_size < 1:
        raise ValueError("bin_size must be greater or equal than 1")

    strength = dict(net.degree(weight=weight))
    ordered_nodes = sorted(strength, key=lambda n: strength[n])

    bucket = {}
    for start in range(0, len(ordered_nodes), bin_size):
        group = ordered_nodes[start:start + bin_size]
        for node in group:
            bucket[node] = group

    return bucket


def lcc_significance(net, A, null_model='degree_match', node_bucket=None,
                     n_iter=1000, bin_size=100, strength_weight='weight'):
    """
    Calculate the statistical significance of the size of the Largest Connected
    Component (LCC) induced by a node set `A`.

    The function compares the observed LCC size against a null distribution
    obtained by repeatedly resampling node sets of the same size under a chosen
    null model.

    Parameters
    ----------
    net : networkx.Graph
        Input network.

    A : iterable
        Node set whose induced LCC is to be evaluated.

    null_model : {'degree_match', 'log_binning', 'strength_binning', 'uniform', 'custom'}, optional
        Null model used for resampling. Default is `'degree_match'`.

    node_bucket : dict, optional
        Required when `null_model='custom'`. Maps each original node to a list
        of admissible proxy nodes.

    n_iter : int, optional
        Number of random iterations. Default is `1000`.

    bin_size : int, optional
        Bin size used when `null_model='log_binning'` or
        `null_model='strength_binning'`. Default is `100`.

    strength_weight : str, optional
        Edge attribute used as weight when `null_model='strength_binning'`.
        Default is `'weight'`.

    Returns
    -------
    dict
        Dictionary with:
        - 'd_mu': mean LCC size under the null model
        - 'd_sigma': standard deviation of null LCC sizes
        - 'z_score': z-score of the observed LCC size
        - 'p_val': empirical upper-tail p-value
        - 'lcc': largest connected component subgraph
        - 'lcc_size': size of the observed LCC
        - 'dist': array of null-model LCC sizes

    Raises
    ------
    ValueError
        If `n_iter <= 0`, if `bin_size < 1` where applicable, if an invalid
        null model is provided, or if `node_bucket` is missing in custom mode.

    Warns
    -----
    UserWarning
        If some nodes in `A` are not present in the network.

    Notes
    -----
    - In `'uniform'` mode, nodes are sampled uniformly from the network,
      excluding the observed nodes in `A`.
    - In `'degree_match'`, `'log_binning'`, `'strength_binning'`, and
      `'custom'` modes, the function delegates sampling to
      `_sample_preserving_degrees`.
    - If the null distribution has zero variance, the z-score is returned as
      `np.nan`.
    """

    sA = set(A)
    set_a = sA & set(net.nodes())

    if len(sA.difference(set_a)) > 0:
        warnings.warn(
            "A contains elements that are not present in the network.",
            UserWarning
        )

    if n_iter <= 0:
        raise ValueError("n_iter must be greater than 0")

    if null_model == 'degree_match':
        bucket = _degree_match_null_model(net)

    elif null_model == 'log_binning':
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
        lower, upper, nodes = _get_degree_binning(net, bin_size=bin_size)
        bucket = _dictionary_from_binning(lower, upper, nodes)

    elif null_model == 'strength_binning':
        bucket = _strength_binning_null_model(
            net,
            bin_size=bin_size,
            weight=strength_weight
        )

    elif null_model == 'uniform':
        bucket = set(net.nodes).difference(set_a)

    elif null_model == 'custom':
        if node_bucket is None:
            raise ValueError("In custom mode, node_bucket must be provided")
        bucket = node_bucket

    else:
        raise ValueError(
            "Null model should be in "
            "['degree_match'|'log_binning'|'strength_binning'|'uniform'|'custom']"
        )

    lcc = extract_lcc(set_a, net)

    distribution = []
    for i in range(n_iter):
        if null_model != 'uniform':
            rs = _sample_preserving_degrees(net, set_a, bucket)
        else:
            rs = random.sample(list(bucket), len(set_a))

        sub = extract_lcc(rs, net)
        distribution.append(len(sub))

        if i % 100 == 0:
            print(f"\rIter {i} of {n_iter}", end="")
    print("")

    l_lcc = len(lcc)
    distribution = np.array(distribution)
    mu = np.mean(distribution)
    sigma = np.std(distribution)
    z = np.nan if sigma == 0 else (l_lcc - mu) / sigma

    S = distribution[distribution >= l_lcc]
    pval = len(S) / len(distribution)

    return {
        'd_mu': mu,
        'd_sigma': sigma,
        'z_score': z,
        'p_val': pval,
        'lcc': lcc,
        'lcc_size': l_lcc,
        'dist': distribution
    }


def _check_dictionary_integrity(network,dictionary):
    newDict = {}
    nodes = set(network.nodes())

    for key,val in dictionary.items():
        sA = set(val)
        set_a = sA & set(nodes)

        if len(sA.difference(set_a)) > 0:
            warnings.warn("WARNING!! {key} contains elements that are not present in the network. Keeping valid elements.", UserWarning)

            if(len(set_a) == 0):
                warnings.warn("WARNING!! {key} does not contain valid nodes and will be avoided.", UserWarning)
        newDict[key] = set_a
    return newDict


def _to_dictionary(results, properties):
    res = {}
    for p in properties:
        res[p] = results[p]

    return res


def _to_dict_tables(results, properties,target_names,source_names):
    table_res = {}
    for p in properties:
        df = pd.DataFrame(columns=target_names, index=source_names)

        for r in results:
            source = r[0]
            target = r[1]
            score = r[2]

            df.loc[source,target] = score[p]
        df = df.apply(pd.to_numeric, errors='coerce')
        table_res[p] = df

    return table_res


def to_dictionary(dataframe, group_names, node_names):
    """
    Converts a DataFrame into a dictionary where the keys are unique group names
    and the values are sets of node names associated with each group.

    Parameters:
    -------------
    dataframe : pandas.DataFrame
        The input DataFrame containing the data to be converted into a dictionary.

    group_names : str
        The column name in the DataFrame that contains the group identifiers.

    node_names : str
        The column name in the DataFrame that contains the node identifiers.

    Returns:
    ----------
    dict
        A dictionary where each key is a unique group name from the `group_names`
        column, and each value is a set of node names from the `node_names` column
        associated with that group.

    Example:
    ----------
    >>> import pandas as pd
    >>> data = {'group': ['A', 'A', 'B', 'B', 'C'],
                'node': ['x', 'y', 'x', 'z', 'y']}
    >>> df = pd.DataFrame(data)
    >>> to_dictionary(df, 'group', 'node')
    {'A': {'x', 'y'}, 'B': {'x', 'z'}, 'C': {'y'}}
    """
    unique_names = list(dataframe[group_names].unique())

    res = {}

    for l in unique_names:
        nodes = dataframe[dataframe[group_names] == l]
        nodes = set(nodes[node_names])

        res[l] = nodes
    return res

def _init_screening_worker(
    shm_name,
    shape,
    dtype_str,
    dm_nodes,
    dm_node_to_idx,
    network,
    sources,
    targets,
    node_bucket,
):
    global _SCREEN_NETWORK, _SCREEN_SOURCES, _SCREEN_TARGETS
    global _SCREEN_NODE_BUCKET, _SCREEN_DISTANCE_MATRIX
    global _SCREEN_SHM, _SCREEN_MATRIX

    _SCREEN_SHM = SharedMemory(name=shm_name)
    dtype = np.dtype(dtype_str)
    _SCREEN_MATRIX = np.ndarray(shape, dtype=dtype, buffer=_SCREEN_SHM.buf)

    D = DistanceMatrix()
    D.nodes = dm_nodes
    D.node_to_idx = dm_node_to_idx
    D.matrix = _SCREEN_MATRIX

    _SCREEN_DISTANCE_MATRIX = D
    _SCREEN_NETWORK = network
    _SCREEN_SOURCES = sources
    _SCREEN_TARGETS = targets
    _SCREEN_NODE_BUCKET = node_bucket

def _calculate_score_worker(task):
    source, target, score, properties, null_model, n_iter, bin_size, symmetric = task

    network = _SCREEN_NETWORK
    sources = _SCREEN_SOURCES
    targets = _SCREEN_TARGETS
    node_bucket = _SCREEN_NODE_BUCKET
    distance_matrix = _SCREEN_DISTANCE_MATRIX

    source_nodes = sources[source]
    target_nodes = targets[target]

    if score == "proximity":
        scores = proximity(
            network,
            source_nodes,
            target_nodes,
            distance_matrix,
            null_model=null_model,
            node_bucket=node_bucket,
            n_iter=n_iter,
            bin_size=bin_size,
            symmetric=symmetric,
        )
        results = _to_dictionary(scores, properties)

    elif score == "separation_z_score":
        scores = separation_z_score(
            network,
            source_nodes,
            target_nodes,
            distance_matrix,
            null_model=null_model,
            node_bucket=node_bucket,
            n_iter=n_iter,
            bin_size=bin_size,
        )
        results = _to_dictionary(scores, properties)

    elif score == "separation":
        scores = separation(network, source_nodes, target_nodes, distance_matrix)
        results = {"raw_separation": scores}

    else:
        raise ValueError("Invalid score.")

    print(f"{source}-{target} finished")
    return (source, target, results)

def _distance_matrix_to_shared(distance_matrix):
    mat = distance_matrix.matrix
    shape = mat.shape
    dtype = mat.dtype

    shm = SharedMemory(create=True, size=mat.nbytes)
    shm_array = np.ndarray(shape, dtype=dtype, buffer=shm.buf)
    shm_array[:] = mat[:]

    return shm, shape, dtype

def screening(
    sources,
    targets,
    network,
    distance_matrix,
    score="proximity",
    properties=["z_score"],
    null_model="degree_match",
    node_bucket=None,
    n_iter=1000,
    bin_size=100,
    symmetric=False,
    n_procs=None,
    chunksize=1,
):
    """
    Screen all source-target group pairs in a network using proximity,
    separation z-score, or raw separation.

    Each entry in `sources` and `targets` defines a named node set. The function
    evaluates every pair `(source_group, target_group)` and returns the requested
    score properties as tables.

    Parameters
    ----------
    sources : dict
        Dictionary mapping source group names to node sets.

    targets : dict
        Dictionary mapping target group names to node sets.

    network : networkx.Graph
        Input graph. The graph must be connected.

    distance_matrix : DistanceMatrix
        Precomputed distance matrix used by the scoring functions.

    score : {'proximity', 'separation_z_score', 'separation'}, optional
        Score to compute for each source-target pair. Default is `'proximity'`.

    properties : list of str, optional
        Statistics to extract from the selected score:
        - for `'proximity'`: `'z_score'`, `'p_value_single_tail'`,
          `'p_value_double_tail'`, `'raw_amspl'`
        - for `'separation_z_score'`: `'z_score'`, `'p_value_single_tail'`,
          `'p_value_double_tail'`, `'raw_separation'`
        - for `'separation'`: `'raw_separation'`

    null_model : {'degree_match', 'log_binning', 'strength_binning',
                  'uniform', 'custom'}, optional
        Null model used for resampling when applicable.

    node_bucket : dict, optional
        Required when `null_model='custom'`.

    n_iter : int, optional
        Number of random iterations for significance calculations. Default is `1000`.

    bin_size : int, optional
        Bin size used by `log_binning` and `strength_binning`. Default is `100`.

    symmetric : bool, optional
        Passed to `proximity` when `score='proximity'`. Default is `False`.

    n_procs : int, optional
        Number of worker processes. If `None`, uses all available CPUs.

    chunksize : int, optional
        Chunk size passed to `ProcessPoolExecutor.map`. Default is `1`.

    Returns
    -------
    dict
        Dictionary mapping each requested property to a pandas DataFrame whose
        rows are source group names and columns are target group names.

    Raises
    ------
    ValueError
        If the network is disconnected, if parameters are invalid, or if
        `node_bucket` is missing in custom mode.

    Notes
    -----
    The function uses multiprocessing with shared memory so that the full
    distance matrix is not duplicated across worker processes.
    """

    if nx.number_connected_components(network) > 1:
        raise ValueError("The network is not connected (it contains more than one connected component)")

    if n_iter <= 0:
        raise ValueError("n_iter must be greater than 0")

    if null_model == "log_binning":
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
    elif null_model == "strength_binning":
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
    elif null_model == "custom":
        if node_bucket is None:
            raise ValueError("In custom mode, node_bucket must be provided")
    else:
        if null_model not in ["degree_match", "uniform"]:
            raise ValueError("Null model should be: 'degree_match'|'log_binning'|'uniform'|'custom'")

    if score not in ["proximity", "separation_z_score", "separation"]:
        raise ValueError("Score should be: 'proximity','separation_z_score','separation'")

    if score == "proximity":
        if not set(properties) <= {"z_score", "p_value_single_tail", "p_value_double_tail", "raw_amspl"}:
            raise ValueError("prop attribute should be: 'z_score','p_value_single_tail','p_value_double_tail','raw_amspl'")
    elif score == "separation_z_score":
        if not set(properties) <= {"z_score", "p_value_single_tail", "p_value_double_tail", "raw_separation"}:
            raise ValueError("prop attribute should be: 'z_score','p_value_single_tail','p_value_double_tail','raw_separation'")
    elif score == "separation":
        if not set(properties) <= {"raw_separation"}:
            raise ValueError("prop attribute should be: 'raw_separation'")

    valid_sources = _check_dictionary_integrity(network, sources)
    valid_targets = _check_dictionary_integrity(network, targets)

    source_names = list(valid_sources.keys())
    target_names = list(valid_targets.keys())

    num_cpus = os.cpu_count() if n_procs is None else n_procs

    shm, shape, dtype = _distance_matrix_to_shared(distance_matrix)

    tasks = [
        (source, target, score, properties, null_model, n_iter, bin_size, symmetric)
        for source in source_names
        for target in target_names
    ]

    ctx = get_context("spawn")

    try:
        with ProcessPoolExecutor(
            max_workers=num_cpus,
            mp_context=ctx,
            initializer=_init_screening_worker,
            initargs=(
                shm.name,
                shape,
                dtype.str,
                distance_matrix.nodes,
                distance_matrix.node_to_idx,
                network,
                valid_sources,
                valid_targets,
                node_bucket,
            ),
        ) as ex:
            results = list(ex.map(_calculate_score_worker, tasks, chunksize=chunksize))

    finally:
        shm.close()
        shm.unlink()

    dict_tables = _to_dict_tables(results, properties, target_names, source_names)
    return dict_tables


def random_walk(G, seed, restart_prob=0.15):
    """
    Perform a random walk with restart (RWR) on a graph using NetworkX PageRank
    with a personalized restart distribution.

    Seed nodes define the personalization vector:
    - if `seed` is a list or set, all valid seed nodes receive equal weight
    - if `seed` is a dictionary, values are interpreted as seed weights and are
      normalized internally

    The function returns both the raw random-walk score and a degree-normalized
    score (`BScore`) that downweights high-degree hubs.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.

    seed : list, set, or dict
        Seed nodes or weighted seed nodes.

    restart_prob : float, optional
        Probability of restarting at a seed node at each step. Default is `0.15`.

    Returns
    -------
    pandas.DataFrame
        DataFrame sorted by descending `BScore`, with columns:
        - 'Node'
        - 'Degree'
        - 'Score'
        - 'BScore'

    Raises
    ------
    ValueError
        If no valid seed nodes are found in the graph.

    TypeError
        If `seed` is not a list, set, or dict.

    Example
    -------
    >>> import networkx as nx
    >>> G = nx.karate_club_graph()
    >>> result = random_walk(G, seed=[0, 33], restart_prob=0.2)
    >>> print(result.head())
    """

    # Convert restart probability to damping factor
    alpha = 1 - restart_prob

    # Create personalization vector
    if isinstance(seed, (list, set)):
        filtered_seed = [s for s in seed if s in G]
        if not filtered_seed:
            raise ValueError("No valid seed nodes found in the graph.")
        weight = 1.0 / len(filtered_seed)
        personalization = {n: weight for n in filtered_seed}

    elif isinstance(seed, dict):
        filtered_seed = {k: v for k, v in seed.items() if k in G}
        if not filtered_seed:
            raise ValueError("No valid seed nodes found in the graph.")
        total = sum(filtered_seed.values())
        personalization = {k: v / total for k, v in filtered_seed.items()}

    else:
        raise TypeError("Seed must be a list, set, or dict.")
    
    # Perform a random walk
    pr = nx.pagerank(G, alpha=alpha, personalization=personalization)

    # Format results
    df = pd.DataFrame({
        'Node': list(pr.keys()),
        'Degree': [G.degree(n) for n in pr.keys()],
        'Score': list(pr.values()), 
    })

    df['BScore'] = df['Score'] / df['Degree']
    df = df.sort_values(by='BScore', ascending=False).reset_index(drop=True)
    return df