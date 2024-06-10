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
- `save_distances`: Saves the precomputed distance matrix to a `pickle` file.
- `load_distances`: Loads a precomputed distance matrix from a `pickle` file.
- `get_amspl`: Calculates the Average Minimum Shortest Path Length between nodes.
- `proximity`: Computes the proximity between two sets of nodes in a graph.
- `separation`: Calculates the separation between two sets of nodes in a network.
- `separation_z_score`: Determines the z-score of the separation between two node sets based on randomized samples.
- `screening`: Screens for proximity/separation between sets of source and target nodes.


Distance metrics.
------------------

When calculating the distance matrix, four information flow metrics are available to the user:

- `shortest_path`: Distance between nodes is based in the length of the path with the least number of edges or lowest total weight that connects two nodes
- `random_walk`: Distance between nodes is based in the probability of reaching one node from another via a random walk.
- `biased_random_walk`: Same as random_walk but compensating the bias induced by the degree of the target node.
- `communicability`: Distance between nodes is based on the concept of communicability, defined as the ability of nodes to communicate or send information through all available paths in a network, considering the indirect and direct connections.


Null models.
-------------

The NetMedPy functions involving statistical analysis allow the user to select among the following null models:

- `degree_match`: selects random samples replicating the original node-set's degree distribution.
- `log_binning`: categorizes the degrees of all nodes within the network into logarithmically sized bins. Samples are then drawn by matching the degree of the original nodes to those within the corresponding bins.
- `strength_binning`: analogous to log_binning, using the strength of the nodes instead of their degrees.
- `uniform`: randomly selects nodes from the entire network, disregarding their degree or strength.
- `custom`: allows users to specify custom null models.


These functions use both exact and approximate methods for degree-preserving and non-degree preserving randomization of node sets. Additionally,
precomputed distance matrices are leveraged for efficient computation.


Required packages:
-------------------

    - networkx
    - numpy
    - pickle
    - multiprocessing
    - random
    - scipy
    - ray


Authors:
---------
    - Andres Aldana Gonzalez (a.aldana@northeastern.edu)
    - Rodrigo Dorantes Gilardi (r.dorantesgilardi@northeastern.edu)


References:
------------
    - Menche, Jörg, et al. "Uncovering disease-disease relationships through the incomplete interactome." Science 347.6224 (2015). DOI 10.1126/science.1257601
    - Guney, Emre, et al.  "Network-based in silico drug efficacy screening." Nature Communications 7,1 (2015). DOI 10.1038/ncomms10331
    - Estrada, Ernesto, and Naomichi Hatano. "Communicability in complex networks." Physical Review E 77.3 (2008): 036111.
    - Masuda, Naoki, Mason A. Porter, and Renaud Lambiotte. "Random walks and diffusion on networks." Physics reports 716 (2017): 1-58.
    - Le, Duc-Hau. "Random walk with restart: A powerful network propagation algorithm in Bioinformatics field." 2017 4th NAFOSTED Conference on Information and Computer Science. IEEE, 2017.
"""

import networkx as nx
import numpy as np
import pickle
from multiprocessing import cpu_count
import random
import ray
import warnings
import os
import pandas as pd
try:
    # This works when the package is installed via pip
    from .DistanceMatrix import DistanceMatrix
except ImportError:
    # This works when using PYTHONPATH
    from DistanceMatrix import DistanceMatrix


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
            if(db != None):
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
            if(db != None):
                min_dist = min(min_dist,db)
        min_distances[idx]=min_dist
        idx+=1
    avg = min_distances.mean()
    return avg


def get_amspl(net,A,B,D):
    """
    Returns the average minimum distance between each node in A and all nodes in B, using the
    distance matrix D to access precomputed distances between nodes.

    Parameters:
    ------------
    net : networkx.Graph
        The input network/graph for which distances need to be computed.

    A : Iterable (list, set, etc.)
        A collection of nodes from which the shortest paths to nodes in B will be computed.

    B : Iterable (list, set, etc.)
        A collection of nodes to which the shortest paths from nodes in A will be computed.

    D : dict of dicts
        A distance matrix where D[i][j] gives the precomputed shortest distance between nodes i and j.
        If there's no path between i and j, D[i][j] should be None.

    Returns:
    ---------
    avg : float
        The average of the minimum distances between each node in A and all nodes in B.

    """
    net_nodes = set(net.nodes)
    sA = set(A)
    sB = set(B)

    #Determine the nodes in A and B that are also in the network
    valid_a = sA & net_nodes
    valid_b = sB & net_nodes

    return _get_amspl_dmatrix(valid_a,valid_b,D)



@ray.remote
def _single_shortest_path(source_nodes, graph, node_to_idx):

    mat_array = np.full((len(source_nodes),len(graph) + 1), np.inf)
    current_row = 0

    for s in source_nodes:
        if nx.get_edge_attributes(graph,'weight'):
            d = nx.shortest_path_length(graph,s,weight='weight')
        else:
            d = nx.shortest_path_length(graph,s)

        source_index = node_to_idx[s]

        mat_array[current_row,0] = source_index

        for target,distance in d.items():
            idx = node_to_idx[target]
            mat_array[current_row,idx+1] = distance
        del d
        current_row += 1

    print(f"Process {os.getpid()} finished")

    return mat_array



def _spl_distance(graph, node_to_idx,num_cpus,n_tasks):
    ray.shutdown()
    ray.init(num_cpus=num_cpus)

    graph_ref = ray.put(graph)
    node_to_idx_ref = ray.put(node_to_idx)

    chunks = _split_into_chunks(list(graph.nodes), n_tasks)

    results = [_single_shortest_path.remote(chunk,graph_ref,node_to_idx_ref)
            for chunk in chunks]

    results = ray.get(results)

    ray.shutdown()
    return results



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


@ray.remote
def _custom_source_distance(source_nodes, graph, node_to_idx, distance, kwargs):
     mat_array = np.full((len(source_nodes),len(graph) + 1), np.inf)
     current_row = 0

     for s in source_nodes:
         res_dict = distance(s,graph,**kwargs)
         source_index = node_to_idx[s]

         mat_array[current_row,0] = source_index

         for target,d in res_dict.items():
             idx = node_to_idx[target]
             mat_array[current_row,idx+1] = d
         del res_dict
         current_row += 1

     print(f"Process {os.getpid()} finished")

     return mat_array



def _custom_all_distance(graph, node_to_idx, distance, num_cpus, n_tasks, kwargs):
    ray.shutdown()
    ray.init(num_cpus=num_cpus)

    graph_ref = ray.put(graph)
    distance_ref = ray.put(distance)
    node_to_idx_ref = ray.put(node_to_idx)

    chunks = _split_into_chunks(list(graph.nodes), n_tasks)

    res = [_custom_source_distance.remote(chunk, graph_ref, node_to_idx_ref, distance_ref, kwargs)
           for chunk in chunks]

    results = ray.get(res)

    ray.shutdown()
    return results


def all_pair_distances(graph,distance="shortest_path",custom_distance=None,
                       reset=0.2, n_processors=None, n_tasks=None, **kwargs):
    """
    Calculates distances between every pair of nodes in a graph according to the specified method and returns
    a DistanceMatrix object. This function supports multiple distance calculation methods, including shortest
    path, various types of random walks, and user-defined methods.

    Parameters:
    ------------
    graph : networkx.Graph
        The input graph for which pairwise distances will be computed. Nodes should be unique and hashable.

    distance : str, optional
        The method used to calculate distances. Options include 'shortest_path', 'random_walk', 'biased_random_walk',
        'communicability', and 'custom'. Default is 'shortest_path'.

    custom_distance : function, optional
        A custom function for distance calculation, used when 'distance' is set to 'custom'. This function should
        have the signature `function_name(a, networkx.Graph, **kwargs)`, where 'a' is the source node, 'networkx.Graph'
        is the graph, and '**kwargs' are additional arguments. It should return a dictionary where each key 'k' is a
        target node and the value is the distance from the source node 'a' to 'k'. The dictionary must include distances
        from 'a' to all other nodes in the graph.

    reset : float, optional
        The reset probability for random walk-based distance calculations. Must be between 0 and 1. Default is 0.2.

    kwargs : dict
        Additional keyword arguments passed to the user-defined distance function.

    Returns:
    ---------
    DistanceMatrix
        A DistanceMatrix object where the value at `matrix[node1][node2]` gives the calculated distance
        from `node1` to `node2` according to the specified method.

    Raises:
    --------
    ValueError
        - If the network is not connected (contains more than one connected component).
        - If the 'distance' parameter is not one of the valid options.
        - If 'reset' is not within the range [0, 1] for random walk-based distances.

    Notes:
    -------
    - The function utilizes multiple processes and should be invoked from the main execution environment only.
    - For user-defined distances, ensure the custom function is correctly structured and returns the expected format.
    """

    #Parameter verification
    if nx.number_connected_components(graph) > 1:
        raise ValueError("The network is not connected (it contains more than one connected component)")

    if not distance in ["shortest_path","random_walk","biased_random_walk","communicability","custom"]:
        raise ValueError("distance must be shortest_path|random_walk|biased_random_walk|communicability|custom")

    if distance=="random_walk" or distance=="biased_random_walk":
        if reset < 0 or reset > 1:
            raise ValueError("Reset for random walks must comply 0 <= reset <= 1")

    if n_processors == None:
        num_cpus = cpu_count()
    else:
        num_cpus = n_processors

    if n_tasks == None:
        n_tasks = num_cpus

    if n_tasks < num_cpus:
        raise ValueError("Number of tasks should be larger or equal the number of processors.")


    D = DistanceMatrix()
    D._from_name_list(list(graph.nodes()))
    del D.matrix
    D.matrix = None

    if distance == "shortest_path":
        spl_res = _spl_distance(graph, D.node_to_idx,num_cpus,n_tasks)
        res_mat = np.zeros((len(graph),len(graph)))

        for sub_mat in spl_res:
            for i in range(len(sub_mat)):
                idx = int(sub_mat[i,0])
                content = sub_mat[i,1:]
                res_mat[idx] = content
        del spl_res
        D.matrix = res_mat

    elif distance == "random_walk":
        M,node_to_index,index_to_node = _RWR_Matrix(graph, reset)

        #  Normalize values
        M = np.log(M)
        ma = np.max(M)
        mi = np.min(M)
        M = 1 - ((M - mi)/(ma - mi))

        D.nodes = graph.nodes()
        D.node_to_idx = node_to_index
        D.matrix = M

    elif distance == "biased_random_walk":
        M,node_to_index,index_to_node = _RWR_Matrix(graph, reset)

        degs = []
        for i in range(len(graph)):
            degs.append(graph.degree(index_to_node[i]))

        degs = np.array(degs)

        #  Remove bias dividing by node degree
        B = M / degs
        B = B.T
        B = B / B.sum(axis=0)
        B = B.T

        #  Normalize values
        B = np.log(B)
        ma = np.max(B)
        mi = np.min(B)
        B = 1 - ((B - mi)/(ma - mi))

        D.nodes = graph.nodes()
        D.node_to_idx = node_to_index
        D.matrix = B

    elif distance == "communicability":
        D = DistanceMatrix()
        D._from_name_list(list(graph.nodes()))
        comm = nx.communicability_exp(graph)

        for a,d in comm.items():
            for b,c in d.items():
                D.put(a,b,c)
        del comm

        B = np.log(D.matrix)
        ma = np.max(B)
        mi = np.min(B)
        B = 1 - ((B - mi)/(ma - mi))

        D.matrix = B

    elif distance == "custom":
        res = _custom_all_distance(graph, D.node_to_idx, custom_distance ,num_cpus,n_tasks, kwargs)
        res_mat = np.zeros((len(graph),len(graph)))

        for sub_mat in res:
            for i in range(len(sub_mat)):
                idx = int(sub_mat[i,0])
                content = sub_mat[i,1:]

                res_mat[idx] = content
        del res
        D.matrix = res_mat
    else:
        pass #This is actually not possible due to the initial verification
    return D



def save_distances(distances, filename):
    """
    Saves the precomputed distance matrix to a file using the `pickle` module.

    This function serializes the given distance matrix and writes it to a specified file. The
    saved file can later be loaded to quickly retrieve the distance matrix without needing to
    recompute the distances.

    Parameters:
    ------------
    distances : DistanceMatrix
        The distance matrix D, where D[a][b] represents the shortest path distance from
        source node a to target node b.

    filename : str
        The path and name of the file to which the distances should be saved. If the file
        already exists, it will be overwritten.

    Notes:
    ------------
    The saved file is in binary format due to the usage of the `pickle` module. Always be cautious
    when loading pickled data from untrusted sources as it can be a security risk.
    """
    with open(filename, 'wb') as file:
        pickle.dump(distances, file)


def load_distances(filename):
    """
     Loads a precomputed distance matrix from a file using the `pickle` module.

     This function deserializes and retrieves a distance matrix (of the DistanceMatrix class)
     that was previously saved to a file. This operation is the inverse of saving the matrix using
     the `pickle` module.

     Parameters:
     ------------
     filename : str
         The path and name of the file from which the distance matrix is to be loaded. The file
         should have been previously saved using the `pickle` module, typically via the `save_distances` function.

     Returns:
     ----------
     DistanceMatrix
         The distance matrix D, where D[a][b] represents the shortest path distance from
         source node a to target node b.

     Notes:
     -------
     - The loaded file is in binary format due to the usage of the `pickle` module. Exercise
       caution when loading pickled data from untrusted sources, as it can pose a security risk.
     """
    with open(filename, 'rb') as file:
        distances = pickle.load(file)
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
    tail = np.abs(tail)
    tail = tail[tail >= d_c]
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
    tail = np.abs(tail)
    tail = tail[tail >= d]
    pval_double = len(tail)/len(distribution)
    return {'d_mu':mu,'d_sigma':sigma,'z_score':z,'p_value_single_tail':pval_single,
            'p_value_double_tail':pval_double, 'raw_amspl':d,'dist':distribution}



def proximity(net,T,S,D,null_model = 'degree_match',node_bucket = None, n_iter=1000,bin_size=100,
              symmetric=False):
    """
   Calculates the proximity between two sets of nodes in a given graph, based on the approach described by Guney et al., 2016.
   The method computes either the average shortest path length (ASPL) or its symmetrical version (SASPL) between two sets of nodes.

   The function first verifies if the network is connected. If it contains more than one connected component, a ValueError is raised.
   It also checks for the existence of all nodes in sets T and S within the network. If any nodes are missing, it issues a warning
   and proceeds with the existing nodes.

   Parameters:
   ------------
   net : networkx.Graph
       The input graph for which pairwise proximities will be computed.

   T : Iterable (list, set, etc.)
       A collection of 'source' nodes.

   S : Iterable (list, set, etc.)
       A collection of 'target' nodes.

   D : DistanceMatrix
       A precomputed distance matrix where D[i][j] provides the shortest distance between nodes i and j.
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
       Number of iterations/samples for assessing significance. Default is 1000.

   bin_size : int, optional
       Determines the size of the logarithmic bins when using the 'log-binning' method. Default is 100.

   symmetric : bool, optional
       If True, computes the symmetrical version of proximity using SASPL; otherwise, uses ASPL. Default is False.

   Returns:
   ---------
   dict
       A dictionary containing various statistics related to proximity, including:
       - 'd_mu': The average distance in the randomized samples.
       - 'd_sigma': The standard deviation of distances in the randomized samples.
       - 'z_score': The z-score of the actual distance in relation to the randomized samples.
       - 'p_value_single_tail': One-tail P-value associated with the proximity z-score
       - 'p_value_double_tail': Two-tail P-value associated with the proximity z-score
       - 'p_val': P-value associated with the z-score.
       - 'raw_amspl': The raw average minimum shortest path length between sets T and S.

       - 'dist': A list containing distances from each randomization iteration.

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
       - If any elements in `T` or `S` do not exist in the network.

   Notes:
   -------
   - The network should be connected to obtain meaningful proximity values. Disconnected components may skew results.
   - Proximity is based on the paper by Guney et al. 2016 (doi:10.1038/ncomms10331).
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
        Determines the size of the logarithmic bins when using 'log_binning'. Relevant only if
        `null_model` is set to 'log_binning'. Default is 100.

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
    tail = np.abs(tail)
    tail = tail[tail >= s]
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


def lcc_significance(net, A, null_model='degree_match', node_bucket = None, n_iter=1000,bin_size=100):
    """
    Calculate the statistical significance of the size of the Largest Connected Component (LCC)
    of a subgraph induced by the node set `A` in the network `net`.

    This function generates a null model distribution for the LCC size by resampling nodes from the
    network while preserving their degrees. The statistical significance of the observed LCC size is
    then determined by comparing it against this null model distribution.

    Parameters:
    -----------
    net : networkx.Graph
        The input network.

    A : list or set
        The set of nodes for which the LCC is to be determined.

    null_model : str, optional (default='degree_match')
        The method used for generating the null model. Can be 'degree_match', 'log_binning',
        'uniform', or 'custom'.

    node_bucket : dictionary, optional
        A collection of nodes to be used in 'custom' mode, mandatory when the null_model is set to 'custom'.
        This parameter should be a dictionary where each key represents a node ('node_k') from the network,
        and the corresponding value is a list of alternative nodes ('proxy_i').
        These alternatives are used by the null model for resampling:
        node_bucket[node_k] = [proxy_1, proxy_2, ..., proxy_m].
        Here, each 'proxy_i' serves as a potential substitute to be sampled in place of 'node_k'.


    n_iter : int, optional (default=1000)
        Number of iterations for generating the null model distribution.

    bin_size : int, optional (default=100)
        Size of bins if 'log_binning' method is used.

    Returns:
    ----------
    dict :
        A dictionary containing:
            - 'd_mu': Mean of the null model LCC size distribution.
            - 'd_sigma': Standard deviation of the null model LCC size distribution.
            - 'z_score': The z-score of the observed LCC size.
            - 'p_val': The p-value corresponding to the z-score.
            - 'lcc': Nodes in the largest connected component of `A`.
            - 'lcc_size': Size of the largest connected component of `A`.
            - 'dist': The null model LCC size distribution.

    Raises:
    ---------
    ValueError:
        - If 'n_iter' is less than or equal to 0.
        - If 'bin_size' is less than 1 when 'log_binning' or 'strength_binning' is used.
        - If 'null_model' is not one of ['degree_match', 'log_binning', 'strength_binning', 'uniform', 'custom'].
        - If 'node_bucket' is not provided when 'null_model' is 'custom'.

    Warnings:
    -----------
    UserWarning:
        - If any elements in `A` do not exist in the network.

    Notes:
    --------
    - Ensure the network does not contain any isolated nodes.
    """

    sA = set(A)
    set_a = sA & set(net.nodes())

    if len(sA.difference(set_a)) > 0:
        warnings.warn("A contains elements that are not present in the network.", UserWarning)
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
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
    elif null_model == 'uniform':
        bucket = set(net.nodes).difference(set_a)
    elif null_model == 'custom':
        bucket = node_bucket
        if node_bucket == None:
            raise ValueError("In custom mode, node_bucket must be provided")
    else:
        raise ValueError("Null model should be in ['degree_match'|'log_binning'|'uniform'|'custom']")

    lcc = extract_lcc(set_a,net)

    distribution = []
    for i in range(n_iter):
        if null_model != 'uniform':
            rs = _sample_preserving_degrees(net, set_a, bucket)
        else:
            rs = random.sample(list(bucket), len(set_a))
        sub = extract_lcc(rs,net)
        distribution.append(len(sub))

        if i%100 == 0:
            print(f"\rIter {i} of {n_iter}",end="")
    print("")

    l_lcc = len(lcc)
    mu = np.mean(distribution)
    sigma = np.std(distribution)
    z = (l_lcc - mu) / sigma
    distribution = np.array(distribution)
    S = distribution[distribution >= l_lcc]
    pval = len(S)/len(distribution)
    return {'d_mu':mu,'d_sigma':sigma,'z_score':z,'p_val':pval,'lcc':lcc, 'lcc_size':l_lcc,'dist':distribution}


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

@ray.remote
def _calculate_score(source, target, sources, targets, ppi, distance_matrix,
                    score, properties, null_model, node_bucket,n_iter,bin_size,symmetric):

    source_nodes = sources[source]
    target_nodes = targets[target]

    if score =="proximity":
        scores = proximity(ppi, source_nodes, target_nodes, distance_matrix,
                        null_model=null_model,node_bucket=node_bucket,n_iter=n_iter,
                        bin_size=bin_size,symmetric=symmetric)

        results = _to_dictionary(scores,properties)
    elif score =="separation_z_score":
        scores = separation_z_score(ppi, source_nodes, target_nodes, distance_matrix,
                        null_model=null_model,node_bucket=node_bucket,n_iter=n_iter,
                        bin_size=bin_size)

        results = _to_dictionary(scores,properties)
    elif score == "separation":
        scores = separation(ppi,source_nodes,target_nodes,distance_matrix)
        results = {"raw_separation":scores}

    print(f"{source}-{target} finished")

    return (source,target,results)


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


def screening(sources,targets, network, distance_matrix, score="proximity", properties=["z_score"],
              null_model= 'degree_match', node_bucket = None, n_iter=1000,bin_size=100,symmetric=False,
              n_procs=None):
    """
    Screens for relationships between sets of source and target nodes within a given network,
    evaluating proximity or separation. This function facilitates drug repurposing and other network
    medicine applications by allowing the assessment of network-based relationships.

    Parameters:
    ------------
    sources : dict
        A dictionary where keys are identifiers (e.g., drug names) and values are sets of nodes
        (e.g., drug target genes) representing source entities.

    targets : dict
        A dictionary where keys are identifiers (e.g., disease names) and values are sets of nodes
        (e.g., disease genes) representing target entities.

    network : networkx.Graph
        The graph representing the network within which the analysis is conducted. Can be weighted or unweighted.

    distance_matrix : DistanceMatrix
        A precomputed matrix providing the shortest distance between node pairs. Should be generated
        using `all_pair_distances` or a similar method.

    score : str, optional
        The metric for comparison. Options: "proximity", "separation_z_score", "separation". Default is "proximity".

    properties : list, optional
        The properties to retrieve from the analysis based on the `score`. Options vary based on the selected `score`:
        - For "proximity": 'z_score', 'p_value_single_tail', 'p_value_double_tail', 'raw_amspl'
        - For "separation_z_score": 'z_score', 'p_value_single_tail', 'p_value_double_tail', 'raw_separation'
        - For "separation": 'raw_separation'
        Default is ["z_score"].

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
        Number of iterations/samples for assessing significance. Default is 1000.

    bin_size : int, optional
        Determines the size of the logarithmic bins when using the 'log-binning' method. Default is 100.

    symmetric : bool, optional
        If True, computes the symmetrical version of proximity using SASPL; otherwise, uses ASPL. Default is False.

    n_procs : int, optional
        Number of processors to use. Defaults to the number of CPUs available if None.

    Returns:
    ---------
    dict of pd.DataFrame
        A dictionary where each key corresponds to a property from the `properties` list and the value is a DataFrame.
        Each DataFrame is indexed by source entities with columns representing target entities, populated with
        the values of the specified `score` for the corresponding property.

    Raises:
    --------
    ValueError
        If the network is not connected, if n_iter <= 0, if bin_size < 1 for certain null models,
        or if the null model or score specified is invalid.

    Notes:
    -------
    This function utilizes Ray for parallel processing to efficiently compute the desired metrics across
    multiple source-target pairs.
    """


    if nx.number_connected_components(network) > 1:
        raise ValueError("The network is not connected (it contains more than one connected component)")

    if n_iter <= 0:
        raise ValueError("n_iter must be greater than 0")


    if null_model == 'log_binning':
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
    elif null_model == 'strength_binning':
        if bin_size < 1:
            raise ValueError("bin_size must be greater or equal than 1")
    elif null_model == 'custom':
        if node_bucket == None:
            raise ValueError("In custom mode, node_bucket must be provided")
    else:
        if null_model not in ['degree_match','uniform']:
            raise ValueError("Null model should be: 'degree_match'|'log_binning'|'uniform'|'custom'")


    if score not in ["proximity","separation_z_score","separation"]:
        raise ValueError("Score should be: 'proximity','separation_z_score','separation'")

    if score == "proximity":
        if not set(properties) <= {"z_score","p_value_single_tail","p_value_double_tail","raw_amspl"}:
            raise ValueError("prop attribute should be: 'z_score','p_value_single_tail','p_value_double_tail','raw_amspl'")
    elif score == "separation_z_score":
        if not set(properties) <= {"z_score","p_value_single_tail","p_value_double_tail","raw_separation"}:
            raise ValueError("prop attribute should be: 'z_score','p_value_single_tail','p_value_double_tail','raw_separation'")


    valid_sources = _check_dictionary_integrity(network, sources)
    valid_targets = _check_dictionary_integrity(network, targets)

    source_names = valid_sources.keys()
    target_names = valid_targets.keys()

    if n_procs == None:
        num_cpus = os.cpu_count()
    else:
        num_cpus = n_procs

    ray.shutdown()
    ray.init(num_cpus = num_cpus)

    net_ref = ray.put(network)
    sources_ref = ray.put(valid_sources)
    targets_ref = ray.put(valid_targets)
    node_bucket_ref = ray.put(node_bucket)
    distance_matrix_ref = ray.put(distance_matrix)

    futures = []

    for source in source_names:
        for target in target_names:
            future = _calculate_score.remote(source, target, sources_ref, targets_ref,
                                     net_ref, distance_matrix_ref, score, properties,
                                     null_model, node_bucket_ref,n_iter,bin_size,
                                     symmetric)
            futures.append(future)

    results = ray.get(futures)

    ray.shutdown()

    dict_tables = _to_dict_tables(results, properties, target_names, source_names)

    return dict_tables
