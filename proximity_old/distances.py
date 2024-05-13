#! /usr/bin/env python3
"""
distances.py contains the functions of network-based distance measures.

creator: rodrigo dorantes gilardi (rodgdor@gmail.com)
date: 03-16-2022
"""
from collections import defaultdict

import numpy as np

try:
    import graph_tool.all as gt
except ModuleNotFoundError:
    pass
try:
    import networkx as nx
except ModuleNotFoundError:
    pass


def proximity(net, T, S, bin_size=100, n_iter=1000) -> dict:
    """Return the proximity between two sets of nodes.
    
    The proximity between two sets of nodes as defined in
    the paper by Guney et al. 2016.

    doi: 10.1038/ncomms10331
    """
    lower, upper, values = get_binning(net, bin_size=bin_size)
    d_c = get_min_shortest_paths(net, T, S)
    distribution = []
    for _ in range(n_iter):
        ran_T  = select_random_nodes(net, T, lower, upper, values)
        ran_S  = select_random_nodes(net, S, lower, upper, values)
        ran_d_c = get_min_shortest_paths(net, ran_T, ran_S)
        distribution.append(ran_d_c)
    mu = np.mean(distribution)
    sigma = np.std(distribution)
    z = (d_c - mu) / sigma

    return {'d_c': d_c, 'z_score': z, 'mu': mu, 'sigma': sigma}


def separation(net, A, B):
    """Calculate the separation [1]_ between node sets A and B.

    Parameters
    ----------
    net: prximity.Network
    A: container
        A subset of nodes in net.Graph
    B: container
        A subset of nodes in net.Graph

    Returns
    -------
    sep: the separation of nodes A and B in net.Graph

    References
    ----
    .. [1] Menche, Jörg, et al.
        Uncovering disease-disease relationships through the
        incomplete interactome.
        Science 347.6224 (2015).
    """
    s_aa = get_avg_shortest_paths(net, A, A)
    s_bb = get_avg_shortest_paths(net, B, B)
    s_ab = get_avg_shortest_paths(net, A, B)

    sep = s_ab - (s_aa + s_bb) / 2

    return sep


def get_binning(net, bin_size=100):
    """Return a histogram of the degrees of the PPI.
    
    The histogram should have bins with size at least equal to
    `bin_size`. For each bin, the bin bounds l and u should be optimized
    such that a bin with bounds l and u - 1 does is of size smaller
    than `bin_size`.

    Note that for a node with degree d to be in a bin with bounds (l, u],
    we have that l < d <= u.

    Parameters
    ----------
    net: proximity.Network
        Usually the protein-protein interaction network. **NOTE** If `net` is a
        gt.Graph instance, it should have a `gt.VertexPropertyMap` with the
        node names caled "ids".
    bin_size: int

    Returns
    -------
    nodes: list
        The nodes of each bin.
    lower: list
        The lower bound of each bin
    upper: list
        The upper bound of each bin.
    """
    graph = net.Graph
    degree2nodes = {}
    if net.module == 'gt':
        try:
            ids = graph.vertex_properties['ids']
        except KeyError:
            raise Exception(
                "The graph should have a vertex property called 'ids'!")
        deg = graph.degree_property_map('out')
        for v in graph.iter_vertices():
            degree2nodes.setdefault(deg[v], list())
            degree2nodes[deg[v]].append(ids[v])

    elif net.module == 'nx':
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
    
    return lower, upper, nodes
    

def select_random_nodes(net, nodes, lower, upper, values) -> list:
    """Return an array with a degree-preserving random selection of nodes.

    **NOTE** If `net` is a gt.Graph instance, it should have a
    `gt.VertexPropertyMap` with the node names caled "ids".

    Parameters
    ----------
    net: proximity.Network
    nodes: list
        Array of reference nodes.
    lower: list
        Lower bounds of histogram of the histogram
    upper: list
        Upper bounds of histogram of the histogram
    values: list
        Nodes in each bin of the histogram
    Returns
    -------
    rand: list
        A list of degree-preserving nodes of size `len(nodes)`
    """
    graph = net.Graph
    graph_tool = 0
    if net.module == 'gt':
        graph_tool = 1
        try:
            ids = graph.vertex_properties['ids']
        except KeyError:
            raise Exception(
                "The graph should have a vertex property called 'ids'!")
        node2id = {}
        for v in graph.iter_vertices():
            node2id[ids[v]] = v

    reference_degrees = defaultdict(set)
    for node in nodes:
        if graph_tool:
            v = graph.vertex(node2id[node])
            d = v.out_degree()
        else:
            d = graph.degree(node)
        reference_degrees[d].add(node)
    sample = []
    for d in sorted(reference_degrees):
        n = len(reference_degrees[d])
        for i in range(len(values)):
            if lower[i] < d <= upper[i]:
                ref = list(set(values[i]) - set(reference_degrees[d]))
                break
        sample.extend(np.random.choice(ref, n, replace=False))
    
    return sample
        

def get_min_shortest_paths(net, T, S) -> float:
    """Get the minimal shortest path lengths between two sets of nodes.
    
    Parameters
    ----------
    net: proximity.Network
    T, S: Container
        Each contains nodes from net.Graph, the drug targets and disease
        genes, respectively.
    Returns
    -------
    d_c: float
    """
    T = set(T)
    S = set(S)
    # Are sets not identical?
    identical_sets = (T == S)
    if identical_sets:
        return 0
    graph = net.Graph
    # Set the function to calculate shortest paths
    if net.module == 'gt':
        ids = graph.vertex_properties['ids']
        labels = {ids[v]: v for v in graph.iter_vertices()}
        # A property map will define the value for paths that
        # does not exist. T 32-bit int has a max value of 2147483647.
        pm = graph.new_vp('int32_t')
        distance = lambda x, y, z: gt.shortest_distance(x, labels[y], labels[z], dist_map=pm)
    elif net.module == 'nx':
        distance = nx.shortest_path_length
    # For each node, get its min distance with nodes from the other group
    min_distance = defaultdict(lambda: float('inf'))
    for a in T:
        for b in S:
            if a != b:
                # Sometimes there is no path between `a` and `b`
                try:
                    spl = distance(graph, a, b)  # gt will return `inf`
                except nx.NetworkXNoPath:
                    continue
            else:
                spl = 0
            # Update min distance
            min_distance[a] = min(spl, min_distance[a])

    # Get average of all min distances
    bad_values = (float('inf'), np.nan, 2147483647)
    min_lengths = [x for x in min_distance.values() if x not in bad_values]
    if min_lengths:
        d_c = np.mean(min_lengths)  # This varname comes the original paper of 2016
    else:
        raise Exception("The two sets are different connected components!")
    
    return d_c


def get_avg_shortest_paths(net, A, B):
    """Get average shortest paths between two sets of nodes.
    
    Parameters
    ----------
    net: proximity.Network
    A, B: Container
        Each contains nodes from net.Graph

    Returns
    -------
    S_AB: float
    """
    A = set(A)
    B = set(B)
    # Are sets not identical?
    different_sets = (A != B)
    graph = net.Graph
    # Set the function to calculate shortest paths
    if net.module == 'gt':
        # A property map will define the value for paths that
        # does not exist. A 32-bit int has a max value of 2147483647.
        pm = graph.new_vp('int32_t')
        distance = lambda x, y, z: gt.shortest_distance(x, y, z, dist_map=pm)
    elif net.module == 'nx':
        distance = nx.shortest_path_length
    # For each node, get its min distance with nodes from the other group
    min_distance = defaultdict(lambda: float('inf'))
    for a in A:
        for b in B:
            # zeros in same-group nodes are ignored
            if different_sets or a != b:
                # Sometimes there is no path between `a` and `b`
                try:
                    spl = distance(graph, a, b)  # gt will return `inf`
                except nx.NetworkXNoPath:
                    continue
                # Update min distance
                min_distance[a] = min(spl, min_distance[a])
                min_distance[b] = min(spl, min_distance[b])
    
    # Get average of all min distances
    bad_values = (float('inf'), np.nan, 2147483647)
    min_lengths = [x for x in min_distance.values() if x not in bad_values]
    S_AB = np.mean(min_lengths)
    
    return S_AB