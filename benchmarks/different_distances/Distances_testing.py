# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 03:36:28 2023

@author: aalda
"""
import networkx as nx
import netmedpy.NetMedPy as network_metrics
import netmedpy.DistanceMatrix as DistanceMatrix


if __name__=="__main__":
    # Create a scale-free network using the Barabási-Albert model
    # Starting with a small number (m0) of initial nodes and adding nodes one by one
    m0 = 5  # Initial number of nodes
    m = 2   # Number of edges to attach from a new node to existing nodes

    # Generate the graph
    G_ba = nx.barabasi_albert_graph(1000, m, seed=42)

    #D1 = network_metrics.all_pair_distances(G_ba, distance="shortest_path")

    #D2 = network_metrics.all_pair_distances(G_ba, distance="random_walk",restart = 0.2)

    D3 = network_metrics.all_pair_distances(G_ba,distance="biased_random_walk", restart = 0.2)

    #D4 = network_metrics.all_pair_distances(G_ba, distance="communicability")

    print(D3.get(0, 1))
