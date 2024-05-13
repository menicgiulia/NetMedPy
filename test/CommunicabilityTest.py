# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 01:45:04 2023

@author: aalda
"""

import networkx as nx
import tools.Cronometer as Cronometer

# Create a graph
G = nx.Graph()
# Add some edges to the graph
G.add_edges_from([('a', 'b'), ('a', 'c'), ('b', 'd'), ('c', 'd')])

# Calculate communicability
comm = nx.communicability(G)

# comm is a dictionary where keys are node pairs and values are communicability scores
#print(comm)

comm_exp = nx.communicability_exp(G)

#print(comm_exp)


# Create a scale-free network using the Barabási-Albert model
# Starting with a small number (m0) of initial nodes and adding nodes one by one
m0 = 5  # Initial number of nodes
m = 2   # Number of edges to attach from a new node to existing nodes

# Generate the graph
G_ba = nx.barabasi_albert_graph(100, m, seed=42)

print("Calculating communicability")

c = Cronometer.Cronometer()

c.tick()

comm_exponential = nx.communicability_exp(G)

c.tock()

print(c.elapsed_milliseconds())
print(c.format_seconds())
# print(comm_exponential[10][20])

print(comm_exponential)
print(list(comm_exponential.keys()))





