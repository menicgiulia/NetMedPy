# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 02:38:55 2023

@author: aalda
"""

import networkx as nx
import tools.Cronometer as Cronometer

# Create a graph
G = nx.Graph()
# Add some edges to the graph
G.add_edges_from([('a', 'b'), ('a', 'c'), ('b', 'd'), ('c', 'd')])




# Create a scale-free network using the Barabási-Albert model
# Starting with a small number (m0) of initial nodes and adding nodes one by one
m0 = 5  # Initial number of nodes
m = 2   # Number of edges to attach from a new node to existing nodes

# Generate the graph
G_ba = nx.barabasi_albert_graph(10, m, seed=42)

print("Calculating random walks")

start_node = 1
restart_prob = 0.2

c = Cronometer.Cronometer()

c.tick()

personalization = {node: 0 for node in G}
personalization[start_node] = 1

rw = nx.pagerank(G_ba, alpha=1-restart_prob, personalization=personalization)

c.tock()

print(c.elapsed_milliseconds())
print(c.format_seconds())
# print(comm_exponential[10][20])

print(rw)















 


