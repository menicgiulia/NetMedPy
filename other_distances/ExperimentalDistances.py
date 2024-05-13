# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:32:12 2023

@author: aalda
"""

import networkx as nx
import numpy as np
import other_distances.ShortestPaths as sp

def get_avg_min_shortest_path_oa_nohalt(net,A,B):
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
        #distances = nx.single_source_shortest_path_length(net,n)
        distances = nx.shortest_path_length(net,n)


        for m in valid_b:
            db = distances.get(m)
            if(db != None):
                min_dist = min(min_dist,db)
        min_distances[idx]=min_dist
        idx+=1
    avg = min_distances.mean()
    return avg




def get_avg_min_shortest_path_oa_halt(net,A,B):
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
        #distances = nx.single_source_shortest_path_length(net,n)
        distances = sp.modified_dijkstra(net, n, valid_b)

        for m in valid_b:
            db = distances.get(m)
            if(db != None):
                min_dist = min(min_dist,db)
        min_distances[idx]=min_dist
        idx+=1
    avg = min_distances.mean()
    return avg
