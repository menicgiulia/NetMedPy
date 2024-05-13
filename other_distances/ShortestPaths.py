# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 00:43:23 2023

@author: aalda
"""

import networkx as nx
import heapq

def modified_dijkstra(graph, start, target_set):
    """
    Implements Dijkstra's algorithm and halts once all nodes in the target_set have been visited.

    Parameters:
    - graph (networkx.Graph): the input graph.
    - start (node): the starting node.
    - target_set (set of nodes): a set of nodes to visit.

    Returns:
    - distances (dict): A dictionary of shortest paths from the start node to all other nodes.
    """

    visited = set()
    distances = {node: float('infinity') for node in graph}
    distances[start] = 0
    priority_queue = [(0, start)]

    # Create a copy of the target set to track remaining nodes to visit
    remaining_targets = target_set.copy()

    while priority_queue:
        current_distance, current_node = heapq.heappop(priority_queue)

        # If the current node has already been visited, skip
        if current_node in visited:
            continue

        visited.add(current_node)
        if current_node in remaining_targets:
            remaining_targets.remove(current_node)

        # If all the nodes in the target_set have been visited, break
        if not remaining_targets:
            break

        for neighbor, edge_attributes in graph[current_node].items():
            distance = edge_attributes.get('weight', 1)  # Assuming a default weight of 1 if no weight is specified
            new_distance = current_distance + distance

            if new_distance < distances[neighbor]:
                distances[neighbor] = new_distance
                heapq.heappush(priority_queue, (new_distance, neighbor))

    return distances

# Example usage:
# if __name__ == "__main__":
#     G = nx.Graph()
#     G.add_weighted_edges_from([(1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1)])
#     start_node = 1
#     target_nodes = {4, 5}
#     distances = modified_dijkstra(G, start_node, target_nodes)
#     print(distances)
