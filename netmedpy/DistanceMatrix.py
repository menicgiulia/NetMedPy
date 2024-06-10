# -*- coding: utf-8 -*-
"""This module provides a DistanceMatrix class designed to help manage distances
between nodes in a network, applicable in graph-theoretical algorithms and network
analysis research.

Required packages:
    - numpy

Authors:
    - Andres Aldana Gonzalez (a.aldana@northeastern.edu)
"""

import numpy as np


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
