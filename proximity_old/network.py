#! /usr/bin/env python3
"""
network.py contains the library-agnostic class Network.

creator: rodrigo dorantes gilardi (rodgdor@gmail.com)
date: 03-16-2022
"""
from proximity.distances import proximity, separation


class Network():
    def __init__(self, graph):
        """A module-agnostic network wrapper
        
        This class creates a network wrapper for graph_tool:graph and
        :networkx:Graph objects.

        Parameters
        ----------
        graph: :graph_tool:`~Graph` or :networkx:`~Graph`

        Raises
        ------
        ValueError if network not an instance of networkx or graph_tool.
        """

        error = ValueError("graph should be an instance of graph_tool.Graph "
                             "or networkx.Graph")
        try:
            module = graph.__module__
        except:
            raise error
        if module == "graph_tool":
            self.module = "gt"
            try:
                graph.vertex_properties['ids']
            except KeyError:
                raise Exception('Graph should have vertex property `ids`!')
        elif module == "networkx.classes.graph":
            self.module = "nx"
        else:
            raise error
        if graph.__class__.__name__ != 'Graph':
            raise error

        self.Graph = graph

    def get_proximity(self, T, S, bin_size=100, n_iter=1000):
        return proximity(self, T, S, bin_size=bin_size, n_iter=n_iter)

    def get_separation(self, A, B):
        """Get Separation of two sets of nodes."""

        return separation(self, A, B)