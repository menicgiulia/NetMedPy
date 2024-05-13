# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 00:24:11 2023

@author: aalda
"""

import tools.LoadContext as LoadContext
import netmedpy.NetMedPy as netMet
import netmedpy.DistanceMatrix as DistanceMatrix
from tools import Cronometer as Cronometer
import networkx as nx
import numpy as np

def log_binning(net):
    graph = net
    degree2nodes = {}

    deg = nx.degree(graph)
    for v, d in deg:
        degree2nodes.setdefault(d, list())
        degree2nodes[d].append(v)


def degree_match(graph):
    degree_dict = {}

    for node in graph.nodes():
        degree = graph.degree(node)

        if degree not in degree_dict:
            degree_dict[degree] = []

        degree_dict[degree].append(node)




if __name__=="__main__":

    CONTEXT = LoadContext.Context()

    PPI = CONTEXT.PPI



    c = Cronometer.Cronometer()


    lb = []

    for i in range(100):
        c.tick()
        log_binning(PPI)
        c.tock()
        lb.append(c.elapsed_milliseconds())


    co = []

    for i in range(100):
        c.tick()
        degree_match(PPI)
        c.tock()
        co.append(c.elapsed_milliseconds())

    lb_time = np.array(lb).mean()
    co_time = np.array(co).mean()
