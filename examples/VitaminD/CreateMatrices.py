# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 23:20:16 2024

@author: aalda
"""
import networkx as nx
import pandas as pd

import netmedpy.NetMedPy as netmedpy
import pickle
import os

ppi = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/channing/ppi932023.pkl"

with open(ppi,'rb') as file:
    ppi = pickle.load(file)

ppi = netmedpy.extract_lcc(ppi.nodes, ppi)

baseMat = "D:/data/matDist/channing/ppi932023/"

#'shortest_path'
types = ["random_walk","biased_random_walk","communicability"]


for t in types:
    D = netmedpy.all_pair_distances(ppi, distance=t,reset=0.2,n_processors=40,n_tasks=40)
    netmedpy.save_distances(D, f"{baseMat}{t}_N2.pkl" )
