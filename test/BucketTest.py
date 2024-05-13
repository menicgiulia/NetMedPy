# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 14:49:52 2023

@author: aalda
"""


import netmedpy.NetMedPy as netMet
import tools.LoadContext as context
import random
import tools.Cronometer as cronometer
import pandas as pd

# if __name__=="__main__":

#     c = context.Context()
#     net = c.PPI
#     bucket = netMet.group_nodes_by_degree(net)

#     cr = cronometer.Cronometer()
#     degrees = dict(net.degree())


#     df = pd.DataFrame(list(degrees.items()), columns=["Node", "Degree"])
#     df = df.sort_values(by='Degree',ascending=False)
#     available = set(net.nodes) - {'NRF1','APP'}


#     for j in range(1000):
#         S = set(random.sample(available, 200))
#         S = S | {'NRF1','APP'}
#         l = []
#         cr.tick()
#         for i in range(10000):
#             if i%1000 == 0:
#                 print(f"\rIteration {i}          ",end="")

#             S2 = netMet._sample_preserving_degrees(net,S,bucket)
#             l.append(len(S2))
#         se = set(l)
#         cr.tock()
#         print(f"Sizes {se}  {cr.elapsed_seconds}")





if __name__=="__main__":

    c = context.Context()
    net = c.PPI

    lower, upper, nodes = netMet._get_degree_binning(net,bin_size=100)
    bucket = netMet._dictionary_from_binning(lower, upper, nodes)

    cr = cronometer.Cronometer()
    degrees = dict(net.degree())


    df = pd.DataFrame(list(degrees.items()), columns=["Node", "Degree"])
    df = df.sort_values(by='Degree',ascending=False)
    available = set(net.nodes) - {'NRF1','APP'}



    for j in range(1000):
        S = set(random.sample(available, 200))
        S = S | {'NRF1','APP'}
        l = []
        cr.tick()
        for i in range(10000):
            if i%1000 == 0:
                print(f"\rIteration {i}          ",end="")

            S2 = netMet._sample_preserving_degrees(net,S,bucket)
            l.append(len(S2))
        se = set(l)
        cr.tock()
        print(f"Sizes {se}  {cr.elapsed_seconds}")
