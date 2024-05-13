# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 15:25:20 2023

@author: aalda
"""

import numpy as np
import matplotlib.pyplot as plt
import tools.LoadContext as context

max_nodes = 18000
edges = (max_nodes**2)*(500000/(18000**2))

set_nodes = 1000
x = np.arange(0, set_nodes,set_nodes/100)
k = (edges + max_nodes)*np.log2(max_nodes)


k = 3

old = k*(x**2)
new = (x**2) + x*k


plt.figure()
plt.plot(x,old,label="$x^2(E + V)log(V)$",color="black")
plt.plot(x,new,label = "$x^2 + (E + V)log(V)x$",color="blue")
plt.xlabel("Nodes in A and B")
plt.ylabel("Operations")
plt.legend(loc="best")
# plt.yscale("log")


