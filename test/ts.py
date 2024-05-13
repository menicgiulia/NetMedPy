# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 18:57:54 2024

@author: aalda
"""



import netmedpy.NetMedPy as netMet
import tools.LoadContext as context
import random
import tools.Cronometer as cronometer
import pandas as pd


if __name__=="__main__":

    c = context.Context()
    net = c.PPI
    GDA = c.GDA

    ra = GDA.query("NewName=='alzheimer disease'")
    ra = ra['HGNC_Symbol']

    ra = set(ra) & set(net.nodes)

    lcc = netMet.extract_lcc(ra, net)

    nl = set(lcc.nodes)

    sig = netMet.lcc_significance(net, nl,null_model='uniform',n_iter=10000)
    print(sig['p_val'])
