# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:25:26 2023

@author: aalda
"""
def get_disease_genes(key,context):
    gda = context.GDA
    ppi_nodes = set(context.PPI.nodes)
    d = gda.query("NewName==@key")

    raw_genes = set(d.HGNC_Symbol.unique())

    ppi_genes = ppi_nodes & raw_genes

    return ppi_genes

def get_disease_names(context):
    gda = context.GDA
    return list(gda.NewName.unique())
