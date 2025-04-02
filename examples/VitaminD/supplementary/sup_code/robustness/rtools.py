import pickle
import networkx as nx

def save(obj, file):
    with open(file, "wb") as file:
        pickle.dump(obj, file)


def load(file):
    with open(file, "rb") as file:
        return pickle.load(file)


def filter_for_ppi(ppi, gene_set):
    res = {}
    nodes = set(ppi.nodes)
    for n,s in gene_set.items():
        new_s = s & nodes
        res[n] = new_s
    return res