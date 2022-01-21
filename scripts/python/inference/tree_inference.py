import random

import numpy as np
import networkx as nx

from Bio import SeqIO
from networkx import Graph
from networkx.algorithms.tree.mst import minimum_spanning_tree
from networkx.algorithms.components import connected_components
from pathlib import PosixPath

from scripts.python.helpers import compute_jc_distance



def compute_MST(fasta: PosixPath):
    #FIXME: Add comments and description of the method

    # get sequences and sequence names from the provided fasta file
    seq_names = []
    sequences = []
    for record in SeqIO.parse(fasta, 'fasta'):
        seq_names.append(record.id)
        sequences.append(record.seq)

    # compute the distance matrix from the sequences(jc distance)
    n_taxa = len(seq_names)
    distance_matrix = np.zeros(shape=(n_taxa, n_taxa))
    for i, seq1 in enumerate(sequences[:-1]):
        for j, seq2 in enumerate(sequences[i+1:]):
            jc_dist = compute_jc_distance(seq1, seq2)
            distance_matrix[[i, j+i+1], [j+i+1, i]] = jc_dist

    # create networkx graph from distance matrix
    graph = nx.from_numpy_matrix(A=distance_matrix, parallel_edges=False)

    #draw_network(graph)
    # rename graph nodes
    name_mapping = {old: new for old, new in zip(range(n_taxa), seq_names)}
    nx.relabel_nodes(G=graph, mapping=name_mapping, copy=False)

    # create MST from graph
    mst = minimum_spanning_tree(graph, algorithm='prim')    # use prim, as the graph is dense
    return mst


def remove_edges(mst: Graph, n_edges: int = None, fraction: float = None, method='random'):
    """ Removes n edges from a graph. The method returns the connected components of the resulting graph (as Graph() objects)
    Either the number of edges (n_edges) or the fraction of edges (fraction) to remove has to be passed. If both are passed, n_edges is being used.
    Possible approaches to remove edges include:
        - 'random': random selection of edges
        - 'random_weighted': random selection of edges. Edges with greater weight have higher probability to be selected
        - 'by_weight': edges with the maximum weight are selected
    :param mst:
    :param n_edges:
    :param fraction:
    :param method: choose 'random', 'random_weighted' or 'by_weight'
    :return: Returns a list of Graph() objects
    """

    # FIXME: Make sure the output graphs are not too small (e.g. at least 3 (or n) taxa per graph)

    if n_edges is None:
        if fraction is None:
            raise ValueError(f"Either 'n_edges' or 'fraction' has to be provided")
        elif not 0 < fraction < 1:
            raise ValueError(f"'fraction' has to lie in the interval (0, 1). But fraction={fraction}")
        else:
            n_edges = round(fraction * len(mst.edges), 0)

    # select edges to remove
    if method == 'random':
        edges_to_remove = random_selection(graph=mst, n=n_edges)
    elif method == 'random_weighted':
        edges_to_remove = weighted_random_selection(graph=mst, n=n_edges)
    elif method == 'by_weight':
        edges_to_remove = select_by_weight(graph=mst, n=n_edges)
    else:
        raise ValueError(f"'method' has to be one of the following choices: ['random', 'random_weighted', 'weight'")

    # remove edges
    for edge in edges_to_remove:
        mst.remove_edge(edge[0], edge[1])

    # create a subgraph from each connected component
    sub_mst = [mst.subgraph(c).copy() for c in nx.connected_components(mst)]

    return sub_mst

def random_selection(graph: Graph, n):
    """ Randomly selects n edges from a graph """
    edges = list(graph.edges)
    return random.sample(edges, n)

def weighted_random_selection(graph: Graph, n):
    """ Randomly selects n edges from a graph.
    The greater the weight of the edge, the greater the probability it will be selected
    """
    idx = np.arange(graph.size())
    size = graph.size(weight='weight')      # sum of all edge weights

    # compute probabilities for each edge
    edges = []
    p = []
    for edge, values in graph.edges.items():
        edges.append(edge)
        p.append(values['weight'] / size)

    # randomly choose edges to remove (weighted by edge weight)
    idx_to_remove = np.random.choice(a=idx, size=n, p=p, replace=False)
    edges_to_remove = [edges[i] for i in idx_to_remove]

    return edges_to_remove

def select_by_weight(graph: Graph, n):
    """ Select n edges with the highest weights """
    # get edges and their weight
    edges = []
    weights = []
    for edge, values in graph.edges.items():
        edges.append(edge)
        weights.append(values['weight'])

    # get indices of the n highest weights
    idx_top = np.argpartition(np.array(weights), -n)[-n:]
    edges_to_remove = [edges[i] for i in idx_top]

    return edges_to_remove