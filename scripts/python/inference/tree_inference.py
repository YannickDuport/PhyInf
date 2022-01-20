import numpy as np
import networkx as nx

from Bio import SeqIO
from networkx.algorithms.tree.mst import minimum_spanning_tree
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