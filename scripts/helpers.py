import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from pathlib import PosixPath
from networkx.algorithms.tree.mst import minimum_spanning_tree


def get_sequences_from_fasta(fasta: PosixPath):
    """ Gets sequences and their name from a given fasta file
    :param fasta:
    :return: A list of sequence names, as well as a list with the corresponding sequence. Both lists are in the same order
    """
    seq_names = []
    seq = []
    for record in SeqIO.parse(fasta, 'fasta'):
        seq_names.append(record.id)
        seq.append(record.seq)
    return seq_names, seq


def compute_jc_distance(seq1, seq2):
    # compute hamming distance between seq1 and seq2
    hamming_dist = sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2)) / len(seq1)

    # compute jc corrected distance
    if hamming_dist >= 0.75:
      return np.inf
    else:
      return -3/4 * np.log(1 - 4/3 * hamming_dist)

def compute_distance_matrix(sequences):

    # Compute distance matrix
    n_seq = len(sequences)
    distance_matrix = np.zeros(shape=(n_seq, n_seq))
    for i, seq1 in enumerate(sequences[:-1]):
        for j, seq2 in enumerate(sequences[i+1:]):
            jc_dist = compute_jc_distance(seq1, seq2)
            distance_matrix[[i, j+i+1], [j+i+1, i]] = jc_dist
    return distance_matrix

def compute_MST(fasta: PosixPath):

    # get sequences and compute the corresponding distance matrix (jc distance)
    seq_names, sequences = get_sequences_from_fasta(fasta)
    n_taxa = len(seq_names)
    dist_matrix = compute_distance_matrix(sequences)


    # create networkx graph from distance matrix
    graph = nx.from_numpy_matrix(A=dist_matrix, parallel_edges=False)

    #draw_network(graph)
    # rename graph nodes
    name_mapping = {old: new for old, new in zip(range(n_taxa), seq_names)}
    nx.relabel_nodes(G=graph, mapping=name_mapping, copy=False)

    # create MST from graph
    mst = minimum_spanning_tree(graph, algorithm='prim')    # use prim, as the graph is dense
    return mst

def draw_network(graph):
    fig = plt.figure()
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos)
    nx.draw_networkx_labels(graph, pos)
    nx.draw_networkx_edge_labels(graph, pos)
    fig.show()

