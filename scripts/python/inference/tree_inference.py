import nxmetis
import pickle

import networkx as nx

from Bio import SeqIO
from pathlib import PosixPath
from networkx.algorithms.tree.mst import minimum_spanning_tree

from numpy.typing import NDArray
from typing import Union, List

from scripts.python import timeit
from scripts.python.helpers import compute_pairwise_jc_distance
from scripts.python.inference.edge_selection import *


@timeit
def compute_distance_matrix(fasta: PosixPath, save: bool = False):
    """ Reads a fasta file and computes the pairwise jukes-cantor distance between the sequences
    If wanted, a tuple of sequence names and distance matrix (seq_names, distance_matrix) can be saved as binary file.
    The file will be called 'distance_matrix.bin' and stored at the same location as the fasta file.

    :param fasta: path to the fasta file
    :param save: set to True if you want to save the distance matrix
    :return: <(seq_names, distance_matrix)>
    """
    # get sequences and sequence names from the provided fasta file
    seq_names = []
    sequences = []
    for record in SeqIO.parse(fasta, 'fasta'):
        seq_names.append(record.id)
        sequences.append(record.seq)

    # compute the distance matrix from the sequences(jc distance)
    n_taxa = len(seq_names)
    distance_matrix = compute_pairwise_jc_distance(n_taxa, sequences)

    # save tuple of sequence names and distance matrix as binary file
    if save:
        outfile = fasta.parent / "distance_matrix.bin"
        print(f"Distance matrix saved to {outfile}")
        pickle.dump((seq_names, distance_matrix), outfile.open('wb'))

    return seq_names, distance_matrix


def compute_MST(distance_matrix: Union[PosixPath, NDArray], seq_names: List[str], model_name: str = None) -> Graph:
    """ Creates a MST from a given distance matrix.
    From the distance matrix an undirected, weighted and fully connected graph is constructed, using NetworkX.
    In the next step a MST is created from that graph, using NetworkX's method 'minimum_spanning_tree()'
    Since the MST is created from a fully connected graph, Prim's algorithm is used for this last step

    :param distance_matrix: Can be an n x n matrix, or a path to a binary file containing an n x n matrix
    :param seq_names: List of strings that will be used to name the nodes (size n)
    :param model_name: String that will be used to name the graph
    :return:
    """

    # if distance matrix is stored in a file, read it in
    if type(distance_matrix) == PosixPath:
        distance_matrix = np.loadtxt(distance_matrix)

    # create networkx graph from distance matrix
    graph = nx.from_numpy_matrix(A=distance_matrix, parallel_edges=False)
    graph.name = model_name

    # rename graph nodes
    n_taxa = distance_matrix.shape[0]
    name_mapping = {old: new for old, new in zip(range(n_taxa), seq_names)}
    nx.relabel_nodes(G=graph, mapping=name_mapping, copy=False)

    # create MST from graph
    mst = minimum_spanning_tree(graph, algorithm='prim')    # use prim, as the graph is dense
    return mst


def remove_edges(mst: Graph, n_edges: int = None, fraction: float = None, method='random', min_size=None) -> List[Graph]:
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
    :param min_size: Minimum number of nodes each resulting connected component must contain
    :return: Returns a list of Graph() objects
    """

    # FIXME: Make sure the output graphs are not too small (e.g. at least 3 (or n) taxa per graph)
    # FIXME: Add warning if method creates a small subgraph (<3 nodes)
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
        edges_to_remove = weighted_random_selection_iteratively(graph=mst, n=n_edges)
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




