import copy
import random
import numpy as np
from networkx import Graph
from networkx.algorithms.components import connected_components
from scripts.python.helpers import timeit






def compute_viable_edges(graph: Graph, size: int = 4, edges=None):
    """ The method detects all edges that connect subgraphs of size 'size' or more
    Uses a brute force approach. Just checks for every single edge whether it's save to remove.
    Slow for large graphs!

    :param graph:
    :param size: The minimum size the subgraphs should have (i.e. minimum number of nodes)
    :param edges: A subset of edges to check. If none are provided, all edges of the graph will be checked
    :return: A list of edges that are save to remove
    """

    if edges is None:
        edges = graph.edges(data=True)

    # iterate through all edges, remove it, check size of the resulting partitions
    edge_pool = []
    for u, v, w in graph.edges(data=True):
        graph_cp = copy.deepcopy(graph)
        graph_cp.remove_edge(u, v)
        if all(len(g) >= size for g in connected_components(graph_cp)):
            edge_pool.append((u, v, w))

    return edge_pool


def random_selection(graph: Graph, n):
    """ Randomly selects n edges from a graph
    Edges that connect a node of degree 1 are removed from the pool of selectable edges """

    edge_pool = compute_viable_edges(graph)
    edge_pool = [(u, v) for u, v, w in edge_pool]
    return random.sample(edge_pool, n)


def weighted_random_selection(graph: Graph, n):
    """ Randomly selects n edges from a graph.
    The greater the weight of the edge, the greater the probability it will be selected
    """
    edge_pool = compute_viable_edges(graph)

    idx = np.arange(len(edge_pool))  # create index for viable each edge
    size = sum([w['weight'] for u, v, w in edge_pool])  # sum of all viable edges

    # compute probabilities for each edge
    edges = []
    p = []
    for u, v, w in edge_pool:
        edges.append((u, v))
        p.append(w['weight'] / size)

    # randomly choose edges to remove (weighted by edge weight)
    idx_to_remove = np.random.choice(a=idx, size=n, p=p, replace=False)
    edges_to_remove = [edges[i] for i in idx_to_remove]

    return edges_to_remove


def select_by_weight(graph: Graph, n):
    """ Select n edges with the highest weights """

    edge_pool = compute_viable_edges(graph)

    # get edges and their weight
    edges = []
    weights = []
    for u, v, w in edge_pool:
        edges.append((u, v))
        weights.append(w['weight'])

    # get indices of the n highest weights
    idx_top = np.argpartition(np.array(weights), -n)[-n:]
    edges_to_remove = [edges[i] for i in idx_top]

    return edges_to_remove

@timeit
def weighted_random_selection_iteratively(graph: Graph, n):
    """ Randomly selects n edges from a graph.
    The greater the weight of the edge, the greater the probability it will be selected

    Similar to weighted_random_selection(), just that the viable edges are computed after each split.
    Note: This has a high chance to fail if the s
    """
    graph_cp = copy.deepcopy(graph)

    edges_to_remove = []
    for i in range(n):
        # compute which edges are safe to remove (i.e. won't create too small partitions)
        edge_pool = compute_viable_edges(graph_cp, 10)  # FIXME: In subsequent iterations, not all edges have to be checked. Only those that were viable in the previous iteration
                                                        # i.e edges that weren't save to remove before, won't magically be ok now

        # Not possible to remove another edge without generating subgraphs smaller than the desired min size
        if edge_pool == []:
            return edges_to_remove

        idx = np.arange(len(edge_pool)) # create index for viable each edge
        size = sum([w['weight'] for u, v, w in edge_pool])  # sum of all viable edges

        # compute probabilities for each edge
        edges = []
        p = []
        for u, v, w in edge_pool:
            edges.append((u, v))
            p.append(w['weight'] / size)

        # randomly choose edges to remove (weighted by edge weight)
        idx_to_remove = np.random.choice(a=idx, size=1, p=p, replace=False)
        edge_to_remove = edges[idx_to_remove[0]]
        edges_to_remove.append(edge_to_remove)
        graph_cp.remove_edge(edge_to_remove[0], edge_to_remove[1])

    return edges_to_remove


""" Deprecaded methods (i.e. methods dont work as intended

def compute_viable_edges_cp(graph: Graph, size: int = 4):
    ''' The method detects all edges that connect subgraphs of size 4 or more
    Procedure:
        1) Detect all nodes of degree 1 in the graph. Let U be the set of nodes of degree 1: U=(u_1, ..., u_n)
        2) Remove node in U from the graph. However, not all are removed.
           If a node v is connected to multiple nodes in U, then only one of those is removed.
           This way, nodes that are connected to many nodes of degree 1 are being preserved.
           E.g. Edges (u_i, v), ..., (u_j, v) exist in the graph. Then for each v only one node u is removed.
        3) Repeat steps 1) and 2) (for a total number of 2)
        4) Find all nodes of degree 1 and remove ALL of them
        5) The remaining edges connect subgraphs of at least size 3

    :param graph:
    :param size: The minimum size the subgraphs should have (i.e. minimum number of nodes)
    :return: A list of edges that are save to remove
    '''
    graph_cp = copy.deepcopy(graph)

    # detect edges that connect subgraphs of size <= 3
    for _ in range(size-2):
        # get all nodes of degree 1
        node_degree1 = [node[0] for node in graph_cp.degree if node[1] == 1]

        # remove nodes of degree 1, but only one per adjacent node
        neighbours = []
        for node in node_degree1:
            if list(graph_cp.edges(node)) != []:  # it could be that it's neighbour was already removed
                neighbour = list(graph_cp.edges(node))[0][1]
                if neighbour not in neighbours:
                    neighbours.append(neighbour)
                    graph_cp.remove_node(node)

    # remove all nodes of degree 1. Edges that connect subgraphs of size >= 4 should be preserved this way
    node_degree1 = [node[0] for node in graph_cp.degree if node[1] == 1]
    graph_cp.remove_nodes_from(node_degree1)

    # the remaining edges all connect subgraphs of size 4 or more
    return list(graph_cp.edges(data=True))
    
    
def compute_viable_edges_cp2(graph: Graph, size: int = 4):
    ''' The method detects all edges that connect subgraphs of size 4 or more
    Procedure:
        Similar to compute_viable_edges_cp.
        However, it doesn't only look at immediate neighbors, but also neighbours within a certain range.
        Still doesn't work.

    :param graph:
    :param size: The minimum size the subgraphs should have (i.e. minimum number of nodes)
    :return: A list of edges that are save to remove
    '''
    graph_cp = copy.deepcopy(graph)
    # S = [graph_cp.subgraph(c).copy() for c in connected_components(graph_cp)]
    # for g in S:
    #     print(len(g.nodes))

    # detect edges that connect subgraphs of size <= 3
    steps_tot = size
    for i in range(size-2):
        steps_tot -= 1
        # get all nodes of degree 1
        node_degree1 = [node[0] for node in graph_cp.degree if node[1] == 1]

        # remove nodes of degree 1, but only one per adjacent node
        neighbours = []
        nodes_to_remove = []
        for node in node_degree1:
            if list(graph_cp.edges(node)) != []:  # it could be that it's neighbour was already removed
                neighbour = list(graph_cp.edges(node))[0][1]
                step = 1
                neighbours_visited = [node, neighbour]
                while (graph_cp.degree[neighbour] == 2) and (step < steps_tot):
                    neighbour_new = list(graph_cp.edges(neighbour))[0][1]
                    if neighbour_new in neighbours_visited:
                        neighbour = list(graph_cp.edges(neighbour))[1][1]
                    else:
                        neighbour = neighbour_new
                    step += 1
                    neighbours_visited.append(neighbour)

                if neighbour not in neighbours:
                    if step < steps_tot:
                        print('urghhh')
                        neighbours.append(neighbour)
                    nodes_to_remove.append(node)
        graph_cp.remove_nodes_from(nodes_to_remove)

    # remove all nodes of degree 1. Edges that connect subgraphs of size >= 4 should be preserved this way
    node_degree1 = [node[0] for node in graph_cp.degree if node[1] == 1]
    graph_cp.remove_nodes_from(node_degree1)

    # the remaining edges all connect subgraphs of size 4 or more
    return list(graph_cp.edges(data=True))

"""