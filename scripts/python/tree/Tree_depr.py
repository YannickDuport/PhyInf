import re
import copy
import warnings

import numpy as np
import matplotlib.pyplot as plt

from copy import deepcopy
from itertools import combinations
from scipy.linalg import expm
from warnings import warn
from queue import Queue
from ast import literal_eval

# Jukes-cantor rate matrix
rate_matrix_jc = np.array([[-1, 1/3, 1/3, 1/3],
                            [1/3, -1, 1/3, 1/3],
                            [1/3, 1/3, -1, 1/3],
                            [1/3, 1/3, 1/3, -1]])

stationary_distribution = [0.25, 0.25, 0.25, 0.25]

class Vertex:
    def __init__(self, name, sequence=None):
        self.name = name
        self.descendants = [] #list with name of descendants
        self.parent = "" #name of parent
        self.degree = 0
        self.indegree = 0
        self.outdegree = 0
        self.distance_to_parent = 0
        self.newick_label = ""
        self.times_visited = 0 #for post order traversal
        self.sequence = sequence
        self.parsimony_set = set()
        self.conditional_likelihood = None
        self.neighbors = {} #the neighbors in the unrooted tree, keys are node name, values are distance
        self.Set = []       # ??? NO IDEA WHAT'S THAT FOR ???


class Tree:
    def __init__(self, name, rate_matrix=None, pi=None):
        self.name = name
        self.nodes = {}  # a dictionary: keys will be names and values the Vertex objects
        self.edges = []  # list with edges: each element is a tuple (parent_name, child_name)
        self.leaves = []  # name of leaves
        self.root = ""  # name of root
        self.newick = ""  # newick string of the tree
        self.vertices_for_post_order_traversal = []  # list of the names of the vertices for a post order traversal (leaves excluded)
        self.vertices_for_pre_order_traversal = []  # list of the names of the vertices for a pre order traversal
        self.rate_matrix = rate_matrix  # rate matrix of the tree
        self.log_likelihood = None  # log-likelihood of the tree
        self.pi = pi  # relative base frequencies
        # self.map = {"A": 0, "C" : 1, "G" : 2, "T": 3}   #a map to remember the index of each base
        # self.bases = list(self.map.keys())              #the bases
        self.nuc_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}  # map from base to index
        self.int_to_nuc = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}  # map from index to base

    def add_vertex(self, vertex, sequence=None):
        """ Adds a Vertex to the tree """
        if vertex in self.nodes:
            warnings.warn(f"Node {vertex} is already in the tree. Choose a different name")
        else:
            self.nodes[vertex] = Vertex(vertex, sequence)

    def add_edge(self, parent, child, distance, undirected=True):
        """ Adds an edge between two nodes """

        # modify the attributes of the nodes
        self.nodes[parent].descendants.append(child)
        self.nodes[child].parent = parent
        self.nodes[child].distance_to_parent = distance
        self.edges.append((parent, child))

        # adjust the degrees
        self.nodes[child].indegree += 1
        self.nodes[parent].outdegree += 1
        self.nodes[parent].degree = self.nodes[parent].indegree + self.nodes[parent].outdegree
        self.nodes[child].degree = self.nodes[child].indegree + self.nodes[child].outdegree

        ##############NEW#####################
        if undirected:
            # create the undirected edges for the unrooted tree later:
            self.nodes[parent].neighbors[child] = distance
            self.nodes[child].neighbors[parent] = distance


    def set_vertices_for_post_order_traversal(self, include_Leaves=False):
        """ Creates an array containing all vertices in the order of post order traversal """
        self.clean_tree()

        # set leaves if not done yet
        if len(self.leaves) == 0:
            self.set_leaves()

        # compute post order traversal
        self.vertices_for_post_order_traversal = []
        vertices_to_visit = deepcopy(self.leaves)
        while len(vertices_to_visit) > 0:
            c = vertices_to_visit.pop()
            p = self.nodes[c].parent
            self.nodes[p].times_visited += 1
            if self.nodes[p].times_visited == 2:
                self.vertices_for_post_order_traversal.append(p)
                if self.nodes[p].indegree == 1:
                    vertices_to_visit.append(p)

        # include the leaves in the list
        if include_Leaves:
            self.vertices_for_post_order_traversal = self.leaves + self.vertices_for_post_order_traversal

    def infer_tree(self, sequences: dict):  # sequences is a dictionary: keys - taxon names, values - sequence strings

        # check number of sequences
        #

        # Clear tree (delete all nodes, edges, etc)
        self.clear_tree()

        # Assign a rate matrix and it's stationary distribution to the tree, if it doesn't have one yet
        if self.rate_matrix is None:
            self.rate_matrix = rate_matrix_jc  # FIXME: How to compute a good guess for a rate matrix?
        if self.pi is None:
            self.computes_stationary_distribution()

        # shuffle sequences
        taxa = list(sequences.keys())
        np.random.shuffle(taxa)

        # INITIALIZE
        # get two nodes and their sequences
        n1 = taxa.pop()
        n2 = taxa.pop()
        s1 = sequences[n1]
        s2 = sequences[n2]

        # Create a two-node tree and optimize branch lenght between them
        self.__initialize_tree((n1, n2), (s1, s2))

        # successively add more sequences
        for i in range(len(taxa)):
            print(f"\n###################################\nADD TAXA {i+1}\n#####################################")
            # add new taxa
            n = taxa.pop()
            s = sequences[n]
            self.add_node(n, s, i)

            # compute ancestral sequence for the newly introduced hidden node
            h = f"h{i}"
            #self.compute_ancestral_sequence(h)
            #self.reroot(h, n, 'root')


    def __initialize_tree(self, nodes, sequences):

        n1, n2 = nodes
        s1, s2 = sequences

        # compute jc-distance between sequences, as an initial guess for edge length
        d12 = self.__compute_jc_distance(s1, s2)

        # add nodes, a root and edges to the tree
        self.add_vertex(n1, s1)
        self.add_vertex(n2, s2)
        self.add_vertex('root')
        self.add_edge('root', n1, 0, undirected=True)
        self.add_edge('root', n2, d12, undirected=True)

        # optimize the edge length
        self.set_root()
        self.optimize_branch_length(n1, n2, d12)

    def add_node(self, n, s, i, threshold=1e-3):

        root = self.root
        self.supress_root()
        edges = self.get_undirected_edges()

        # array to store results [(<edge>, <x,y,z>, <likelihood>)]
        results = []
        # add node to each edge in the tree
        for edge in edges:
            tree = copy.deepcopy(self)
            (n1, n2), d12 = edge
            s1 = tree.nodes[n1].sequence
            s2 = tree.nodes[n2].sequence

            # compute initial guesses for edge weights
            # x is the distance from the new node n to the newly introduced hidden node
            # y is the distance from node n1 to the newly introduced hidden node
            d13 = tree.__compute_jc_distance(s1, s)
            d23 = tree.__compute_jc_distance(s2, s)
            d12 = tree.__compute_jc_distance(s1, s2)
            x = (d13 + d23 - d12) / 2
            y = d13 - x
            z = d23 - x
            if z < 0 or y < 0:
                print("z or y is negative!!!!!!!!!")

            x_init, y_init, z_init = x, y, z
            #print(x, y, z)

            # add hidden node on the edge
            h = f"h{i}"
            tree.add_vertex_on_edge(h, (n1, n2), w=(y, z))

            # add new node n to the tree
            tree.add_vertex(n, s)
            tree.add_edge(n, h, distance=x, undirected=True)

            # compute initial likelihood after adding the new node
            tree.reroot(n, h, 'root')
            likelihood = tree.compute_likelihood()

            # optimize branch length x, y, z
            print(f"\nOptimize x")
            x = tree.optimize_branch_length(n1=n, n2=h, t=x, threshold=1e-5)
            print(f"\nOptimize y")
            y = tree.optimize_branch_length(n1=n1, n2=h, t=y, threshold=1e-5)
            print(f"\nOptimize z")
            z = tree.optimize_branch_length(n1=n2, n2=h, t=z, threshold=1e-5)
            likelihood_new = tree.log_likelihood

            # do a 2nd round of branch length optimization
            print(f"\nOptimize x")
            x = tree.optimize_branch_length(n1=n, n2=h, threshold=1e-5)
            print(f"\nOptimize y")
            y = tree.optimize_branch_length(n1=n1, n2=h, threshold=1e-5)
            print(f"\nOptimize z")
            tree.optimize_branch_length(n1=n2, n2=h, threshold=1e-5)

            tree.supress_root()
            x = tree.nodes[n].neighbors[h]
            y = tree.nodes[n1].neighbors[h]
            z = tree.nodes[n2].neighbors[h]
            print(f"initial (x, y, z): {(x_init, y_init, z_init)}")
            print(f"final (x, y, z): {(x, y, z)}")
            print(f"y+z (init): {self.nodes[n1].neighbors[n2]}; y+z (final): {y+z}")

            results.append([tree.log_likelihood, (n1, n2), (x, y, z)])

        optimal_topology = max(results, key=lambda x: x[0])
        # add hidden node and new node to actual tree
        _, (u, v), (x, y, z) = optimal_topology
        self.add_vertex_on_edge(h, (u, v), (y, z))
        self.nodes[h].neighbors[v] = z
        self.add_vertex(n, s)
        self.add_edge(h, n, y, undirected=True)
        self.compute_ancestral_sequence(h, u, v)
        self.reroot(h, n, 'root')

    ######   LIKELIHOOD COMPUTATION       #####################
    ######   BRANCH LENGTH OPTIMIZATION   #####################

    def compute_likelihood(self):
        """ Computes the log-likelihood of the complete tree """

        # set leaves, root, and the array containing vertices for post order traversal
        # FIXME: Can probably be removed
        self.set_leaves()
        self.set_root()
        self.set_vertices_for_post_order_traversal()

        # check whether all leaves have sequences assigned, and whether sequences have same length
        seq_len = []
        for leaf in self.leaves:
            if self.nodes[leaf].sequence is None:
                raise AttributeError(f"Leaf {leaf} has no sequence assigned")
            seq_len.append(len(self.nodes[leaf].sequence))
        if len(np.unique(seq_len)) != 1:
            raise ValueError(f"Leaf sequences have not the same lengths: {seq_len}")

        # set the conditional likelihoods in each node to 0
        for node in self.nodes.values():
            node.conditional_likelihood = np.zeros(shape=(seq_len[0], 4))

        # compute log_likelihood by iterating over the sites, computing the site's likelihood and adding up all site likelihoods
        log_likelihood = 0
        for site in np.arange(seq_len[0]):
            likelihood = self.compute_likelihood_site(site)
            log_likelihood += np.log(self.compute_likelihood_site(site))

        # assign the log-likelihood to the tree's attribute and return it
        self.log_likelihood = log_likelihood
        return log_likelihood

    def compute_likelihood_site(self, site: int):
        """ computes the likelihood of the tree for one position in the sequence """

        # set the conditional likelihoods in the leaves
        for leaf in self.leaves:
            leaf = self.nodes[leaf]
            state = leaf.sequence[site]
            leaf.conditional_likelihood[site] = [1 if state == nuc else 0 for nuc in ['A', 'C', 'G', 'T']]

        # compute the conditional likelihoods of each hidden node by going from the nodes towards the rood (post order traversal)
        for node in self.vertices_for_post_order_traversal:
            node = self.nodes[node]
            conditional_likelihood = np.array([1, 1, 1, 1])

            for child in node.descendants:
                child = self.nodes[child]
                t = float(child.distance_to_parent)
                p = expm(self.rate_matrix * t)
                l = np.dot(p, child.conditional_likelihood[site])
                conditional_likelihood = conditional_likelihood * l
            node.conditional_likelihood[site] = conditional_likelihood

        # final site likelihood is the dot-product of the stationary distribution and the conditional likelihood in the root
        root = self.nodes[self.root]
        likelihood = np.dot(self.pi, root.conditional_likelihood[site])

        return likelihood

    def optimize_branch_length(self, n1, n2, t=None, threshold=1e-5, max_iter=20):
        """ optimizes the edge length between two nodes n1 and n2.
            The method expects an initial guess t for the edge length """

        if t is None:
            t = self.nodes[n1].neighbors[n2]
        root = self.root

        # reroot the tree along the edge (n1, n2) and compute the initial log-likelihood
        self.supress_root()
        self.reroot(n1, n2, root)
        self.compute_likelihood()
        print(f"step: 0; t: {t}; log_likelihood: {self.log_likelihood}")

        # optimize edge length
        diff = np.inf
        step = 0
        while diff > threshold and step < max_iter:

            # compute new edge length
            t_next = self.compute_next_branch_length(n1, n2, root, t)

            # update edge length
            self.nodes[n1].neighbors[root] = t_next
            self.nodes[n1].distance_to_parent = t_next

            # update log likelihood and conditional likelihood (only the conditional likelihoods in the root need to be updated!)
            log_likelihood = 0
            for site in range(len(self.nodes[self.leaves[0]].sequence)):
                conditional_likelihood = np.array([1, 1, 1, 1])
                for n in [n1, n2]:
                    dist = float(self.nodes[n].distance_to_parent)
                    p = expm(self.rate_matrix * dist)
                    l = np.dot(p, self.nodes[n].conditional_likelihood[site])
                    conditional_likelihood = conditional_likelihood * l
                self.nodes[root].conditional_likelihood[site] = conditional_likelihood
                log_likelihood += np.log(np.dot(self.pi, conditional_likelihood))
            self.log_likelihood = log_likelihood

            diff = abs(t - t_next)
            t = t_next
            step += 1

            print(f"step: {step}; t: {t_next}; log_likelihood: {self.log_likelihood}")

        return t_next
        # # unroot tree so it can be rerooted again
        # self.supress_root()

        # self.reroot(n1, n2, root)

    def compute_next_branch_length(self, n1, n2, root, t):
        """Newton-Raphson to compute next branch length"""

        # get vertices
        n1 = self.nodes[n1]
        n2 = self.nodes[n2]
        root = self.nodes[root]

        # compute p1 (Q*e^(Qt)) and p2 (Q^2*e^(Qt))
        p_1 = np.dot(self.rate_matrix, expm(self.rate_matrix * t))
        p_2 = np.dot(np.dot(self.rate_matrix, self.rate_matrix), expm(self.rate_matrix * t))

        ### compute derivatives of log likelihoods ###

        derivative_1 = 0
        derivative_2 = 0
        # iterate over sites
        for i in range(len(self.nodes[self.leaves[0]].sequence)):
            derivative_site_likelihood_1 = 0
            derivative_site_likelihood_2 = 0
            # compute site likelihood
            site_likelihood = np.dot(self.pi, root.conditional_likelihood[i])
            # compute 1st and 2nd derivative of site likelihoods
            derivative_site_likelihood_1 += np.dot([self.pi * n2.conditional_likelihood[i]],
                                                   np.dot(p_1, n1.conditional_likelihood[i]))
            derivative_site_likelihood_2 += np.dot([self.pi * n2.conditional_likelihood[i]],
                                                   np.dot(p_2, n1.conditional_likelihood[i]))
            # update log likelihood derivatives
            derivative_1 += (1 / site_likelihood) * derivative_site_likelihood_1
            derivative_2 += (1 / site_likelihood) * derivative_site_likelihood_2 - \
                            np.power((1 / site_likelihood) * derivative_site_likelihood_1, 2)

        # compute t_next
        t_next = (t - derivative_1 / derivative_2).item(0)

        # make sure t doesn't diverge or is negative
        if t_next > 5:  # no idea what a reasonable upper bound for edge length is. 5 is probably way too large
            warn(
                f"edge length between {n1.name} and {n2.name} is diverging: t(i+1) = {t_next}. Edge length ({n1.name}, {n2.name}) is set to 1e-4")
            t_next = 1e-5  # from slides
        if t_next < 0:  # not sure what to do in case of negative edge lengths
            # raise ValueError(
            #     f"New edge length between {n1.name} and {n2.name} is negative: {t_next}. Use a better initial guess for edge ({n1.name}, {n2.name})")
            warn(f"New edge length between {n1.name} and {n2.name} is negative: {t_next}. New initial guess is 1e-4)")
            if t <= 1e-5:
                t_next = 1e-4
            else:
                t_next = t/10
        return t_next


    def compute_ancestral_sequence(self, n, u, v):
    # FIXME: Not sure whether this makes sense.
    # Idea is, that since all nodes already have sequences, only the immediate neighbors have to be checked
    # Then it should be ok to just take the 'consensus'
        if self.nodes[n].sequence is not None:
            raise ValueError(f"Node {n} already has a sequence assigned")

        #self.supress_root()
        neighbors = self.nodes[n].neighbors.keys()
        sequences = [self.nodes[u].sequence, self.nodes[v].sequence]
        new_seq = ""
        for site in range(len(sequences[0])):
            new_seq += {seq[site] for seq in sequences}.pop()
        self.nodes[n].sequence = new_seq

    def __compute_jc_distance(self, seq1: str, seq2: str):
        """ Computes the jukes-cantor distance between two sequences
            If the hamming distance is too large (i.e. >= 0.75) the functions returns infinity
        """
        # compute hamming distance between seq1 and seq2
        hamming_dist = sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2)) / len(seq1)
        # # compute jc corrected distance
        if hamming_dist >= 0.75:
            return np.inf
        else:
            return -3 / 4 * np.log(1 - 4 / 3 * hamming_dist)

    def computes_stationary_distribution(self):
        """ Computes the stationary distribution of the tree's rate matrix and assigns it to the variable 'pi'. """

        if self.rate_matrix is None:
            raise AttributeError("No rate matrix was assigned")

        rate_matrix = deepcopy(self.rate_matrix)
        # define the right hand side of the equation
        b = np.array([[0], [0], [0], [0], [1]])
        # transpose the rate matrix
        Q_t = rate_matrix.transpose()
        # add a new row for the distribution condition
        newrow = [1, 1, 1, 1]
        Q_t = np.vstack([Q_t, newrow])
        # solve the linear equation system using numpy (finds the least-squares solution)
        pi_stationary = np.linalg.lstsq(Q_t, b, rcond=None)[0]

        self.pi = pi_stationary.transpose()[0]

    def set_leaves(self):
        """ Finds all the leaves in the tree and assigns them to the leaves attribute"""
        self.leaves = []
        for node in self.nodes:
            if self.nodes[node].outdegree == 0:
                self.leaves.append(node)

    def set_root(self):
        """ Finds the root of a tree and assigns it to the attribute 'root' """
        for node in self.nodes:
            if self.nodes[node].indegree == 0:
                self.root = node
                break

    def clean_tree(self):
        """ Resets the attributes 'set' and 'times_visited' of each node in the tree """
        for n in self.nodes:
            self.nodes[n].Set = []
            self.nodes[n].times_visited = 0

    def clear_tree(self):
        """ Clears the tree, i.e. deletes all nodes, edges, etc.. Keeps rate matrix"""
        self.nodes = {}
        self.edges = []
        self.leaves = []
        self.root = ""
        self.newick = ""
        self.vertices_for_post_order_traversal = []
        self.vertices_for_pre_order_traversal = []
        self.log_likelihood = None

    def get_undirected_edges(self):
        undirected_edges = []
        for n in self.nodes:
            for neighbor, weight in self.nodes[n].neighbors.items():
                if ((neighbor, n), weight) not in undirected_edges:
                    undirected_edges.append(((n, neighbor), weight))
        return undirected_edges

        ######   UN- AND REROOTING TREE   #####################

    def delete_directions(self):
        """
        Deletes the directed structure of the tree
        """
        for node in self.nodes:
            self.nodes[node].parent = ""
            self.nodes[node].descendants = []
            self.nodes[node].distance_to_parent = 0
            self.nodes[node].indegree = 0
            self.nodes[node].outdegree = 0
        self.edges = []  # delete directed edges

    def supress_root(self):
        """
        Supresses the root and turns the tree into an unrooted, undirected tree
        """
        self.set_root()

        # the structures for the unrooted tree are already in place --> delete directed structure
        new_dist = 0
        for c in self.nodes[self.root].descendants:
            del self.nodes[c].neighbors[self.root]  # delete the root from the neighbors
            new_dist += self.nodes[c].distance_to_parent  # get the new distance between the two nodes

        c1, c2 = self.nodes[self.root].descendants
        self.nodes[c1].neighbors[c2] = new_dist
        self.nodes[c2].neighbors[c1] = new_dist

        del self.nodes[self.root]  # remove the root from the set of nodes
        self.delete_directions()  # delete all the directions from the tree
        self.root = ""

    def restore_directions(self, new_root):
        """
        Creates a directed tree from the given undircted structure
        """
        vertices_to_change = [new_root]

        while len(vertices_to_change) > 0:
            current = vertices_to_change.pop()
            for n in self.nodes[current].neighbors:
                if n != self.nodes[current].parent:  # if not parent, add as child
                    d = self.nodes[current].neighbors[n]  # get the distance to the neighbor
                    self.add_edge(parent=current, child=n, distance=d, undirected=False)  # only add a directed edge
                    vertices_to_change.append(n)

    def reroot(self, v1, v2, root_name):
        """
        Reroots the tree on the edge between v1 and v2 --> restores the directed structure
        """
        self.add_vertex(root_name)  # add the new root to the structure
        dist = self.nodes[v1].neighbors[v2]  # get the distance between the two nodes
        del self.nodes[v1].neighbors[v2]
        del self.nodes[v2].neighbors[v1]  # delete the edge between v1 and v2
        self.nodes[root_name].neighbors = {v1: dist, v2: 0}  # add them to the neighbors witht the correct edge weight
        self.nodes[v1].neighbors[root_name] = dist
        self.nodes[v2].neighbors[root_name] = 0
        self.restore_directions(new_root=root_name)
        self.root = root_name

    def add_vertex_on_edge(self, node, edge, w):
        """ Adds a new vertex onto a preexisting edge
        :param node: name of the new node
        :param edge: edge, tuple of nodes (v1, v2)
        :param w: distance from v1 to the new node
        """
        v1, v2 = edge
        w1, w2 = w
        dist = self.nodes[v1].neighbors[v2]
        # if w > dist:
        #     raise ValueError(
        #         f"The passed branch length is greater than the branch length between {v1} and {v2}: {w} > {dist}")

        self.add_vertex(node)  # add the new node to the structure
        dist = self.nodes[v1].neighbors[v2]  # get the distance between the two nodes
        del self.nodes[v1].neighbors[v2]
        del self.nodes[v2].neighbors[v1]  # delete the edge between v1 and v2
        self.nodes[node].neighbors = {v1: w1, v2: w2}  # add them to the neighbors witht the correct edge weight
        self.nodes[v1].neighbors[node] = w1
        self.nodes[v2].neighbors[node] = w2
