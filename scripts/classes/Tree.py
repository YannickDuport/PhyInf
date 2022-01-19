import numpy as np
from copy import deepcopy
import re
from scipy.linalg import expm
import math

from scripts.classes.Vertex import Vertex

class Tree:
    def __init__(self, name, rate_matrix=None, pi_rho=None, sequence_length=None):
        self.name = name
        self.nodes = {}  # a dictionary: keys will be names and values the Vertex objects
        self.edges = []  # list with edges: each element is a tuple (parent_name, child_name)
        self.leaves = []  # name of leaves
        self.root = ""  # name of root
        self.Newick = ""  # newick format
        self.verticesForPostOrderTraversal = []  # list of the names of the vertices for a post order traversal (leaves excluded)
        self.verticesForPreOrderTraversal = []  # list of the names of the vertices for a pre order traversal
        self.Q = rate_matrix  # rate matrix of the model
        self.pi = pi_rho  # root base compisition /state frequency
        self.map = {"A": 0, "C": 1, "G": 2, "T": 3}  # a map to remember the index of each base
        self.bases = list(self.map.keys())  # the bases
        self.ParsimonyScore = 0
        self.LogLikelihood = 0
        self.n = sequence_length  # length of the sequences_mini to be evolved

    def AddVertex(self, vertex):
        if vertex in self.nodes:
            print("Node ", vertex, " is already in the tree")
        else:
            self.nodes[vertex] = Vertex(vertex)

    def AddEdge(self, parent, child, distance, undirected=True):
        # modify the attributes of the nodes
        self.nodes[parent].children.append(child)
        self.nodes[child].parent = parent
        self.nodes[child].distanceToParent = distance
        self.edges.append((parent, child))
        # adjust the degrees
        self.nodes[child].indegree += 1
        self.nodes[parent].outdegree += 1
        self.nodes[parent].degree = self.nodes[parent].indegree + self.nodes[parent].outdegree
        self.nodes[child].degree = self.nodes[child].indegree + self.nodes[child].outdegree

        ##############NEW#####################
        if undirected:
            # create the undirected edges for the unrooted tree later:
            self.nodes[parent].Neighbors[child] = distance
            self.nodes[child].Neighbors[parent] = distance

    def SetLeaves(self):
        # fills the leaves list
        self.leaves = []
        for node in self.nodes:
            if self.nodes[node].outdegree == 0:
                self.leaves.append(node)

    def SetVerticesForPostOrderTraversal(self, include_Leaves=False):
        self.CleanTree()
        # computes a post order traversal
        if len(self.leaves) == 0:
            self.SetLeaves()

        self.verticesForPostOrderTraversal = []
        verticesToVisit = deepcopy(self.leaves)
        while len(verticesToVisit) > 0:
            c = verticesToVisit.pop()
            p = self.nodes[c].parent
            self.nodes[p].timesVisited += 1
            if self.nodes[p].timesVisited == 2:
                self.verticesForPostOrderTraversal.append(p)
                if self.nodes[p].indegree == 1:
                    verticesToVisit.append(p)
        if include_Leaves:  # include the leaves in the list
            self.verticesForPostOrderTraversal = self.leaves + self.verticesForPostOrderTraversal

    def SetVerticesForPreOrderTraversal(self):
        self.SetVerticesForPostOrderTraversal(include_Leaves=True)
        num_vertices = len(self.nodes)
        for i in range(num_vertices):
            self.verticesForPreOrderTraversal.append(self.verticesForPostOrderTraversal[num_vertices - 1 - i])

    def SetRoot(self):
        ### sets the root
        for node in self.nodes:
            if self.nodes[node].indegree == 0:
                self.root = node
                break

    def CleanTree(self):
        for n in self.nodes:
            self.nodes[n].Set = []
            self.nodes[n].timesVisited = 0

    def ComputeNewickFormat(self):
        # first initialize the leaves
        self.SetLeaves()
        for leaf in self.leaves:
            self.nodes[leaf].newickLabel = leaf + ":" + str(self.nodes[leaf].distanceToParent)

        # set the vertices for the post order traversal
        self.SetVerticesForPostOrderTraversal()
        for v in self.verticesForPostOrderTraversal:
            if self.nodes[v].indegree != 0:  # if not root
                c1, c2 = self.nodes[v].children[0], self.nodes[v].children[1]
                # generate the newick label for the vertex by combining the newick label of its children
                self.nodes[v].newickLabel = "(" + self.nodes[c1].newickLabel + "," + self.nodes[
                    c2].newickLabel + "):" + str(self.nodes[v].distanceToParent)
            else:  # if root
                c1, c2 = self.nodes[v].children[0], self.nodes[v].children[1]
                # if the root is reached, put a semicolon at the end
                self.nodes[v].newickLabel = "(" + self.nodes[c1].newickLabel + "," + self.nodes[c2].newickLabel + ");"

        self.SetRoot()
        self.Newick = self.nodes[self.root].newickLabel
        return self.Newick

    def ReadNewickFile(self, newickFile="", path=""):
        ###reads in a newick file
        if path != "":
            file = open(path, "r")
            newick_string = file.readline()
        else:
            newick_string = newickFile
        rx = r'\([^()]+\)'
        numTries = 0
        index = 1
        while "(" in newick_string:
            numTries += 1
            match = re.search(rx, newick_string)  # match the string to the regular expression
            string_match = match.group()  # get the match as a string
            siblings = string_match[1:-1]  # remove the two brackets ( and ) --> get the siblings
            if len(siblings.split(',')) == 3:
                # In some trees the root has three children. Group two together and add a new parent with distance 0
                to_replace = ','.join(siblings.split(',')[:2])
                newick_string = newick_string.replace(to_replace, f"({to_replace}):0")
                continue
            left, right = siblings.split(",")  # split into the left and right half
            left_name, left_distance = left.split(
                ":")  # split into the name of the left child and its distance to the parent
            right_name, right_distance = right.split(
                ":")  # split into the name of the right child and its distance to the parent
            if left_name not in self.nodes:
                self.AddVertex(left_name)
            if right_name not in self.nodes:
                self.AddVertex(right_name)
            hidden_v = "h" + str(index)  # hidden vertex is created
            index += 1
            self.AddVertex(hidden_v)
            # add adges to the new hidden vertex
            self.AddEdge(hidden_v, left_name, float(left_distance))
            self.AddEdge(hidden_v, right_name, float(right_distance))
            newick_string = newick_string.replace(string_match,
                                                  hidden_v)  # replace the leaf pair with the hidden vertex name

    def EvolveSequences(self):
        # evolve a sequence along the root
        self.SetVerticesForPreOrderTraversal()

        for node in self.verticesForPreOrderTraversal:
            if self.nodes[node].indegree == 0:  # if root
                self.nodes[node].sequence = "".join(np.random.choice(a=self.bases, size=self.n, p=self.pi))
            else:
                t = self.nodes[node].distanceToParent
                P = expm(self.Q * t)
                parent = self.nodes[node].parent
                for base in self.nodes[parent].sequence:
                    distr = P[self.map[base]]  # get the right row as prrob distribution
                    self.nodes[node].sequence += np.random.choice(a=self.bases, p=distr)  # sample new character

    def GetMRCA(self, node1, node2):
        # 4 conditions:
        # 1. node1 or node2 is root --> MRCA: root
        if node1 == self.root or node2 == self.root:
            mrca = self.root
        # 2. node is a child of node2 --> MRCA: node2
        elif node1 in self.nodes[node2].children:
            mrca = node2
        # 3. node2 is a child of node1 --> MRCA: node1
        elif node2 in self.nodes[node1].children:
            mrca = node1
        else:
            ancestors1 = [node1]  # list all the ancestors of node 1
            current = node1
            while current != self.root:  # traverse the parents untill the root is reached
                p = self.nodes[current].parent
                ancestors1.append(p)
                current = p
            current = node2  # find an ancestor of node 2 that is in ancestors1
            while current != self.root:
                p = self.nodes[current].parent
                if p in ancestors1:
                    mrca = p
                    break
                else:
                    current = p
        return (mrca)

    def ComputeTreeDistance(self, node1, node2):
        # find the path lenght between u and v
        mrca = self.GetMRCA(node1, node2)
        # wpl : weighted path length
        # if either of the vertex is the parent of the other,
        # you get the edge lenght between them
        if (node1 == self.nodes[node2].parent):
            wpl = self.nodes[node2].distanceToParent
        elif (node2 == self.nodes[node1].parent):
            wpl = self.nodes[node1].distanceToParent
        # otherwise, iterate over the parents of the node until you
        # reach to the mrca of the nodes
        else:
            wpl_to_mrca1 = 0
            wpl_to_mrca2 = 0
            # note that current changes after each iteration
            current = node1
            while (current != mrca):
                wpl_to_mrca1 += self.nodes[current].distanceToParent
                current = self.nodes[current].parent
            current = node2
            while (current != mrca):
                wpl_to_mrca2 += self.nodes[current].distanceToParent
                current = self.nodes[current].parent
            # sum the distances to the mrca to get the path length
            wpl = wpl_to_mrca1 + wpl_to_mrca2
        return (wpl)

    def FitchHartigan(self, site):
        """
        performs the fitch hartigan algorithm for one site of the sequences_mini
        Paramters:
        ----------
        site int
          The site of the sequence

        Returns:
        --------
        parsimony_score int
          The compute parsimony score for the site
        """
        # first initialize the root base sets
        for leaf in self.leaves:
            self.nodes[leaf].Set = [self.nodes[leaf].sequence[site]]
        score = 0
        for node in self.verticesForPostOrderTraversal:
            c1, c2 = self.nodes[node].children
            # set the set of the parent note
            intersection = list(set(self.nodes[c1].Set) & set(self.nodes[c2].Set))  # get intersection of children sets
            if len(intersection) == 0:
                score += 1
                self.nodes[node].Set = list(
                    set(self.nodes[c1].Set + self.nodes[c2].Set))  # set the current set as union of the children sets
            else:
                self.nodes[node].Set = intersection  # set the current set as intersection of the children sets

        return score

    def ComputeParsimonyScore(self):
        """
        Performs the Fitch Hartigan algorithm for each site and computes the total parsimony score for the whole tree
        Returns:
        --------
        parsimony_socres list of int
          The parsimony score for each site. To get the total sum of the tree, just sum the list up
        """
        self.CleanTree()
        self.SetLeaves()
        self.SetRoot()
        self.SetVerticesForPostOrderTraversal()
        n = len(self.nodes[self.leaves[0]].sequence)  # length of leaf sequences_mini
        parsimony_scores = []
        for i in range(n):  # go through each site
            # do fitch for the site
            current_score = self.FitchHartigan(site=i)
            parsimony_scores.append(current_score)
        self.ParsimonyScore = np.sum(parsimony_scores)
        return parsimony_scores

    def ComputeLikelihood(self):
        self.SetLeaves()
        self.SetRoot()
        # clean up old data in case there is any
        for n in self.nodes:
            self.nodes[n].LikelihoodVector = []
        # initialize the leave vectors

        for leaf in self.leaves:
            for site in self.nodes[leaf].sequence:
                vec = {'A': 0, 'C': 0, 'G': 0, 'T': 0}  # conditional likelihood vector to be filled
                vec[site] = 1
                self.nodes[leaf].LikelihoodVector.append(vec)

        self.SetVerticesForPostOrderTraversal()

        # perform a post-order traversal to compute the likelikood vectors for each node
        for node in self.verticesForPostOrderTraversal:
            c1, c2 = self.nodes[node].children  # get the chidren of the current node
            n = len(self.nodes[c1].LikelihoodVector)
            for i in range(n):
                # get the transition matrices to the children nodes
                P1 = expm(self.Q * self.nodes[c1].distanceToParent)
                P2 = expm(self.Q * self.nodes[c2].distanceToParent)
                vec = {}
                vec1, vec2 = self.nodes[c1].LikelihoodVector[i], self.nodes[c2].LikelihoodVector[i]
                for x in self.bases:
                    sum1 = [P1[self.map[x]][self.map[y]] * vec1[y] for y in self.bases]
                    sum2 = [P2[self.map[x]][self.map[y]] * vec2[y] for y in self.bases]
                    vec[x] = np.sum(sum1) * np.sum(sum2)
                self.nodes[node].LikelihoodVector.append(vec)

                # compute the log likelihood
        loglike = 0
        # iterate through each site
        for i in range(len(self.nodes[self.root].LikelihoodVector)):
            l_i = 0
            for b in self.bases:
                l_i += self.pi[self.map[b]] * self.nodes[self.root].LikelihoodVector[i][b]
            loglike += np.log(l_i)
        self.LogLikelihood = loglike
        return loglike

    def DeleteDirections(self):
        """
        Deletes the directed structure of the tree
        """
        for node in self.nodes:
            self.nodes[node].parent = ""
            self.nodes[node].children = []
            self.nodes[node].distanceToParent = 0
            self.nodes[node].indegree = 0
            self.nodes[node].outdegree = 0
        self.edges = []  # delete directed edges

    def SurpressRoot(self):
        """
        Supresses the root and turns the tree into an unrooted, undirected tree
        """
        self.SetRoot()
        # the structures for the unrooted tree are already in place --> delete directed structure
        new_dist = 0
        for c in self.nodes[self.root].children:
            del self.nodes[c].Neighbors[self.root]  # delete the root from the neighbors
            new_dist += self.nodes[c].distanceToParent  # get the new distance between the two nodes

        c1, c2 = self.nodes[self.root].children
        self.nodes[c1].Neighbors[c2] = new_dist
        self.nodes[c2].Neighbors[c1] = new_dist

        del self.nodes[self.root]  # remove the root from the set of nodes
        self.DeleteDirections()  # delete all the directions from the tree
        self.root = ""

    def RestoreDirections(self, new_root):
        """
        Creates a directed tree from the given undircted structure
        """
        verticesToChange = [new_root]

        while len(verticesToChange) > 0:
            current = verticesToChange.pop()
            for n in self.nodes[current].Neighbors:
                if n != self.nodes[current].parent:  # if not parent, add as child
                    d = self.nodes[current].Neighbors[n]  # get the distance to the neighbor
                    self.AddEdge(parent=current, child=n, distance=d, undirected=False)  # only add a directed edge
                    verticesToChange.append(n)

    def Reroot(self, v1, v2, root_name):
        """
        Reroots the tree on the edge between v1 and v2 --> restores the directed structure
        """
        self.AddVertex(root_name)  # add the new root to the structure
        dist = self.nodes[v1].Neighbors[v2]  # get the distance between the two nodes
        del self.nodes[v1].Neighbors[v2]
        del self.nodes[v2].Neighbors[v1]  # delete the edge between v1 and v2
        self.nodes[root_name].Neighbors = {v1: dist, v2: 0}  # add them to the neighbors witht the correct edge weight
        self.nodes[v1].Neighbors[root_name] = dist
        self.nodes[v2].Neighbors[root_name] = 0
        self.RestoreDirections(new_root=root_name)
        self.root = root_name

    def PerformNNI(self, u, v):
        """
        Performs 2 NNI moves and returns the new trees
        """
        if self.root != "":
            self.SurpressRoot()  # if the root hasn't bene supressed yet, do that here

        # create two new tree objects --> these are two unrooted, undirected trees
        # always make sure that u and v each have 3 neighbors, otherwise you can't interchange them
        # so the edge needs to be an inner edge
        if self.nodes[u].degree != 3 or self.nodes[v].degree != 3:
            print("Please choose an inner edge!")
            return
        t1 = deepcopy(self)
        t2 = deepcopy(self)
        NN_u = deepcopy(self.nodes[u].Neighbors)
        del NN_u[v]  # two neighbors to switch: u1, u2
        u1, u2 = list(NN_u.keys())
        NN_v = deepcopy(self.nodes[v].Neighbors)
        del NN_v[u]  # two neighbors to switch: v1, v2
        v1, v2 = list(NN_v.keys())

        # change 1: one side has u1 & v1 (let's say the u side) and the other u2 & v2 (the v side)
        # u1 is already in there, so we need to add v1 and delete u2
        t1.nodes[u].Neighbors[v1] = NN_v[v1]  # keep the right distance
        del t1.nodes[u].Neighbors[u2]  # delete u2
        # v2 is already a neighbor of v, add u2 and delete v1
        t1.nodes[v].Neighbors[u2] = NN_u[u2]
        del t1.nodes[v].Neighbors[v1]

        # adjust the neighbor sets of the neighbors that got switched
        t1.nodes[v1].Neighbors[u] = NN_v[v1]
        del t1.nodes[v1].Neighbors[v]
        t1.nodes[u2].Neighbors[v] = NN_u[u2]
        del t1.nodes[u2].Neighbors[u]

        # change 2: one side has u1 & v2 and the other u2 & v1
        t2.nodes[u].Neighbors[v2] = NN_v[v2]  # keep the right distance
        del t2.nodes[u].Neighbors[u2]  # delete u2
        # v1 is already a neighbor of v, add u2 and delete v2
        t2.nodes[v].Neighbors[u2] = NN_u[u2]
        del t2.nodes[v].Neighbors[v2]

        # adjust the neighbor sets of the neighbors that got switched
        t2.nodes[v2].Neighbors[u] = NN_v[v2]
        del t2.nodes[v2].Neighbors[v]
        t2.nodes[u2].Neighbors[v] = NN_u[u2]
        del t2.nodes[u2].Neighbors[u]

        return t1, t2

    def GetLikelihoods(self, v1, v2, t, site):
        """
        Computes the site-wise likelihoods and their derivatives needed to compute the
        partial derivatives of the log likelihood
        Parameters:
        -----------
        v1 str
          left child of the root, has distance 0 to the root

        v2 str
          right child of the root, has the original edge length as distance

        t float
          the current distance to (gets updated during the NW iterations)

        site int
          the site we are currently on

        Returns:
        --------
        L_i float
          The sitewise likelihood

        L_i_first float
          The first derivative of L_i

        L_i_second float
          The second derivative of L_i
        """
        L_i = 0
        L_i_first = 0  # first derivative of L_i
        L_i_second = 0  # second derivative of L_i
        tm = expm(self.Q * t)  # get the transition matrix: e^(Qt)
        tm_1 = self.Q.dot(tm)  # Q*e^(Qt)
        tm_2 = self.Q.dot(tm_1)  # Q*Q*e^(Qt)
        for x in self.bases:
            # get the likelihood vectors
            vector_v1 = self.nodes[v1].LikelihoodVector[site]
            vector_v2 = self.nodes[v2].LikelihoodVector[site]

            # left side of the sum --> the same for all three formuals
            # use the child with distance 0 --> v2
            left_side = self.pi[self.map[x]] * vector_v2[x]

            # compute the right side of the sum: sum_z matrix(z|x) * L^(i)_(v) (z)
            # the right side for the derivatives is nearly the same, just the matrix changes
            # attention: use the child whose distance to the root is zero for the left side and the other for the right side
            # here: v2 has distance 0 --> use v1 for the right side
            right_side = [tm[self.map[x]][self.map[z]] * vector_v1[z] for z in self.bases]  # right side for L_i
            right_side1 = [tm_1[self.map[x]][self.map[z]] * vector_v1[z] for z in
                           self.bases]  # right side for the first derivative of L_i
            right_side2 = [tm_2[self.map[x]][self.map[z]] * vector_v1[z] for z in
                           self.bases]  # right side for the second derivative of L_i
            # add to the likelihoods
            L_i += left_side * np.sum(right_side)
            L_i_first += left_side * np.sum(right_side1)
            L_i_second += left_side * np.sum(right_side2)

        return L_i, L_i_first, L_i_second

    def GetDerivatives(self, v1, v2, t):
        """
        Computes the first and second derivatives of the log-likelihood function

        Parameters:
        -----------
        v1 str
          left child of the root, has distance 0 to the root

        v2 str
          right child of the root, has the original edge length as distance

        t float
          the current distance to (gets updated during the NW iterations)

        Returns:
        -------
        l_first float
          The first order partial derivative of the log-likelihood score

        l_second float
          The second order partial derivative of the log-likelihood score
        """
        l_first = 0  # first partial derivative
        l_second = 0  # second partial derivative
        for i in range(self.n):
            L_i, L_i_first, L_i_second = self.GetLikelihoods(v1=v1, v2=v2, t=t,
                                                             site=i)  # get the sitewise lieklihoods and its derivatives
            l_first += (1 / L_i) * L_i_first
            l_second += ((1 / L_i) * L_i_second) - math.pow(((1 / L_i) * L_i_first), 2)

        return l_first, l_second

    def NewtonRaphson(self, v1, v2, threshold, max_iter):
        """
        Optimizes one edge length using the Newton-Raphson method

        Parameters:
        -----------
        v1 str
          The parent node of the edge that needs to optimized

        v2 str
          The child node of the edge that needs to be optimized

        threshold float
          The threshold for convergence

        max_iter int
          The maximum number of iterations

        Returns:
        --------
        t_i float
          The optimized branch length for the edge
        """
        t_i = self.nodes[v1].distanceToParent  # v2 has distance 0
        diff = 0.01
        iteration = 0
        while np.abs(diff) > threshold and iteration < max_iter:
            print("Iteration: ", iteration + 1)
            l_first, l_second = self.GetDerivatives(v1=v1, v2=v2, t=t_i)
            diff = l_first / l_second
            t_i -= diff
            iteration += 1
        return t_i

    def OptimizeEdgeLengths(self, threshold=0.000001, max_iter=10):
        """
        Optimizes the branch lengths of the tree using the Newton-Raphson method

        Paramters:
        ----------
        threshold float
          The threshold for convergence

        max_iter int
          The maximum number of iterations

        Returns:
        --------
        optimized_edges dic
          Contains the optimized branch lengths. The keys are tuples (edges) and the values are the new distances
        """
        self.CleanTree()
        self.SetRoot()
        self.SetLeaves()

        NR_tree = deepcopy(self)
        edges = deepcopy(self.edges)
        c1, c2 = self.nodes[self.root].children
        edges.append((c1, c2))
        # iterate through the edges
        optimized_edges = {}
        index = 1
        for (u, v) in edges:
            if u == self.root or v == self.root:  # skip the root
                continue
            print("Optimizing the edge between ", u, " and ", v)
            print("Original length ", self.nodes[v].distanceToParent)
            NR_tree.SurpressRoot()
            root_name = "root" + str(index)
            NR_tree.Reroot(v1=u, v2=v, root_name=root_name)  # set a new root
            NR_tree.ComputeLikelihood()  # compute the likelihood vectors
            optimized_length = NR_tree.NewtonRaphson(v1=u, v2=v, threshold=threshold, max_iter=max_iter)
            optimized_edges[(u, v)] = optimized_length
            print("Optimized length ", optimized_length)

        return optimized_edges