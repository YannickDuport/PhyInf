class Vertex:
  def __init__(self,name):
    self.name = name
    self.children = [] #list with name of children
    self.parent = "" #name of parent
    self.degree = 0
    self.indegree = 0
    self.outdegree = 0
    self.distanceToParent = 0
    self.newickLabel = ""
    self.timesVisited = 0 #for post order traversal
    self.sequence = ""
    self.Set = []
    self.LikelihoodVector = [] #one dictionary per site
    self.Neighbors = {} #the neighbors in the unrooted tree, keys are node name, values are distance
