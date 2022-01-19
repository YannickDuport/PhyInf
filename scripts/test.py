from scripts.classes.Tree import Tree
from scripts import TREE_MINI_PATH, TREE_PATH
from scripts.generate_sequences.functions import create_parameter_file

path = TREE_MINI_PATH
trees = []
for file in path.rglob("*.newick"):
    tree_id = file.parent.name
    tree = Tree(tree_id)
    n = tree.ReadNewickFile(path=file)
    print(tree_id, tree.ComputeNewickFormat())

create_parameter_file(TREE_MINI_PATH, TREE_MINI_PATH)
