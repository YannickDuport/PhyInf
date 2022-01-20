import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from pathlib import PosixPath


# currently not used
def create_parameter_file(tree_path: PosixPath, out_path: PosixPath):
    """ Reads all files containing model parameters (log_0.txt) from a given path.
    Parameters are stored in a txt-file in the format:
    <Tree-ID>,<MODEL>,<STATE_FREQUENCY_A>,<STATE_FREQUENCY_A>,<STATE_FREQUENCY_A>,<STATE_FREQUENCY_A>,
    <RATE_MATRIX_VALUE_A_TO_C>,<RATE_MATRIX_VALUE_A_TO_G>,<RATE_MATRIX_VALUE_A_TO_T>,
    <RATE_MATRIX_VALUE_C_TO_G>,<RATE_MATRIX_VALUE_C_TO_T>,<RATE_MATRIX_VALUE_G_TO_T>

    :param tree_path: Path containing model/tree directories
    :param out_path: Out path for final text file
    :return:
    """

    out_strings = []
    out_file = out_path / "model_parameters.txt"

    # iterate over all log-files
    for logfile in tree_path.rglob('log_0.txt'):
        tree_id = logfile.parent.name
        with logfile.open('r') as f:
            for line in f:
                if line.startswith("rate A <-> C"):
                    rate_A_C = float(line.split(':')[1])
                elif line.startswith("rate A <-> G"):
                    rate_A_G = float(line.split(':')[1])
                elif line.startswith("rate A <-> T"):
                    rate_A_T = float(line.split(':')[1])
                elif line.startswith("rate C <-> G"):
                    rate_C_G = float(line.split(':')[1])
                elif line.startswith("rate C <-> T"):
                    rate_C_T = float(line.split(':')[1])
                elif line.startswith("rate G <-> T"):
                    rate_G_T = float(line.split(':')[1])
                elif line.startswith("freq pi(A)"):
                    pi_A = float(line.split(':')[1])
                elif line.startswith("freq pi(C)"):
                    pi_C = float(line.split(':')[1])
                elif line.startswith("freq pi(G)"):
                    pi_G = float(line.split(':')[1])
                elif line.startswith("freq pi(T)"):
                    pi_T = float(line.split(':')[1])
        out_strings.append(
            f"{tree_id};GTR;{rate_A_C};{rate_A_G};{rate_A_T};{rate_C_G};{rate_C_T};{rate_G_T};{pi_A};{pi_C};{pi_G};{pi_T}\n"
        )

        with out_file.open('w') as out:
            for line in out_strings:
                out.write(line)


def compute_jc_distance(seq1, seq2):
    # compute hamming distance between seq1 and seq2
    hamming_dist = sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2)) / len(seq1)

    # compute jc corrected distance
    if hamming_dist >= 0.75:
      return np.inf
    else:
      return -3/4 * np.log(1 - 4/3 * hamming_dist)

def draw_network(graph):
    fig = plt.figure()
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos)
    nx.draw_networkx_labels(graph, pos)
    nx.draw_networkx_edge_labels(graph, pos)
    fig.show()

