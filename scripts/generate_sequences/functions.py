import pandas as pd
from pathlib import Path, PosixPath
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




