from pathlib import Path
import time

ROOT_PATH = Path(__file__).parent.parent.parent
DATA_PATH = ROOT_PATH / "data"
TREE_PATH = DATA_PATH / "trees"
TREE_MINI_PATH = DATA_PATH / "trees_mini"
SEQ_MINI_PATH = DATA_PATH / 'sequences_mini'
SEQ_PATH = DATA_PATH / "sequences"

def timeit(method):
    """Timing decorator"""

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('{:20}  {:8.4f} [s]'.format(method.__name__, (te - ts)))
        return result

    return timed