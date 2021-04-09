import os
import pickle

import sys
sys.path.append('../python_data/')
from ed_result_wrapper import EDResult


def get_result(param, datadir, dump_dir='results_dump/'):
    """ This either reads from a result from the HDF5 file, or a pickle dump
        from a previous run. If new results from a new run should be loade, the
        appropriate pickle dump must be cleared, otherwise it defaults to that.
        This is done for performance reasons.
    """
    tag = "_{:s}_{:d}x{:d}x{:d}_ws{:d}{:d}{:d}".format(param['gauge_particles'], *param['L'], *param['winding_sector'])
    result_file = f"{dump_dir}/result_dump{tag}.pickle"

    result = None
    if os.path.isfile(result_file):
        result = pickle.load(open(result_file, 'rb'))
    else:
        result = EDResult(param, datadir=datadir)
        pickle.dump(result, open(result_file, 'wb'))

    return result, tag
