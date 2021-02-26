import h5py
import numpy as np
from parse import parse

def rl(rf):
    return " ".join(rf.readline().split())


with h5py.File('datasets.hdf5', 'w') as f:
    with open('datasets.dat', 'r') as rf:
        while(1):
            line = rl(rf)
            if not line or not line.startswith('#START_PACKET'):
                break
            cycle_id = parse('#CYCLE_ID: {:d}', rl(rf))[0]
            slater_matrix_dim = parse('#SLATER_MATRIX_DIM: {:d}', rl(rf))[0]
            nupdates = parse('#NUPDATES: {:d}', rl(rf))[0]
            assert(rf.readline().startswith('#SLATER_MATRIX'))

            # Read matrices
            slater_matrix = np.zeros((slater_matrix_dim,slater_matrix_dim))
            slater_inverse = np.zeros((slater_matrix_dim,slater_matrix_dim))
            for i in range(slater_matrix_dim*slater_matrix_dim):
                res = parse('({i:d},{j:d}) {sla:e} {inv:e}', rl(rf))
                slater_matrix[res['i']-1, res['j']-1] = res['sla']
                slater_inverse[res['i']-1, res['j']-1] = res['inv']

            # Read updates
            col_update_index = np.zeros(nupdates, dtype='i')
            updates = np.zeros((nupdates, slater_matrix_dim))
            for n in range(nupdates):
                col_update_index[n] = parse('#COL_UPDATE_INDEX: {:d}', rl(rf))[0]
                for i in range(slater_matrix_dim):
                    res = parse('#COL_UPDATE_COMP_({i:d}): {x:e}', rl(rf))
                    updates[n][res['i']-1] = res['x']

            assert(rf.readline().startswith('#END_PACKET'))
            rf.readline()

            cycle = f.create_group('cycle_{}'.format(cycle_id))
            cycle.create_dataset("slater_matrix_dim", data=slater_matrix_dim)
            cycle.create_dataset("nupdates", data=nupdates)
            cycle.create_dataset("slater_matrix", data=slater_matrix, compression='gzip')
            cycle.create_dataset("slater_inverse", data=slater_inverse, compression='gzip')
            cycle.create_dataset("col_update_index", data=col_update_index)
            cycle.create_dataset("updates", data=updates, compression='gzip')
