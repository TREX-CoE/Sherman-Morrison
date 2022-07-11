import numpy as np
import h5py

h5file = h5py.File('dataset_zeropadded_cm.hdf5', 'r')
print(f"#cycle, det, cond_2, cond_f, norm_2, norm_f")
print(f"# 1, 2, 3, 4, 5, 6")
for key in h5file.keys():
    cycle = h5file.get(key)
    slater_matrix = cycle.get('slater_matrix')
    slater_matrix = np.array(slater_matrix )
    slater_inverse_t = cycle.get('slater_inverse_t')

    slater_inverse_t = np.array(slater_inverse_t)
    slater_inverse = slater_inverse_t.transpose()
    slater_matrix_sq = slater_matrix[:, 0:21]
    slater_inverse_sq = slater_inverse[0:21, :]

    det = np.linalg.det(slater_matrix_sq)
    cond_2 = np.linalg.cond(slater_inverse_sq, p=2)
    cond_f = np.linalg.cond(slater_inverse_sq, p='fro')

    Id_appr = np.matmul(slater_matrix_sq, slater_inverse_sq)
    Err  = Id_appr - np.identity(Id_appr.shape[0])
    normf = np.linalg.norm(Err, ord='fro')
    norm2 = np.linalg.norm(Err, ord=2)
    print(f"{key[6:]}, {det}, {cond_2}, {cond_f}, {norm2}, {normf}")
h5file.close()
