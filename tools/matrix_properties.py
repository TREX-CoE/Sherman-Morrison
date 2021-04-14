#!/usr/bin/env python
import h5py
import sys
import numpy as np
import math

fname = sys.argv[1]
dataset =h5py.File(fname, "r")
cycle_num = sys.argv[2]
cycle  = "cycle_" + cycle_num
cgroup = dataset[cycle]

slater_matrix = cgroup["slater_matrix"][()]
slater_inverse = cgroup["slater_inverse"][()]
nupdates = cgroup["nupdates"][()]
update_idxs = cgroup["col_update_index"][()]
updates = cgroup["updates"][()]
dim = cgroup["slater_matrix_dim"][()]
Id = np.identity(dim)

print(f'========================\n== UPDATE CYCLE: {cycle_num} ==\n========================')
print()
print(f'Number of updates: {nupdates}')
print(f'Columns that need to be updated (replaced): {update_idxs}\n')

# det_sm = np.linalg.det(slater_matrix)
# det_si = np.linalg.det(slater_inverse)
# eigen_values_sm, eigen_vectors_sm = np.linalg.eig(slater_matrix)
# eigen_values_si, eigen_vectors_si = np.linalg.eig(slater_inverse)

# SSi = np.matmul(slater_matrix, slater_inverse)
# detSSi = np.linalg.det(SSi)
# delta = Id - SSi
# res_max = np.amax(abs(delta))
# res_max2 = res_max*res_max
# res_fro = np.linalg.norm(delta, 'fro')
# res_fro2 = res_fro*res_fro

# print(f'Determinant of Slater-matrix: {det_sm}')
# print(f'Determinant of inverse Slater-matrix: {det_si}\n')

# print(f'Smalles eigenvalue (modulus) of Slater-matrix: {np.amin(abs(eigen_values_sm))}')
# print(f'Smalles eigenvalue (modulus) of inverse Slater-matrix: {np.amin(abs(eigen_values_si))}\n')

# print(f'Det(S*Sinv): {detSSi}')
# print(f'Residual (Id - S*Sinv) (Max norm): {res_max}')
# print(f'Residual (Max norm sq.): {res_max2}')
# print(f'Residual (Frob. norm): {res_fro}')
# print(f'Residual (Frob. norm sq.): {res_fro2}\n')

print(f'+-----------------------------------------------+\n| Updating Slater-matrix...                     |')
slater_matrix_new = slater_matrix
index = 0
for update in updates:
    column  = update_idxs[index] - 1
    slater_matrix_new[column] = update
    index = index + 1

print(f'| Computing inverse of updated Slater-matrix... |\n+-----------------------------------------------+\n')
slater_inverse_new = np.linalg.inv(slater_matrix_new)

det_sm = np.linalg.det(slater_matrix_new)
det_si = np.linalg.det(slater_inverse_new)
eigen_values_sm, eigen_vectors_sm = np.linalg.eig(slater_matrix_new)
eigen_values_si, eigen_vectors_si = np.linalg.eig(slater_inverse_new)

SSi = np.matmul(slater_matrix_new, slater_inverse_new)
detSSi = np.linalg.det(SSi)
delta = Id - SSi
res_max = np.amax(abs(delta))
res_max2 = res_max*res_max
res_fro = np.linalg.norm(delta, 'fro')
res_fro2 = res_fro*res_fro

print(f'(DUSM) Determinant of UPDATED Slater-matrix: {cycle_num}, {det_sm}')
print(f'(DUSI) Determinant of UPDATED inverse Slater-matrix: {det_si}\n')

print(f'(SEVUSM) Smalles eigenvalue (modulus) of UPDATED Slater-matrix: {np.amin(abs(eigen_values_sm))}')
print(f'(SEVUSI) Smalles eigenvalue (modulus) of UPDATED inverse Slater-matrix: {np.amin(abs(eigen_values_si))}\n')

print(f'(DSSI) det(S*Sinv): {detSSi}')
print(f'(RMN) Residual (Id - S*Sinv) (Max norm): {res_max}')
print(f'(RMN2) Residual (Max norm sq.): {res_max2}')
print(f'(RFN) Residual (Frob. norm): {res_fro}')
print(f'(RFN2) Residual (Frob. norm sq.): {res_fro2}\n')

dataset.close()