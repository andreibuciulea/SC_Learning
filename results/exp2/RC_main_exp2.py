# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 17:27:42 2023

@author: Andrei
"""

import gudhi
import numpy as np
from scipy.io import loadmat
from scipy.io import savemat

nme = 'corr_data_exp2.mat'
annots = loadmat(nme)
Cors = annots['all_C_X0_X1']


#%

nG = 100
M = 9
N = 20
all_H1 = np.zeros((N,N,nG,M))
all_H2 = np.zeros((N,N**2,nG,M))
all_B2 = np.zeros((3,1350,nG,M))
all_T2_vals = np.zeros((1350,nG,M))
T2_vals = np.zeros((1350))
B2 = np.zeros((3,1350))
for g in range(nG):
    for m in range(M): 
        correlation_matrix = Cors[:,:,g,m]
        distance_matrix = 1 - correlation_matrix
        rips_complex = gudhi.RipsComplex(distance_matrix=distance_matrix, max_edge_length=1.0,sparse=2)
        
        simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
        result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
            repr(simplex_tree.num_simplices()) + ' simplices - ' + \
            repr(simplex_tree.num_vertices()) + ' vertices.'
        print(result_str)
        fmt = '%s -> %.2f'
        H1 = np.zeros((N,N))
        H2 = np.zeros((N,N**2))
        i = 0
        for filtered_value in simplex_tree.get_filtration():
            link = filtered_value[0]
            if len(link) == 2:
                H1[link[0],link[1]] = filtered_value[1]        
                H1[link[1],link[0]] = filtered_value[1]   
            elif len(filtered_value[0])==3:
                link = np.array(link)-1
                H2[link[0],N*link[1]+link[2]] = filtered_value[1]
                H2[link[0],N*link[2]+link[1]] = filtered_value[1]
                H2[link[1],N*link[0]+link[2]] = filtered_value[1]
                H2[link[1],N*link[2]+link[0]] = filtered_value[1]  
                H2[link[2],N*link[0]+link[1]] = filtered_value[1]
                H2[link[2],N*link[1]+link[0]] = filtered_value[1]  
                
                B2[:,i] = link+2
                T2_vals[i] = filtered_value[1]
                i+=1
                
        all_H1[:,:,g,m] = H1
        all_H2[:,:,g,m] = H2
        all_B2[:,:,g,m] = B2
        all_T2_vals[:,g,m] = T2_vals
        

res_nme = 'RC_results_exp2_sparse.mat'
savemat(res_nme, {'all_H1': all_H1, 'all_H2':all_H2,'all_T2_hat':all_B2,'all_T2_vals':all_T2_vals})

#%%

# User defined correlation matrix is:
# |1     0.06    0.23    0.01    0.89|
# |0.06  1       0.74    0.01    0.61|
# |0.23  0.74    1       0.72    0.03|
# |0.01  0.01    0.72    1       0.7 |
# |0.89  0.61    0.03    0.7     1   |
#correlation_matrix=np.array([[1., 0.06, 0.23, 0.01, 0.89],
#                            [0.06, 1., 0.74, 0.01, 0.61],
#                            [0.23, 0.74, 1., 0.72, 0.03],
#                            [0.01, 0.01, 0.72, 1., 0.7],
#                            [0.89, 0.61, 0.03, 0.7, 1.]], float)

nG = 100
M = 9
N = 20
all_H1 = np.zeros((N,N,nG,M))
all_H2 = np.zeros((N,N**2,nG,M))
all_B2 = np.zeros((3,1350,nG,M))
all_T2_vals = np.zeros((1350,nG,M))
T2_vals = np.zeros((1350))
B2 = np.zeros((3,1350))
for g in range(nG):
    for m in range(M): 
        correlation_matrix = Cors[:,:,g,m]
        distance_matrix = 1 - correlation_matrix
        rips_complex = gudhi.RipsComplex(distance_matrix=distance_matrix, max_edge_length=1.0)
        
        simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
        result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
            repr(simplex_tree.num_simplices()) + ' simplices - ' + \
            repr(simplex_tree.num_vertices()) + ' vertices.'
        print(result_str)
        fmt = '%s -> %.2f'
        H1 = np.zeros((N,N))
        H2 = np.zeros((N,N**2))
        i = 0
        for filtered_value in simplex_tree.get_filtration():
            link = filtered_value[0]
            if len(link) == 2:
                H1[link[0],link[1]] = filtered_value[1]        
                H1[link[1],link[0]] = filtered_value[1]   
            elif len(filtered_value[0])==3:
                link = np.array(link)-1
                H2[link[0],N*link[1]+link[2]] = filtered_value[1]
                H2[link[0],N*link[2]+link[1]] = filtered_value[1]
                H2[link[1],N*link[0]+link[2]] = filtered_value[1]
                H2[link[1],N*link[2]+link[0]] = filtered_value[1]  
                H2[link[2],N*link[0]+link[1]] = filtered_value[1]
                H2[link[2],N*link[1]+link[0]] = filtered_value[1]  
                
                B2[:,i] = link+2
                T2_vals[i] = filtered_value[1]
                i+=1
                
        all_H1[:,:,g,m] = H1
        all_H2[:,:,g,m] = H2
        all_B2[:,:,g,m] = B2
        all_T2_vals[:,g,m] = T2_vals
        
res_nme = 'RC_results_exp2.mat'
savemat(res_nme, {'all_H1': all_H1, 'all_H2':all_H2,'all_T2_hat':all_B2,'all_T2_vals':all_T2_vals})


