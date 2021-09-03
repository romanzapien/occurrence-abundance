#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (make data of Fig 6B-D - source code)
@author: Roman Zapien-Campos - 2021
"""

# Import functions
from sc_reduction import *


# Choose the panel to compute
panel = '6D.1' # Options are Type 1: '6B.1', '6C.1', '6D.1'; Type 2: '6B.2', '6C.2', '6D.2'


### Parameters ###

## Fraction of immigrants
# Gini 0.6
p_ = np.array([0.00038154, 0.00067294, 0.00077313, 0.00084858, 0.00126614, 0.00137603, 0.00170499, 0.00173883, 0.00193053, 0.00216253, 0.00437584, 0.00527042, 0.00728132, 0.00754192, 0.00775608, 0.00803994, 0.01048105, 0.01178707, 0.01260198, 0.01373537, 0.01404538, 0.01485447, 0.01632661, 0.01833767, 0.02257611, 0.02528464, 0.02540005, 0.02571126, 0.02691791, 0.03328889, 0.03412534, 0.03504708, 0.03723251, 0.0373183, 0.0412504, 0.06049704, 0.07532528, 0.07580609, 0.09674053, 0.18218821])

## Growth rates
# SD 0.1
gR_ = np.array([1.15442723, 0.85911785, 1.11160448, 1.02066184, 0.92910304, 0.89754369, 0.81035959, 1.0647798, 0.98090045, 1.04132257, 0.85088964, 1.04953219, 1.05582736, 0.82319735, 1.09151151, 1.15998802, 0.82941692, 1.0370813, 0.95696431, 0.99065505, 0.91642347, 1.02392111, 1.02467767, 0.96302918, 1.01826739, 1.04124979, 1.07662959, 1.20783799, 0.97987992, 1.14226125, 1.00054386, 0.94077612, 1.00486504, 0.96320097, 1.14489223, 1.02131939, 0.91426701, 0.88496494, 0.89403731, 1.06832477])

## Death rates
# SD 0.1
dR_ = np.array([0.97997334, 0.97178519, 1.01221396, 1.07464284, 0.89822532, 1.05447841, 0.90908804, 0.93517044, 1.11642025, 1.0223452, 0.96517128, 1.00578925, 0.9853043, 0.82862568, 0.98050893, 1.06342287, 1.11187164, 1.11135766, 1.16597829, 1.12204569, 0.96069496, 1.01766923, 0.97712838, 1.07470396, 0.85933591, 1.09789077, 0.94291763, 0.87396482, 0.9675674, 1.19007246, 1.09704941, 0.91754573, 0.88923551, 0.97910369, 1.00315385, 1.01284054, 1.00109989, 0.76639059, 1.0966333, 0.97279744])

# Overall prameters
N = int(1E3)
m = 1E-1 * N
S = 40
p = p_
gR = gR_
dR = dR_   

# Vectors to store the computed values
P0_nneutral = np.zeros(10)
mean_freq_nneutral = np.zeros(10)

# Values of migration, growth, death to test
p_change = np.logspace(-3,0,11)[:-1]
gR_change = np.linspace(0.8,1.2,10)
dR_change = np.linspace(0.8,1.2,10)

# print(panel)

# Compute the observables of each type
for i in range(10):

    # Focal type
    if panel[-1] == '1':
        j = 2
    elif panel[-1] == '2':
        j = 14
    
    # Panel specific value to test
    if panel[:2] == '6B':
        p = p_ / (1. - p_[j]) * (1. - p_change[i])
        p[j] = p_change[i]
    elif panel[:2] == '6C':
        gR[j] = gR_change[i]        
    elif panel[:2] == '6D':
        dR[j] = dR_change[i]
            
    # print(i)
        
    # Compute the probability distribution of type i
    P_ni_nneutral, n_j_nneutral = reduction_prob_dist(m, N, p, gR, dR, S, j)
    
    # Compute the occurrence and mean frequencies 
    P0_nneutral[i] = prob_occurrence(P_ni_nneutral)
    mean_freq_nneutral[i] = mean_abundance(P_ni_nneutral, N) / N
    
    # print(n_j_nneutral[:,0].sum(), '\n')

# Saved the compute values of the chosen panel
if panel[:2] == '6B': np.savez_compressed('data/%s.npz'%panel, N = N, m = m, S = S, p = p_, gR = gR_, dR = dR_, P0_nneutral = P0_nneutral, mean_freq_nneutral = mean_freq_nneutral, p_change = p_change)
elif panel[:2] == '6C': np.savez_compressed('data/%s.npz'%panel, N = N, m = m, S = S, p = p_, gR = gR_, dR = dR_, P0_nneutral = P0_nneutral, mean_freq_nneutral = mean_freq_nneutral, gR_change = gR_change)
elif panel[:2] == '6D': np.savez_compressed('data/%s.npz'%panel, N = N, m = m, S = S, p = p_, gR = gR_, dR = dR_, P0_nneutral = P0_nneutral, mean_freq_nneutral = mean_freq_nneutral, dR_change = dR_change)