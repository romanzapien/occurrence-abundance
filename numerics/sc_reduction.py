#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (PDF - source code)
@author: Roman Zapien-Campos - 2021
"""

### Import packages ###
import numpy as np
from  scipy.optimize import least_squares


### Gini index of a vector ###
def gini(p_n):
    
    # number of types
    s = len(p_n)
    
    # normalization
    p = np.array(p_n) / p_n.sum()
    
    gini_index = 0
    
    # compute gini index of p
    for i1 in range(s):
        
        for i2 in range(s):
            
            gini_index += abs(p[i1] - p[i2])
            
    gini_index = gini_index / (2 * (s - 1))
    
    return gini_index



### Recurrence equation (for 2 types) ###

# Transition rate of increase
def R_up(m, N, p, gr, dr, n_i, i):
        
    t_up = 0
            
    j = int(not i)
    
    t_up += (m * p[i] + gr[i] * n_i) * dr[j] * (N - n_i) / N
    
    return t_up


# Transition rate of decrease
def R_down(m, N, p, gr, dr, n_i, i):
        
    t_down = 0
    
    j = int(not i)
        
    t_down += (m * (1. - p[i]) + gr[j] * (N - n_i)) * dr[i] * n_i / N
    
    return t_down


# Compute the probability distribution using the recurrence equation
def recurrence_prob_dist(m, N, p, gr, dr, i):
    
    abundance = np.arange(int(N + 1))
    
    # vector to store the computed probabilities
    prob = np.zeros(int(N + 1))
    # arbitrary selection of a boundary value
    prob[0] = 1
    
    # recurrence towards the opposite boundary
    for n_i in abundance[1:]:
    
        prob[n_i] = prob[n_i - 1] * ( R_up(m, N, p, gr, dr, n_i - 1, i) / R_down(m, N, p, gr, dr, n_i, i))
    
    # normalization
    prob = prob / prob.sum()
    
    return prob



### Model reduction (for S types) ###

# Rate of increase of type i
def R_i_up(m, p, gR, n_i, i):
    
    return m * p[i] + gR[i] * n_i


# Rate of decrease of type i
def R_i_down(dR, N, n_i, i):
    
    return dR[i] * n_i / N


# Rate of increase of type j
def R_j_up(m, p, gR, n_j, j):
    
    return m * p[j] + gR[j] * n_j


# Rate of decrease of type j
def R_j_down(dR, N, n_j, j):
    
    return dR[j] * n_j / N
    

# Step of expected abundace of others computation
def sys_eq_nj(n_j_tofit, n_i, P_ni, n_j, i, m, N, p, gR, dR, S):
    
    j = np.arange(S) != i
    
    if n_i == N:
    
        # Set of equations
        eq_list = n_j_tofit * R_j_down(dR, N, n_j_tofit, j).sum() - R_j_down(dR, N, n_j_tofit, j)
        
    else:
        
        # Probability balance equation (constant)
        C_pb = - P_ni[n_i] * ( R_i_up(m, p, gR, n_i, i) * R_j_down(dR, N, n_j[:, n_i], j).sum() + R_i_down(dR, N, n_i, i) * R_j_up(m, p, gR, n_j[:, n_i], j).sum() ) + P_ni[n_i + 1] * R_i_down(dR, N, n_i + 1, i) * R_j_up(m, p, gR, n_j[:, n_i + 1], j).sum()
        
        # Mean conditioned value equation (constant)
        C_mv = - n_j[:, n_i] * P_ni[n_i] * ( R_i_up(m, p, gR, n_i, i) * R_j_down(dR, N, n_j[:, n_i], j).sum() + R_i_down(dR, N, n_i, i) * R_j_up(m, p, gR, n_j[:, n_i], j).sum() ) + P_ni[n_i + 1] * R_i_down(dR, N, n_i + 1, i) * ( n_j[:, n_i + 1] * R_j_up(m, p, gR, n_j[:, n_i + 1], j).sum() + R_j_up(m, p, gR, n_j[:, n_i + 1], j) ) 
    
        # Set of equations
        eq_list = C_mv - C_pb * ( n_j_tofit - R_j_down(dR, N, n_j_tofit, j) / R_j_down(dR, N, n_j_tofit, j).sum() )
    
    return eq_list


# Step of probability computation
def P_compute(n_i, P_ni, n_j, i, m, N, p, gR, dR, S):
    
    j = np.arange(S) != i
    
    if n_i == N:
        
        # Probability balance equation (constant)
        C_pb = - P_ni[n_i] * R_i_down(dR, N, n_i, i) * R_j_up(m, p, gR, n_j[:, n_i], j).sum()
        
    else:
        
        # Probability balance equation (constant)
        C_pb = - P_ni[n_i] * ( R_i_up(m, p, gR, n_i, i) * R_j_down(dR, N, n_j[:, n_i], j).sum() + R_i_down(dR, N, n_i, i) * R_j_up(m, p, gR, n_j[:, n_i], j).sum() ) + P_ni[n_i + 1] * R_i_down(dR, N, n_i + 1, i) * R_j_up(m, p, gR, n_j[:, n_i + 1], j).sum()

    # Probability balance equation (unknowns)
    return - C_pb / ( R_i_up(m, p, gR, n_i - 1, i) * R_j_down(dR, N, n_j[:, n_i - 1], j).sum() )


# Compute the probability distribution using the model reduction method
def reduction_prob_dist(m, N, p, gR, dR, S, i):
    
    # vector to store the computed probabilities of type i
    P_ni = np.zeros(N+1)
    # vector to store the expected abundances of other types, conditioned on type i
    n_j = np.zeros((S-1, N+1))
    
    # selection of the boundary n_i = N
    starting_point = N
    n_i = starting_point
    
    # recurrence equation towards the boundary n_i = 0
    while n_i > 1:
        
        # try starting from the value n_i = N
        try:
            
            # Set the values at the boundary            
            P_ni[starting_point+1:] = 0
            P_ni[starting_point] = 1
        
            # Iterate towards n_i = 0
            for n_i in np.arange(starting_point, 0, -1):
            
                # Initial guess
                initial_guess = n_j[:, n_i] + 1 / (S - 1)
            
                # Solve system of equations
                n_j[:, n_i - 1] = least_squares(sys_eq_nj, x0 = initial_guess, args = (n_i, P_ni, n_j, i, m, N, p, gR, dR, S), jac = 'cs', loss = 'cauchy').x
                
                # Compute the probability
                P_ni[n_i - 1] = P_compute(n_i, P_ni, n_j, i, m, N, p, gR, dR, S)
                  
            # Normalization                      
            P_ni = P_ni/P_ni.sum()
        
        # otherwise try from a boundary n_i < N
        except:
            
            starting_point -= 1
            # print(starting_point)
            
    return P_ni, n_j



### Compute observables ###
    
# Probability of occurrence of type i
def prob_occurrence(P_ni):

    # Probability of occurrence
    p_occurrence = 1 - P_ni[0]
    
    return p_occurrence


# Expected abundance of type i
def mean_abundance(P_ni, N):
    
    # Mean abundance
    mean_abundance = (P_ni * np.arange(N+1)).sum()
    
    return mean_abundance