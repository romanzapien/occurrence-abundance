#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (recurrence eq. - source code)
@author: Roman Zapien-Campos - 2021
"""

### Import packages ###
import numpy as np


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