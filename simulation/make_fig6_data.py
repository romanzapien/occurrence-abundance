#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (make data for Fig 6)
@author: Roman Zapien-Campos - 2021
"""

# Import packages

from sc_gillespie import gillespie

import numpy as np


### Compute Data ###

# Please comment/uncomment the following lines according to your needs

# Fig6B.1
panel = '6B.1'
i = 11
from par_gillespie_fig5C import *

# # Fig6B.2
# panel = '6B.2'
# i = 37
# from par_gillespie_fig5C import *

# # Fig6C.1
# panel = '6C.1'
# i = 11
# from par_gillespie_fig5C import *

# # Fig6C.2
# panel = '6C.2'
# i = 37
# from par_gillespie_fig5C import *

# # Fig6D.1
# panel = '6D.1'
# i = 11
# from par_gillespie_fig5C import *

# # Fig6D.2
# panel = '6D.2'
# i = 37
# from par_gillespie_fig5C import *


# Values of migration, growth, death rate to test
p_change = np.logspace(-3,0,11)[:-1]
gR_change = np.linspace(0.8,1.2,10)
dR_change = np.linspace(0.8,1.2,10)

# Vectors to store the computed values
P0_nneutral = np.zeros(10)
mean_freq_nneutral = np.zeros(10)

p_ = np.copy(p)

# Compute the observables of each type
for k in range(10):

    # Specific parameter to test
    if panel[:2] == '6B':
        p = p_ / (1. - p_[i]) * (1. - p_change[k])
        p[i] = p_change[k]
        
    elif panel[:2] == '6C':
        gR[i] = gR_change[k]        
    
    elif panel[:2] == '6D':
        dR[i] = dR_change[k]
    
    # Arguments for the gillespie simulation
    args = (N, m, time_sim, N_taxa, N_hosts, n_timepoints, p, gR, dR)
    
    # Gillespie simulation
    timeseries_time, timeseries_data = gillespie(args)
    
    # Mean frequency at equilibrium
    mean_freq_nneutral[k] = timeseries_data[:,-1,i].sum(0)/(N_hosts * N)
    
    # Probability of occurrence at equilibrium
    P0_nneutral[k] = (timeseries_data[:,-1,i]>0).sum(0)/N_hosts

# Save data according to its kind
if panel[:2] == '6B': np.savez_compressed('../data/%s.npz'%panel, N = N, m = m, time_sim = time_sim, S = N_taxa, N_hosts = N_hosts, n_timepoints = n_timepoints, gR = gR, dR = dR, P0_nneutral = P0_nneutral, mean_freq_nneutral = mean_freq_nneutral, p_change = p_change)
elif panel[:2] == '6C': np.savez_compressed('../data/%s.npz'%panel, N = N, m = m, time_sim = time_sim, S = N_taxa, N_hosts = N_hosts, n_timepoints = n_timepoints, p = p, dR = dR, P0_nneutral = P0_nneutral, mean_freq_nneutral = mean_freq_nneutral, gR_change = gR_change)
elif panel[:2] == '6D': np.savez_compressed('../data/%s.npz'%panel, N = N, m = m, time_sim = time_sim, S = N_taxa, N_hosts = N_hosts, n_timepoints = n_timepoints, p = p, gR = gR, P0_nneutral = P0_nneutral, mean_freq_nneutral = mean_freq_nneutral, dR_change = dR_change)