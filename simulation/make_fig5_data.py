#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (make data for Fig 5)
@author: Roman Zapien-Campos - 2021
"""

# Import packages

from sc_gillespie import gillespie

import numpy as np


### Compute Data ###

# Please comment/uncomment the following lines according to your needs

# Fig5A
panel = '5A'
from par_gillespie_fig5A import *

# # Fig5B
# panel = '5B'
# from par_gillespie_fig5B import *

# # Fig5C
# panel = '5C'
# from par_gillespie_fig5C import *

# # Fig5D
# panel = '5D'
# from par_gillespie_fig5D import *

# Arguments for the gillespie simulation
args = (N, m, time_sim, N_taxa, N_hosts, n_timepoints, p, gR, dR)

# Gillespie simulation
timeseries_time, timeseries_data = gillespie(args)

# Mean frequency at equilibrium
mean_freq_nneutral = timeseries_data[:,-1,:].sum(0)/(N_hosts * N)

# Probability of occurrence at equilibrium
P0_nneutral = (timeseries_data[:,-1,:]>0).sum(0)/N_hosts

# Save data
np.savez_compressed('../data/%s.npz'%panel, N = N, m = m, time_sim = time_sim, S = N_taxa, N_hosts = N_hosts, n_timepoints = n_timepoints, p = p, gR = gR, dR = dR, timeseries_time = timeseries_time, timeseries_data = timeseries_data, P0_nneutral = P0_nneutral, mean_freq_nneutral = mean_freq_nneutral)