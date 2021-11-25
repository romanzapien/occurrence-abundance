#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (parameters for Fig 4D - source code)
@author: Roman Zapien-Campos - 2021
"""

# Import packages

import numpy as np


### Define parameters ###

# General parameters

# Number of microbes within a host
N = 1E3

# Number of microbial taxa
N_taxa = 6E0

# Migration rate
m = 1E-1 * N

# Number of hosts
N_hosts = 1E3

# Time to be simulated
time_sim = 1E3

# Number of timepoints to store
n_timepoints = 1E2

# Fraction of immigrants
p = np.array([0.33443333, 0.00168333, 0.12375, 0.08033333, 0.4176, 0.0422])

# Growth rate
gR = np.array([0.93815818, 1.03693635, 0.82649531, 1.0933327, 1.04066052, 0.80158416])
  
# Death rate
dR = np.array([1.00512998, 0.82321306, 1.00487796, 1.00530296, 0.94770017, 0.84279409])