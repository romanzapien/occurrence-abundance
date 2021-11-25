#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (Gillespie - source code)
@author: Roman Zapien-Campos - 2021
"""

### Import packages ###

import numpy as np

from numba import jit

from random import random


# Compute the Gini index of a vector
def gini(p_n):
    
    s = len(p_n)
    
    p = np.array(p_n) / p_n.sum()
    
    gini_index = 0
    
    for i1 in range(s):
        
        for i2 in range(s):
            
            if i1 != i2: gini_index += abs(p[i1] - p[i2])
            
    gini_index = gini_index / (2 * (s - 1))
    
    return gini_index


# Get a vector of frequencies in the source pool numerators (pn)
def get_pn_sample(gini_target, S):

    resolution = int(1E4)
    
    # Initial condition    
    
    probs = np.random.random(S)
    
    probs = probs / probs.sum()
    
    p_n = np.random.multinomial(S * resolution, probs)
        
    # Initial Gini index
        
    gini_index = gini(p_n)
    
    # Iterate until the Gini target is achieved
    
    while gini_index < 0.99 * gini_target or gini_index > 1.01 * gini_target:
        
        index_up = np.random.choice(np.arange(S)[p_n < resolution * S])
        
        index_down = np.random.choice(np.arange(S)[p_n > 0])
        
        p_n[index_down] -= 1
        
        p_n[index_up] += 1
        
        if abs(gini_target - gini(p_n)) < abs(gini_target - gini_index): 
            
            gini_index = gini(p_n)
    
        else: 
            
            p_n[index_down] += 1
        
            p_n[index_up] -= 1
                    
    return p_n/(S*resolution)


# Compute the rates of the current state: all hosts (ah)
#@jit(nopython = True)
def compute_rates_ah(n, growth_rate, death_rate, N, m, p, N_taxa, N_hosts):
    
    death = np.transpose(death_rate * (n / N), (0, 2, 1))
    
    immigration = m * p
    
    birth = growth_rate * n
    
    tnRs = np.einsum('ijk, ijk -> ijk', death, immigration + birth)
        
    return tnRs


# Compute the rates of the current state: single host (sh)
@jit(nopython = True)
def compute_rates_sh(n, growth_rate, death_rate, N, m, p, N_taxa, host):
    
    death = (death_rate * (n[host] / N)).T
    
    immigration = m * p
    
    birth = growth_rate * n[host]
    
    tnRs = np.outer(death, immigration + birth)
        
    return tnRs


# Compute time and choice parameters
@jit(nopython = True)
def compute_time_n_choice_par(tnRs):

    time_par = 1. / tnRs.sum()
    
    choice_par = tnRs * time_par

    return time_par, choice_par


# Sample the time
@jit(nopython = True)
def time_sample(time_par):

    time_sampled = np.random.exponential(time_par)
    
    return time_sampled


# Sample the rates
@jit(nopython = True)
def sample_discrete(time_par, choice_par, N_taxa):

    q = random()
    
    p_sum = 0
            
    taxon_decrease = -1
    
    taxon_increase = -1
        
    i = 0

    while p_sum < q:
        
        taxon_decrease = i // (N_taxa)
        
        taxon_increase = i % (N_taxa)
            
        p_sum += choice_par[taxon_decrease, taxon_increase]
        
        i += 1
        
    return taxon_decrease, taxon_increase


# Gillespie algorithm: all hosts simultaneously
def gillespie(args):

    # Unpack parameters

    [N, m, time_sim] = [np.float(arg) for arg in args[:3]]
    
    [N_taxa, N_hosts, n_timepoints] = [np.int(arg) for arg in args[3:6]]
    
    p, growth_rate, death_rate = args[6:]
    
    p = p / p.sum()
        
    # Initial condition
    
    init_cond = np.zeros(N_taxa)
        
    for n in range(int(N)):
        
        init_cond[np.random.choice(N_taxa,1,p=p)] += 1
    
    # Create and store data if requeried
    
    if n_timepoints != 0:
        
        # Create arrays to store results
    
        timeseries_time = np.logspace(0, np.log10(time_sim), n_timepoints - 1, dtype = np.float32)
        
        timeseries_time = np.hstack((0, timeseries_time))
        
        timeseries_data = np.zeros((N_hosts, n_timepoints, N_taxa), dtype = np.float32)
                        
        step = np.ones(N_hosts, dtype = np.int)
        
        # Write initial conditions in the storing arrays
        
        timeseries_data[:, 0, :] = init_cond
    
    # Create operational arrays

    time = np.zeros(N_hosts)
    
    n = init_cond * np.ones((N_hosts, 1, N_taxa), dtype = np.float32)
    
    sim_hosts = np.arange(N_hosts, dtype = np.int)
        
    # Compute initial transition rates for all hosts
    
    tnRs = compute_rates_ah(n, growth_rate, death_rate, N, m, p, N_taxa, N_hosts)
    
    tnRs[:, range(N_taxa), range(N_taxa)] = 0 # Do not modify, it causes problems !!!
    
    # While loop to compute the system dynamics until time_sim is reached

    while sim_hosts.shape[0] > 0:
        
        for host in sim_hosts:
        
            time_par, choice_par = compute_time_n_choice_par(tnRs[host])
            
            taxon_decrease, taxon_increase  = sample_discrete(time_par, choice_par, N_taxa)
            
            time[host] += time_sample(time_par)
                                    
            if n_timepoints != 0:
                
                while time[host] > timeseries_time[step[host]]:
                                        
                    timeseries_data[host, step[host], :] = n[host]
                    
                    step[host] += 1
                    
                    if time[host] >= time_sim: 
                        
                        sim_hosts = sim_hosts[sim_hosts != host]
                        
                        timeseries_data[host, step[host]:, :] = n[host]
                        
                        break
                    
            n[host, 0, taxon_decrease] -= 1
                
            n[host, 0, taxon_increase] += 1
    
            tnRs[host] = compute_rates_sh(n, growth_rate, death_rate, N, m, p, N_taxa, host)
            
            tnRs[host][np.diag_indices(N_taxa)] = 0
            
    return timeseries_time, timeseries_data