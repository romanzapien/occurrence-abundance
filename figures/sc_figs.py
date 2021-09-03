#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: occurrence-abundance pattern (Figures - source code)
@author: Roman Zapien-Campos - 2021
"""

# Import packages
import numpy as np
import matplotlib.pyplot as mp
from matplotlib import colors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import LogLocator, NullFormatter

import sys
sys.path.insert(1, '../numerics/')
from sc_reduction import *


### Figure 1 ###
# Expected equilibrium of a type if rates in the community are neutral.
def fig1():
    
    fig, axes = mp.subplots(nrows = 1, ncols = 3, figsize = (11.1, 3.2))
    
    ## Fig 1A
    
    focal_type = 0
    
    N = 1E3
    p = [1/2, 1/2]
    gr = [1., 1.]
    dr = [1., 1.]
    
    abundance = np.arange(int(N + 1)) / N
    
    # Low migration
    
    m = 1E-3 * N
    prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    axes[0].plot(abundance, prob, '-', color = 'lightseagreen', label = r'$10^{-3}$')
    
    # Large migration
    
    m = 1E-1 * N
    prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    axes[0].plot(abundance, prob, '-', color = 'darkorchid', label = r'$10^{-1}$')
    
    # Annotate
    
    axes[0].legend(title=r'immigration ($m/N$)')

    axes[0].set_xlabel('frequency', fontsize = 16)
    axes[0].set_ylabel('probability density', fontsize = 16)

    ## Fig 1B-C
    
    abundance = np.arange(int(N + 1))
    P = np.logspace(-4,0,50)
    
    # Low migration
    
    m = 1E-3 * N
    
    P0 = np.zeros(len(P))
    mean_freq = np.zeros(len(P))
    
    for p_ in P:
    
        p = [p_, 1 - p_]
    
        prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
        
        P0[p_ == P] = 1. - prob[0]
        mean_freq[p_ == P] = sum(prob * abundance/N)

    mean_freq[-1] = 1

    axes[1].loglog(P, mean_freq, '-', color = 'lightseagreen', label = r'$10^{-3}$')
    axes[2].semilogx(P, P0, '-', color = 'lightseagreen', label = r'$10^{-3}$')

    # Large migration

    m = 1E-1 * N
    
    P0 = np.zeros(len(P))
    mean_freq = np.zeros(len(P))
    
    for p_ in P:
    
        p = [p_, 1 - p_]
    
        prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
            
        P0[p_ == P] = 1. - prob[0]
        mean_freq[p_ == P] = sum(prob * abundance/N)

    mean_freq[-1] = 1

    axes[1].loglog(P, mean_freq, '-', color = 'darkorchid', label = r'$10^{-1}$')
    axes[2].semilogx(P, P0, '-', color = 'darkorchid', label = r'$10^{-1}$')    

    # Annotate

    fig.subplots_adjust(wspace=0.4)

    axes[0].ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)
    axes[0].tick_params(axis='both', direction='in', which='both', right=True, top=True)

    for ax in [axes[1], axes[2]]:

        ax.set_xlabel('fraction of migrants', fontsize = 16)
        ax.legend(title=r'immigration ($m/N$)')
        ax.set_xlim(7E-4,1.5)
        ax.tick_params(axis='both', direction='in', which='both', right=True, top=True)
    
    axes[1].set_ylabel('mean frequency', fontsize = 16) 
    axes[1].set_ylim(7E-4,1.5)
    
    axes[2].set_ylabel('occurrence freq.', fontsize = 16)
    axes[2].set_ylim(-0.05,1.05) 

    axes[0].text(-0.25, 1.03, 'A', transform = axes[0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1].text(-0.25, 1.03, 'B', transform = axes[1].transAxes, fontsize = 25, fontweight = 'bold')
    axes[2].text(-0.25, 1.03, 'C', transform = axes[2].transAxes, fontsize = 25, fontweight = 'bold')
    
    # Save
    
    mp.savefig('fig1.pdf', dpi = 300, bbox_inches='tight')



### Figure 2 ###
# Effect of non-neutral growth rates on the equilibrium of a community with two types.
def fig2():
    
    fig, axes = mp.subplots(nrows = 2, ncols = 2, figsize = (6.6, 6.6))
    
    focal_type = 0
    
    N = 1E3
    dr = [1., 1.]
    
    abundance = np.arange(int(N + 1))
    P = np.logspace(-4,0,50)
    
    # Growth rates and colors and labels associated to them
    Gr_i = np.array([0.8, 0.9, 1., 1.1, 1.2])
    color_gr = np.array(['lightseagreen', 'darkorchid', 'black', 'crimson', 'darkorange'])
    label_gr = np.array(['0.8', '0.9', 'neutral', '1.1', '1.2'])
    
    ## Fig 2A and 2C (low migration)
    
    m = 1E-3 * N

    for gr_i in Gr_i:
        
        P0 = np.ones(len(P))*np.nan
        mean_freq = np.ones(len(P))*np.nan
        
        for p_ in P:
        
            p = [p_, 1. - p_]
            gr = [gr_i, 1.]
            
            prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    
            P0[p_ == P] = 1. - prob[0]
            mean_freq[p_ == P] = sum(prob * abundance/N)
            
        mean_freq[mean_freq == 0] = np.nan

        axes[0,0].loglog(P, mean_freq, '-', color = color_gr[gr_i == Gr_i][0], label = label_gr[gr_i == Gr_i][0])
        axes[1,0].semilogx(P, P0, '-', color = color_gr[gr_i == Gr_i][0], label = label_gr[gr_i == Gr_i][0])    
        
    ## Fig 2B and 2D (large migration)

    m = 1E-1 * N
    
    for gr_i in Gr_i:
    
        P0 = np.ones(len(P))*np.nan
        mean_freq = np.ones(len(P))*np.nan
        
        for p_ in P:
        
            p = [p_, 1. - p_]
            gr = [gr_i, 1.]
            
            prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    
            P0[p_ == P] = 1. - prob[0]
            mean_freq[p_ == P] = sum(prob * abundance/N)
            
        mean_freq[mean_freq == 0] = np.nan
    
        axes[0,1].loglog(P, mean_freq, '-', color = color_gr[gr_i == Gr_i][0], label = label_gr[gr_i == Gr_i][0])
        axes[1,1].semilogx(P, P0, '-', color = color_gr[gr_i == Gr_i][0], label = label_gr[gr_i == Gr_i][0])

    # Annotate

    axes[1,1].legend(title=r'$f_1$', fontsize = 10, bbox_to_anchor=(1.05, 1.54), loc='upper left')
            
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    
    for ax in [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]:
    
        ax.set_xlabel('fraction of migrants', fontsize = 16)    
        ax.set_xlim(7E-4,1.5)
        ax.tick_params(axis='both', direction='in', which='both', right=True, top=True)
    
        if ax == axes[0,0] or ax == axes[0,1]: 
            ax.set_ylabel('mean frequency', fontsize = 16)
            ax.set_ylim(1E-6,1.5)
            
        else:
            ax.set_ylabel('occurrence freq.', fontsize = 16)
            ax.set_ylim(-0.05,1.05)

    axes[0,0].text(-0.25, 1.03, 'A', transform = axes[0,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[0,1].text(-0.25, 1.03, 'B', transform = axes[0,1].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,0].text(-0.25, 1.03, 'C', transform = axes[1,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,1].text(-0.25, 1.03, 'D', transform = axes[1,1].transAxes, fontsize = 25, fontweight = 'bold')

    axes[0,0].text(0.5, 1.15, 'small immigration', transform = axes[0,0].transAxes, ha='center', fontsize = 16, fontweight = 'bold')
    axes[0,1].text(0.5, 1.15, 'large immigration', transform = axes[0,1].transAxes, ha='center', fontsize = 16, fontweight = 'bold')

    # Save

    fig.savefig('fig2.pdf', dpi = 300, bbox_inches='tight')

    

### Figure 3 ###
# Effect of non-neutral death rates on the equilibrium of a community with two types.
def fig3():
    
    fig, axes = mp.subplots(nrows = 2, ncols = 2, figsize = (6.6, 6.6))
    
    focal_type = 0
    
    N = 1E3
    gr = [1., 1.]
    
    abundance = np.arange(int(N + 1))
    P = np.logspace(-4,0,50)

    # Death rates and colors and labels associated to them
    Dr_i = np.array([0.8, 0.9, 1., 1.1, 1.2])
    color_dr = np.array(['darkorange', 'crimson', 'black', 'darkorchid', 'lightseagreen'])
    label_dr = np.array(['0.8', '0.9', 'neutral', '1.1', '1.2'])
    
    ## Fig 3A and 3C (low migration)
    
    m = 1E-3 * N

    for dr_i in Dr_i:
        
        P0 = np.ones(len(P))*np.nan
        mean_freq = np.ones(len(P))*np.nan
        
        for p_ in P:
        
            p = [p_, 1. - p_]
            dr = [dr_i, 1.]
            
            prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    
            P0[p_ == P] = 1. - prob[0]
            mean_freq[p_ == P] = sum(prob * abundance/N)

        mean_freq[mean_freq == 0] = np.nan

        axes[0,0].loglog(P, mean_freq, '-', color = color_dr[dr_i == Dr_i][0], label = label_dr[dr_i == Dr_i][0])
        axes[1,0].semilogx(P, P0, '-', color = color_dr[dr_i == Dr_i][0], label = label_dr[dr_i == Dr_i][0])
        
    ## Fig 3B and 3D (large migration)

    m = 1E-1 * N
    
    for dr_i in Dr_i:
    
        P0 = np.ones(len(P))*np.nan
        mean_freq = np.ones(len(P))*np.nan
        
        for p_ in P:
        
            p = [p_, 1. - p_]
            dr = [dr_i, 1.]
            
            prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    
            P0[p_ == P] = 1. - prob[0]
            mean_freq[p_ == P] = sum(prob * abundance/N)
            
        mean_freq[mean_freq == 0] = np.nan

        axes[0,1].loglog(P, mean_freq, '-', color = color_dr[dr_i == Dr_i][0], label = label_dr[dr_i == Dr_i][0])
        axes[1,1].semilogx(P, P0, '-', color = color_dr[dr_i == Dr_i][0], label = label_dr[dr_i == Dr_i][0])

    # Annotate

    axes[1,1].legend(title=r'$\phi_1$', fontsize = 10, bbox_to_anchor=(1.05, 1.54), loc='upper left')
    
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    
    for ax in [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]:
    
        ax.set_xlabel('fraction of migrants', fontsize = 16)    
        ax.set_xlim(7E-4,1.5)
        ax.tick_params(axis='both', direction='in', which='both', right=True, top=True)
    
        if ax == axes[0,0] or ax == axes[0,1]: 
            ax.set_ylabel('mean frequency', fontsize = 16)
            ax.set_ylim(1E-6,1.5)
            
        else:
            ax.set_ylabel('occurrence freq.', fontsize = 16)
            ax.set_ylim(-0.05,1.05)
    
    axes[0,0].text(-0.25, 1.03, 'A', transform = axes[0,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[0,1].text(-0.25, 1.03, 'B', transform = axes[0,1].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,0].text(-0.25, 1.03, 'C', transform = axes[1,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,1].text(-0.25, 1.03, 'D', transform = axes[1,1].transAxes, fontsize = 25, fontweight = 'bold')

    axes[0,0].text(0.5, 1.15, 'small immigration', transform = axes[0,0].transAxes, ha='center', fontsize = 16, fontweight = 'bold')
    axes[0,1].text(0.5, 1.15, 'large immigration', transform = axes[0,1].transAxes, ha='center', fontsize = 16, fontweight = 'bold')

    # Save

    fig.savefig('fig3.pdf', dpi = 300, bbox_inches='tight')

    

### Figure 4 ###
# Occurrence-abundance pattern in general non-neutral communities.
def fig4():
    
    fig, axes = mp.subplots(nrows = 1, ncols = 3, figsize = (11.1, 3.2))

    N = int(1E3)
    S = 40    

    ## Fig 4A
    
    m = 1E-1 * N

    # Non-neutral points
    
    data = np.load('data/5C.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
        
    axes[0].semilogx(mean_freq_nneutral, P0_nneutral, '.', color = 'black')
    
    # Fig 4B and Fig 4C
        
    S = 6
    p = [0.33443333, 0.00168333, 0.12375, 0.08033333, 0.4176, 0.0422]
    gR = [1.97227108, 1.13511873, 1.0345976, 1.07541209, 0.78716513, 1.31372745]
    dR = [1.00512998, 0.82321306, 1.00487796, 1.00530296, 0.94770017, 0.84279409]
    gr = np.ones(2)
    dr = np.ones(2)
    focal_type = 0
    
    P = np.logspace(-4,0,50)
    
    ## Fig 4B
    
    m = 1E-3 * N
    
    # Neutral curve
    
    P0_neutral = np.zeros(len(P))
    mean_freq_neutral = np.zeros(len(P))
    
    for p_ in P:
        
        prob = recurrence_prob_dist(m, N, [p_, 1 - p_], gr, dr, focal_type)
        
        P0_neutral[p_ == P] = prob_occurrence(prob)
        mean_freq_neutral[p_ == P] = mean_abundance(prob, N) / N
                    
    axes[1].semilogx(mean_freq_neutral, P0_neutral, color = 'black')
    
    # Neutral points

    P0_neutral = np.zeros(S)
    mean_freq_neutral = np.zeros(S)

    for i in range(S):
                
        p_ = [p[i], 1 - p[i]]
    
        prob = recurrence_prob_dist(m, N, p_, gr, dr, 0)
        
        P0_neutral[i] = prob_occurrence(prob)
        mean_freq_neutral[i] = mean_abundance(prob, N) / N
                
    axes[1].semilogx(mean_freq_neutral, P0_neutral, '.', color = '#4778F7')

    # Non-neutral points
        
    P0_nneutral = np.zeros(S)
    mean_freq_nneutral = np.zeros(S)
    
    for i in range(S):
                                        
        P_ni_nneutral, n_j_nneutral = reduction_prob_dist(m, N, p, gR, dR, S, i)

        P0_nneutral[i] = prob_occurrence(P_ni_nneutral)
        mean_freq_nneutral[i] = mean_abundance(P_ni_nneutral, N) / N
                
    axes[1].semilogx(mean_freq_nneutral, P0_nneutral, '.', color = 'black')
    
    # Assumption from data
    
    P0_assumed = np.zeros(S)
    mean_freq_assumed = np.zeros(S)

    for i in range(S):
        
        p_ = [mean_freq_nneutral[i], 1 - mean_freq_nneutral[i]]
    
        prob = recurrence_prob_dist(m, N, p_, gr, dr, 0)
        
        P0_assumed[i] = prob_occurrence(prob)
        mean_freq_assumed[i] = mean_abundance(prob, N) / N
        
    axes[1].semilogx(mean_freq_assumed, P0_assumed, '.', color = '#FFB000')
    
    # Differences btw points
    
    for i in range(S):
        
        axes[1].annotate('', xy = [mean_freq_nneutral[i], P0_nneutral[i]], xytext = [mean_freq_neutral[i], P0_neutral[i]], color = '#4778F7', arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color = '#4778F7'))
        axes[1].annotate('', xy = [mean_freq_nneutral[i], P0_nneutral[i]], xytext = [mean_freq_assumed[i], P0_assumed[i]], color = '#FFB000', arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color = '#FFB000'))  
        
    ss_assumed = sum(np.sqrt((mean_freq_nneutral - mean_freq_assumed)**2 + (P0_nneutral - P0_assumed)**2)**2)
    ss_real = sum(np.sqrt((mean_freq_nneutral - mean_freq_neutral)**2 + (P0_nneutral - P0_neutral)**2)**2)
    
    axes[1].text(0.32, 0.8, 'sum of sq. = %.3f'%ss_assumed, color = '#FFB000', horizontalalignment='center', verticalalignment='center', transform = axes[1].transAxes)
    axes[1].text(0.32, 0.9, 'sum of sq. = %.3f'%ss_real, color = '#4778F7', horizontalalignment='center', verticalalignment='center', transform = axes[1].transAxes)
    
    ## Fig 4C
    
    m = 1E-1 * N
        
    # Neutral curve
    
    P0_neutral = np.zeros(len(P))
    mean_freq_neutral = np.zeros(len(P))
    
    for p_ in P:
    
        prob = recurrence_prob_dist(m, N, [p_, 1 - p_], gr, dr, focal_type)
        
        P0_neutral[p_ == P] = prob_occurrence(prob)
        mean_freq_neutral[p_ == P] = mean_abundance(prob, N) / N
            
    axes[2].semilogx(mean_freq_neutral, P0_neutral, color = 'black')
    
    # Neutral points

    P0_neutral = np.zeros(S)
    mean_freq_neutral = np.zeros(S)

    for i in range(S):
        
        p_ = [p[i], 1 - p[i]]
    
        prob = recurrence_prob_dist(m, N, p_, gr, dr, 0)
        
        P0_neutral[i] = prob_occurrence(prob)
        mean_freq_neutral[i] = mean_abundance(prob, N) / N
        
    axes[2].semilogx(mean_freq_neutral, P0_neutral, '.', color = '#4778F7')

    # Non-neutral points
        
    P0_nneutral = np.zeros(S)
    mean_freq_nneutral = np.zeros(S)
    
    for i in range(S):
                                   
        P_ni_nneutral, n_j_nneutral = reduction_prob_dist(m, N, p, gR, dR, S, i)

        P0_nneutral[i] = prob_occurrence(P_ni_nneutral)
        mean_freq_nneutral[i] = mean_abundance(P_ni_nneutral, N) / N
                
    axes[2].semilogx(mean_freq_nneutral, P0_nneutral, '.', color = 'black')
    
    # Assumption from data
    
    P0_assumed = np.zeros(S)
    mean_freq_assumed = np.zeros(S)

    for i in range(S):
        
        p_ = [mean_freq_nneutral[i], 1 - mean_freq_nneutral[i]]
    
        prob = recurrence_prob_dist(m, N, p_, gr, dr, 0)
        
        P0_assumed[i] = prob_occurrence(prob)
        mean_freq_assumed[i] = mean_abundance(prob, N) / N
        
    axes[2].semilogx(mean_freq_assumed, P0_assumed, '.', color = '#FFB000')
    
    # Differences btw points
    
    for i in range(S):
          
        axes[2].annotate('', xy = [mean_freq_nneutral[i], P0_nneutral[i]], xytext = [mean_freq_neutral[i], P0_neutral[i]], color = '#4778F7', arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color = '#4778F7'))
        axes[2].annotate('', xy = [mean_freq_nneutral[i], P0_nneutral[i]], xytext = [mean_freq_assumed[i], P0_assumed[i]], color = '#FFB000', arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color = '#FFB000'))
    
    ss_assumed = sum(np.sqrt((mean_freq_nneutral - mean_freq_assumed)**2 + (P0_nneutral - P0_assumed)**2)**2)
    ss_real = sum(np.sqrt((mean_freq_nneutral - mean_freq_neutral)**2 + (P0_nneutral - P0_neutral)**2)**2)
    
    axes[2].text(0.68, 0.1, 'sum of sq. = %.3f'%ss_assumed, color = '#FFB000', horizontalalignment='center', verticalalignment='center', transform = axes[2].transAxes)
    axes[2].text(0.68, 0.2, 'sum of sq. = %.3f'%ss_real, color = '#4778F7', horizontalalignment='center', verticalalignment='center', transform = axes[2].transAxes)
    
    # Annotate

    fig.subplots_adjust(wspace=0.4)

    for ax in [axes[0], axes[1], axes[2]]:
    
        ax.set_xlabel('mean frequency', fontsize = 16)    
        ax.set_ylabel('occurrence freq.', fontsize = 16)
        ax.set_xlim(7E-5,1.5)
        ax.set_ylim(-0.05,1.05)
        ax.tick_params(axis='both', direction='in', which='both', right=True, top=True)
        ax.xaxis.set_major_locator(LogLocator(base=10,numticks=6))
        ax.xaxis.set_minor_locator(LogLocator(base=10,subs=np.arange(0.2,1.,0.1),numticks=6))
        ax.xaxis.set_minor_formatter(NullFormatter())
    
    axes[0].text(-0.25, 1.03, 'A', transform = axes[0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1].text(-0.25, 1.03, 'B', transform = axes[1].transAxes, fontsize = 25, fontweight = 'bold')
    axes[2].text(-0.25, 1.03, 'C', transform = axes[2].transAxes, fontsize = 25, fontweight = 'bold')

    axes[1].text(0.5, 1.1, 'small immigration', transform = axes[1].transAxes, ha='center', fontsize = 18, fontweight = 'bold')
    axes[2].text(0.5, 1.1, 'large immigration', transform = axes[2].transAxes, ha='center', fontsize = 18, fontweight = 'bold')
 
    # Save
    
    mp.savefig('fig4.pdf', dpi = 300, bbox_inches='tight')

    

### Figure 5 ###
# Occurrence-abundance pattern for different levels of asymmetry in the parameters.
def fig5():
    
    fig, axes = mp.subplots(nrows = 2, ncols = 2, figsize = (6.6, 6.6))
            
    N = 1E3
    m = 1E-1 * N
    
    # Neutral curves (5A, 5B, 5C, 5D)

    P = np.logspace(-4,0,50)
    gr = np.ones(2)
    dr = np.ones(2)
    focal_type = 0
    
    P0_neutral = np.zeros(len(P))
    mean_freq_neutral = np.zeros(len(P))
    
    for p_ in P:
        
        prob = recurrence_prob_dist(m, N, [p_, 1 - p_], gr, dr, focal_type)
        
        P0_neutral[p_ == P] = prob_occurrence(prob)
        mean_freq_neutral[p_ == P] = mean_abundance(prob, N) / N
                    
    for ax in [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]:
        
        ax.semilogx(mean_freq_neutral, P0_neutral, color = 'black', alpha = 0.2)
    
    ## Fig 5A
        
    data = np.load('data/5A.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
        
    axes[0,0].scatter(mean_freq_nneutral, P0_nneutral, s = 15, c = np.sqrt(abs(1.-data['gR'])**2 + abs(1.-data['dR'])**2), cmap = 'viridis', vmin = 0)
    
    ## Fig 5B
        
    data = np.load('data/5B.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
        
    axes[0,1].scatter(mean_freq_nneutral, P0_nneutral, s = 15, c = np.sqrt(abs(1.-data['gR'])**2 + abs(1.-data['dR'])**2), cmap = 'viridis', vmin = 0)
    
    ## Fig 5C
        
    data = np.load('data/5C.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
        
    axes[1,0].scatter(mean_freq_nneutral, P0_nneutral, s = 15, c = np.sqrt(abs(1.-data['gR'])**2 + abs(1.-data['dR'])**2), cmap = 'viridis', vmin = 0)
        
    ## Fig 5D
    
    data = np.load('data/5D.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
        
    axes[1,1].scatter(mean_freq_nneutral, P0_nneutral, s = 15, c = np.sqrt(abs(1.-data['gR'])**2 + abs(1.-data['dR'])**2), cmap = 'viridis', vmin = 0)
    
    # Annotate
    
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    
    for ax in [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]:
    
        ax.set_xlabel('mean frequency', fontsize = 16)    
        ax.set_ylabel('occurrence freq.', fontsize = 16)
        ax.set_xlim(7E-5,1.5)
        ax.set_ylim(-0.05,1.05)
        ax.tick_params(axis='both', direction='in', which='both', right=True, top=True)
        ax.xaxis.set_major_locator(LogLocator(base=10,numticks=6))
        ax.xaxis.set_minor_locator(LogLocator(base=10,subs=np.arange(0.2,1.,0.1),numticks=6))
        ax.xaxis.set_minor_formatter(NullFormatter())
    
    axes[0,0].text(-0.25, 1.03, 'A', transform = axes[0,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[0,1].text(-0.25, 1.03, 'B', transform = axes[0,1].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,0].text(-0.25, 1.03, 'C', transform = axes[1,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,1].text(-0.25, 1.03, 'D', transform = axes[1,1].transAxes, fontsize = 25, fontweight = 'bold')

    axes[0,0].text(0.5, 1.15, r"low variance in $f_i$, $\phi_i$", transform = axes[0,0].transAxes, ha='center', fontsize = 16, fontstyle = 'italic')
    axes[0,1].text(0.5, 1.15, r"high variance in $f_i$, $\phi_i$", transform = axes[0,1].transAxes, ha='center', fontsize = 16, fontstyle = 'italic')

    axes[0,0].text(-0.45, 0.5, r"high symmetry in $p_i$", transform = axes[0,0].transAxes, va='center', fontsize = 15.5, fontstyle = 'italic', rotation=90)
    axes[1,0].text(-0.45, 0.5, r"low symmetry in $p_i$", transform = axes[1,0].transAxes, va='center', fontsize = 15.5, fontstyle = 'italic', rotation=90)

    # Save

    fig.savefig('fig5.pdf', dpi = 300, bbox_inches='tight')

    

### Figure 6 ###
# Effect of growth, death, and immigration at the level of types.
def fig6():
    
    fig = mp.figure(figsize = (8.7, 3.2))
    
    gs = gridspec.GridSpec(22, 56, wspace=10, hspace=10)
    
    axes00 = mp.subplot(gs[:20, :20]) 
    axes01 = mp.subplot(gs[:10, 22:32])
    axes02 = mp.subplot(gs[:10, 34:44])  
    axes03 = mp.subplot(gs[:10, 46:])  
    axes11 = mp.subplot(gs[12:, 22:32])
    axes12 = mp.subplot(gs[12:, 34:44])
    axes13 = mp.subplot(gs[12:, 46:])
        
    N = int(1E3)
    m = 1E-1 * N

    # Neutral curve (6A, 6B, 6C, 6D)

    P = np.logspace(-4,0,50)
    gr = np.ones(2)
    dr = np.ones(2)
    focal_type = 0
    
    P0_neutral = np.zeros(len(P))
    mean_freq_neutral = np.zeros(len(P))
    
    for p_ in P:
        
        prob = recurrence_prob_dist(m, N, [p_, 1 - p_], gr, dr, focal_type)
        
        P0_neutral[p_ == P] = prob_occurrence(prob)
        mean_freq_neutral[p_ == P] = mean_abundance(prob, N) / N
                    
    for ax in [axes00, axes01, axes02, axes03, axes11, axes12, axes13]:
        ax.semilogx(mean_freq_neutral, P0_neutral, color = 'black', alpha = 0.2)
    
    ## Fig 6A
        
    data = np.load('data/5C.npz')
    
    P0_nneutral_ref = data['P0_nneutral']
    
    mean_freq_nneutral_ref = data['mean_freq_nneutral']
        
    axes00.semilogx(mean_freq_nneutral_ref, P0_nneutral_ref, '.', color = 'black')
    
    axes00.scatter(mean_freq_nneutral_ref[2], P0_nneutral_ref[2], s=100, facecolors='none', edgecolors='#8E888D', linewidth = 1.5)
    
    axes00.scatter(mean_freq_nneutral_ref[14], P0_nneutral_ref[14], s=100, facecolors='none', edgecolors='#FF8917', linewidth = 1.5)

    ## Fig 6B.1

    data = np.load('data/6B.1.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
    
    axes11.scatter(mean_freq_nneutral, P0_nneutral, s=20, c = data['p_change'], cmap = 'viridis', norm=colors.LogNorm(vmin=1E-4, vmax=1))
    
    axes11.scatter(mean_freq_nneutral_ref[2], P0_nneutral_ref[2], s=80, facecolors='none', edgecolors='#8E888D', linewidth = 1.5)
        
    ## Fig 6B.2
    
    data = np.load('data/6B.2.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
    
    axes01.scatter(mean_freq_nneutral, P0_nneutral, s=20, c = data['p_change'], cmap = 'viridis', norm=colors.LogNorm(vmin=1E-4, vmax=1))
    
    axes01.scatter(mean_freq_nneutral_ref[14], P0_nneutral_ref[14], s=80, facecolors='none', edgecolors='#FF8917', linewidth = 1.5)
    
    ## Fig 6C.1
    
    data = np.load('data/6C.1.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
    
    axes12.scatter(mean_freq_nneutral, P0_nneutral, s=20, c = data['gR_change'], cmap = 'coolwarm')
    
    axes12.scatter(mean_freq_nneutral_ref[2], P0_nneutral_ref[2], s=80, facecolors='none', edgecolors='#8E888D', linewidth = 1.5)
    
    ## Fig 6C.2
    
    data = np.load('data/6C.2.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
    
    axes02.scatter(mean_freq_nneutral, P0_nneutral, s=20, c = data['gR_change'], cmap = 'coolwarm')
    
    axes02.scatter(mean_freq_nneutral_ref[14], P0_nneutral_ref[14], s=80, facecolors='none', edgecolors='#FF8917', linewidth = 1.5)
        
    ## Fig 6D.1
    
    data = np.load('data/6D.1.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
    
    axes13.scatter(mean_freq_nneutral, P0_nneutral, s=20, c = data['dR_change'], cmap = 'coolwarm_r')
    
    axes13.scatter(mean_freq_nneutral_ref[2], P0_nneutral_ref[2], s=80, facecolors='none', edgecolors='#8E888D', linewidth = 1.5)
        
    ## Fig 6D.2
    
    data = np.load('data/6D.2.npz')
    
    P0_nneutral = data['P0_nneutral']
    
    mean_freq_nneutral = data['mean_freq_nneutral']
    
    axes03.scatter(mean_freq_nneutral, P0_nneutral, s=20, c = data['dR_change'], cmap = 'coolwarm_r')
    
    axes03.scatter(mean_freq_nneutral_ref[14], P0_nneutral_ref[14], s=80, facecolors='none', edgecolors='#FF8917', linewidth = 1.5)

    # Annotate
    
    # Add colorbars
    cax11 = fig.add_axes([0.435, 0., 0.12, 0.02])
    cax12 = fig.add_axes([0.605, 0., 0.12, 0.02])
    cax13 = fig.add_axes([0.775, 0., 0.12, 0.02])
    cax11.axis('off')
    cax12.axis('off')
    cax13.axis('off')
    cbar11 = cax11.imshow(np.array([[1E-4, 1]]), norm=colors.LogNorm(vmin=1E-4, vmax=1), cmap='viridis')
    cbar12 = cax12.imshow(np.array([[0.8, 1.2]]), cmap='coolwarm')
    cbar13 = cax13.imshow(np.array([[0.8, 1.2]]), cmap='coolwarm_r')
    cbar11.set_visible(False)
    cbar12.set_visible(False)
    cbar13.set_visible(False)
    cbar11 = fig.colorbar(cbar11, ax = cax11, orientation = 'horizontal', fraction = 1.)
    cbar12 = fig.colorbar(cbar12, ax = cax12, orientation = 'horizontal', fraction = 1.)
    cbar13 = fig.colorbar(cbar13, ax = cax13, orientation = 'horizontal', fraction = 1.)
    
    # Add triangles indicating the values on panel 6A
    p = data['p']
    cax11.annotate(r'$\blacktriangledown$', xy = (30 - 13.5 * np.log10(p[2]*1E4),0), color = '#8E888D', fontsize = 14, xycoords = ("axes fraction", "axes fraction"))
    cax11.annotate(r'$\blacktriangledown$', xy = (30 - 13.5 * np.log10(p[14]*1E4),0), color = '#FF8917', fontsize = 14, xycoords = ("axes fraction", "axes fraction"))
    gR = data['gR']
    cax12.annotate(r'$\blacktriangledown$', xy = (30 - 54 * ((gR[2]-0.8)/0.4),0), color = '#8E888D', fontsize = 14, xycoords = ("axes fraction", "axes fraction"))
    cax12.annotate(r'$\blacktriangledown$', xy = (30 - 54 * ((gR[14]-0.8)/0.4),0), color = '#FF8917', fontsize = 14, xycoords = ("axes fraction", "axes fraction"))
    dR = data['dR']
    cax13.annotate(r'$\blacktriangledown$', xy = (30 - 54 * ((dR[2]-0.8)/0.4),0), color = '#8E888D', fontsize = 14, xycoords = ("axes fraction", "axes fraction"))
    cax13.annotate(r'$\blacktriangledown$', xy = (30 - 54 * ((dR[14]-0.8)/0.4),0), color = '#FF8917', fontsize = 14, xycoords = ("axes fraction", "axes fraction"))

    for ax in [axes00, axes01, axes02, axes03, axes11, axes12, axes13]:
        ax.set_xlim(7E-5,1.5)
        ax.set_ylim(-0.05,1.05)    
        ax.tick_params(axis='both', direction='in', which='both', right=True, top=True)
        ax.xaxis.set_major_locator(LogLocator(base=10,numticks=3))
        
    for ax in [axes01, axes02, axes03]:
        ax.tick_params(color='#FF8917', labelcolor='#FF8917')
        mp.setp(ax.spines.values(), color='#FF8917')
        
    for ax in [axes11, axes12, axes13]:
        ax.tick_params(color='#8E888D', labelcolor='#8E888D')
        mp.setp(ax.spines.values(), color='#8E888D')

    axes00.set_xlabel('mean frequency', fontsize = 16)    
    axes00.set_ylabel('occurrence freq.', fontsize = 16)
    axes00.xaxis.set_major_locator(LogLocator(base=10,numticks=6))
    axes00.xaxis.set_minor_locator(LogLocator(base=10,subs=np.arange(0.2,1.,0.1),numticks=6))
    axes00.xaxis.set_minor_formatter(NullFormatter())
    
    axes00.text(-0.25, 1.13, 'A', transform = axes00.transAxes, fontsize = 25, fontweight = 'bold')
    axes01.text(-0.25, 1.3, 'B', transform = axes01.transAxes, fontsize = 25, fontweight = 'bold')
    axes02.text(-0.25, 1.3, 'C', transform = axes02.transAxes, fontsize = 25, fontweight = 'bold')
    axes03.text(-0.25, 1.3, 'D', transform = axes03.transAxes, fontsize = 25, fontweight = 'bold')
    
    axes01.text(0.5, 1.1, 'fr. migrants', transform = axes01.transAxes, ha='center', fontsize = 12, fontweight = 'bold')
    axes02.text(0.5, 1.1, 'growth r.', transform = axes02.transAxes, ha='center', fontsize = 12, fontweight = 'bold')
    axes03.text(0.5, 1.1, 'death r.', transform = axes03.transAxes, ha='center', fontsize = 12, fontweight = 'bold')
 
    # Save
    
    mp.savefig('fig6.pdf', dpi = 300, bbox_inches='tight')
    


### Figure 7 ###
# Variance and second raw moment of the frequency.     
def fig7():
    
    fig, axes = mp.subplots(nrows = 2, ncols = 2, figsize = (6.6, 6.6))
    
    focal_type = 0
    
    N = 1E3
    m = 1E-1 * N
    
    abundance = np.arange(int(N + 1))
    P = np.logspace(-4,0,50)
    
    ## Fig 7A and 7C (differences in growth rate)
    
    dr = [1., 1.]

    # Growth rates and colors and labels associated to them    
    Gr_i = np.array([0.8, 0.9, 1., 1.1, 1.2])
    color_gr = np.array(['lightseagreen', 'darkorchid', 'black', 'crimson', 'darkorange'])
    label_gr = np.array(['0.8', '0.9', 'neutral', '1.1', '1.2'])

    for gr_i in Gr_i:
    
        second_central_moment = np.ones(len(P))*np.nan
        second_raw_moment = np.ones(len(P))*np.nan
        
        for p_ in P:
        
            p = [p_, 1. - p_]
            gr = [gr_i, 1.]
            
            prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    
            mean_freq = sum(prob * abundance/N)  
            second_central_moment[p_ == P] = sum(prob * (np.arange(N+1)/N - mean_freq)**2)
            second_raw_moment[p_ == P] = sum(prob * (abundance/N)**2)
            
        second_central_moment[second_central_moment == 0] = np.nan
        second_raw_moment[second_raw_moment == 0] = np.nan

        axes[0,0].loglog(P, second_central_moment, '-', color = color_gr[gr_i == Gr_i][0], label = label_gr[gr_i == Gr_i][0])
        axes[0,1].loglog(P, second_raw_moment, '-', color = color_gr[gr_i == Gr_i][0], label = label_gr[gr_i == Gr_i][0])
        
    ## Fig 7B and 7D (differences in death rate)
        
    gr = [1., 1.]

    # Death rates and colors and labels associated to them
    Dr_i = np.array([0.8, 0.9, 1., 1.1, 1.2])
    color_dr = np.array(['darkorange', 'crimson', 'black', 'darkorchid', 'lightseagreen'])
    label_dr = np.array(['0.8', '0.9', 'neutral', '1.1', '1.2'])
    
    for dr_i in Dr_i:
    
        second_central_moment = np.ones(len(P))*np.nan
        second_raw_moment = np.ones(len(P))*np.nan
        
        for p_ in P:
        
            p = [p_, 1. - p_]
            dr = [dr_i, 1.]
            
            prob = recurrence_prob_dist(m, N, p, gr, dr, focal_type)
    
            mean_freq = sum(prob * abundance/N)  
            second_central_moment[p_ == P] = sum(prob * (np.arange(N+1)/N - mean_freq)**2)
            second_raw_moment[p_ == P] = sum(prob * (abundance/N)**2)
            
        second_central_moment[second_central_moment == 0] = np.nan
        second_raw_moment[second_raw_moment == 0] = np.nan

        axes[1,0].loglog(P, second_central_moment, '-', color = color_dr[dr_i == Dr_i][0], label = label_dr[dr_i == Dr_i][0])
        axes[1,1].loglog(P, second_raw_moment, '-', color = color_dr[dr_i == Dr_i][0], label = label_dr[dr_i == Dr_i][0])

    # Annotate

    axes[0,1].legend(title=r'$f_1$', fontsize = 10, bbox_to_anchor=(1.05, 0.5), loc='center left')
    axes[1,1].legend(title=r'$\phi_1$', fontsize = 10, bbox_to_anchor=(1.05, 0.5), loc='center left')
    
    fig.subplots_adjust(wspace=0.4, hspace=0.4)
    
    for ax in [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]:
    
        ax.set_xlabel('fraction of migrants', fontsize = 16)
        ax.set_xlim(7E-4,1.5)
        ax.set_ylim(1E-6,1.5)        
        ax.tick_params(axis='both', direction='in', which='both', right=True, top=True)
        
        if ax == axes[0,0] or ax == axes[1,0]: ax.set_ylabel('frequency variance', fontsize = 16)
        else: ax.set_ylabel(r'freq. 2$^{nd}$ raw mom.', fontsize = 16)
    
    axes[0,0].text(-0.25, 1.03, 'A', transform = axes[0,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[0,1].text(-0.25, 1.03, 'B', transform = axes[0,1].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,0].text(-0.25, 1.03, 'C', transform = axes[1,0].transAxes, fontsize = 25, fontweight = 'bold')
    axes[1,1].text(-0.25, 1.03, 'D', transform = axes[1,1].transAxes, fontsize = 25, fontweight = 'bold')

    # Save

    fig.savefig('fig7.pdf', dpi = 300, bbox_inches='tight')