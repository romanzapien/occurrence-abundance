# occurrence-abundance

This repository includes the code used to produce the results presented 
in *Zapien, et al.*, 2021.

The aim is to study the effect of growth and death rates on the 
occurrence-abundance pattern of communities.

## Organization

1. Data (**occurrence-abundance/data**).

Directory storing empirical data and to store simulated data.

2. Figures (**occurrence-abundance/figures**).

Scripts to produce all the figures.

3. Numerics (**occurrence-abundance/numerics**).

Scripts to solve the model numerically using the master equation.

4. Simulations (**occurrence-abundance/simulation**).

Scripts to simulate communities using the Gillespie algorithm.

## Usage

Scripts are structured in the following way:

`sc_*.py` contains the source code.

`par_*.py` contains the parameters.

`exe_*.py` produces the output.

`make_*.py` generates data.

## Requirements

All the code has been developed in Python 3.6.
