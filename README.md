# Mistletoe–host tree–bird model

MATLAB codes associated with the manuscript
“Dynamics of a mistletoe–host tree–bird system: parasitism, mutualism and multistability”.

## Files

- `Tree_Mistletoe_bird_model.m`  
  Generates the numerical simulations for the three parameter values
  \(d_1=0.6582\), \(d_2=0.6555\), and \(d_3=0.6654\), and computes the
  equilibria and their local stability.

- `Monodromy_Floquet.m`  
  Computes the monodromy matrix and Floquet multipliers for the stable
  periodic orbit at \(d_2=0.6555\).

- `Mistletoe_PRCC.m`  
  Performs the Latin Hypercube Sampling and PRCC sensitivity analysis
  conditioned on the oscillatory regime.

## Software

The numerical simulations were performed in MATLAB.
Numerical continuation and bifurcation analysis were carried out using MATCONT.

## Notes

Minor graphical formatting of the figures included in the manuscript was adjusted
interactively in MATLAB.