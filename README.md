# Example of optimisation for control of HPAI

This analysis uses a toy model of HPAI transmission to optimise control interventions: preventive culling and vaccination coverage.

 ## Model
 The model is an SIRV model implemented with SimInf. It is kept purposefully simple to test the optimisation package and provide an example for future work. 
It includes two control strategies:
 * Vaccination of susceptible individuals
 * Preventive culling, implemented at random

The aim is to optimise the decision variables related to these two interventions, vaccination coverage and proportion of preventive culling.

## Objective function
We test two objective functions to minimise:
 * The total number of infected individuals at the end of the epidemic $n_{infected}$
 * The costs of the epidemic, defined as $Cost = C_{vaccination} n_{vaccinated} + C_{PrevCulling} n_{culled} + C_{epidemic} n_{infected}$

## Optimization
We use the package `nloptr` to optimise the control of HPAI. It takes as input the simulation model, the objective function calculated using the result of the simulation for given values of the decision variables and the constraints on the decision variables. It returns the optimal values of the decision variables. 


