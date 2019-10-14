#include <stdio.h>

void Objective(int n,
	       int nsteps, 
	       double dt,
	       double *all_state, 
	       double *all_beta,
	       double *objective);

void AdjointObjective(int n, 
		      int id, 
		      double dt,
		      double *all_state, 
		      double *all_beta, 
		      double *all_adj_state, 
		      double *all_adj_beta);

void Update(int n,
	    int id,
            double F0,
            double h,
            double alphab,
            double *state_b,
            double *all_state,  
            double *all_beta);

void AdjointUpdate(int n,
		   int id,
	           double F0,
		   double h,
	           double alphab,
	           double *state_b,
		   double *all_state,
		   double *all_beta,
	           double *all_adj_state, 
	           double *all_adj_beta);
