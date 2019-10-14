#include "2DHeat_ForwEuler.h"
#include <math.h>

void AdjointObjective(int n,
	              int id,
	              double dt,
	              double *all_state,
	              double *all_beta,
	              double *all_adj_state,
	              double *all_adj_beta){

	double u_match;

	int i, j;

	double *this_adj_state = all_adj_state + n*n*id;
	double *this_adj_beta  = all_adj_beta  + n*id;

	for(j=0; j<n; j++){
		this_adj_beta[j] = 0.0;
		for(i=0; i<n; i++){
			this_adj_state[j*n+i] = 0.0;
		}
	}
	
	u_match        = 2.2 * (exp(-0.1*id*dt*id*dt)-1.0);
	this_adj_state[(n-1)*n + (n+1)/2] = 2 * (all_state[id*n*n + (n-1)*n + (n+1)/2] - u_match);

}
