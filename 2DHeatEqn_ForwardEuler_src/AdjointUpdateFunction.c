#include "2DHeat_ForwEuler.h"

void AdjointUpdate(int n, 
		   int id,
	           double F0,
		   double h,
	           double alphab,
	           double *state_b,
		   double *all_state,
		   double *all_beta,
	           double *all_adj_state, 
	           double *all_adj_beta){

	int i, j;

	double *last_state = all_state + n*n*(id-1);
	double *this_beta  = all_beta  + n*id;

	double *this_adj_state = all_adj_state + n*n*id;
	double *last_adj_state = all_adj_state + n*n*(id-1);
	double *this_adj_beta  = all_adj_beta  + n*id;

	for(j=0; j<n; j++){
		last_adj_state[j*n+1] += 4./3.*this_adj_state[j*n];
		last_adj_state[j*n+2] -= 1./3.*this_adj_state[j*n];
		last_adj_state[(j+1)*n-2] += 4./3.*this_adj_state[(j+1)*n-1];
		last_adj_state[(j+1)*n-3] -= 1./3.*this_adj_state[(j+1)*n-1];
	}

	for(i=0; i<n; i++){
		last_adj_state[  n+i] += 4./3.*this_adj_state[i];
		last_adj_state[2*n+i] -= 1./3.*this_adj_state[i];
		last_adj_state[i] -= 2./3. * h * this_beta[i] * alphab * this_adj_state[i];
		this_adj_beta[i] += 2./3. * h * alphab * (state_b[i] - last_state[i]) * this_adj_state[i];
		last_adj_state[(n-2)*n+i] += 4./3.*this_adj_state[(n-1)*n+i];
		last_adj_state[(n-3)*n+i] -= 1./3.*this_adj_state[(n-1)*n+i];
	}	

	for(i=1; i<n-1; i++)
		for(j=1; j<n-1; j++){
			last_adj_state[(j+1)*n+i] += F0 * this_adj_state[j*n+i];
			last_adj_state[(j-1)*n+i] += F0 * this_adj_state[j*n+i];
			last_adj_state[j*n+(i+1)] += F0 * this_adj_state[j*n+i];
			last_adj_state[j*n+(i-1)] += F0 * this_adj_state[j*n+i];
			last_adj_state[j*n+i] += (1.-4.*F0) * this_adj_state[j*n+i];
	}

}
