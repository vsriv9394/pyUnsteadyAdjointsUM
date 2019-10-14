#include "2DHeat_ForwEuler.h"

void Update(int n,
	    int id,
	    double F0,
	    double h,
	    double alphab,
	    double *state_b,
	    double *all_state,
	    double *all_beta){

	int i, j;

	double *this_state = all_state + n*n*id;
	double *last_state = all_state + n*n*(id-1);
	double *this_beta  = all_beta  + n*id;

	for(j=1; j<n-1; j++)
		for(i=1; i<n-1; i++)
			this_state[j*n+i] = F0 * (last_state[(j+1)*n+i] +
					          last_state[(j-1)*n+i] +
						  last_state[j*n+(i+1)] +
						  last_state[j*n+(i-1)] +
				     (1./F0-4.) * last_state[j*n+i]);

	for(j=0; j<n; j++){
		this_state[j*n] = (4.*last_state[j*n+1] - last_state[j*n+2])/3.;
		this_state[(j+1)*n-1] = (4.*last_state[(j+1)*n-2] - last_state[(j+1)*n-3])/3.;
	}

	for(i=0; i<n; i++){
		this_state[i] = (4.*last_state[n+i] - last_state[2*n+i] + 2.*h*this_beta[i]*alphab*(state_b[i] - last_state[i]))/3.;
		this_state[(n-1)*n+i] = (4.*last_state[(n-2)*n+i] - last_state[(n-3)*n+i])/3.;
	}	

}
