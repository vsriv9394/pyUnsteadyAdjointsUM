#include "2DHeat_ForwEuler.h"
#include <math.h>

void Objective(int n,
	       int nsteps,
	       double dt,
	       double *all_state,
	       double *all_beta,
	       double *objective){

	int i, j, k;
	double u_match;
	*objective = 0.0;
	
	for(k=0; k<nsteps; k++){
	
		u_match     = 2.2*(exp(-0.1*k*dt*k*dt)-1.0);
		*objective += pow(all_state[k*n*n + (n-1)*n + (n+1)/2] - u_match, 2);
	
	}	

}
