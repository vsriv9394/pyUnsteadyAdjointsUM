import numpy as np
from Heat2D_FE import *

npoints = 65
nsteps  = 30000

start_opt       = 200
end_opt         = 300
opt_step_length = 0.5

if start_opt==0:
    beta = np.ones((npoints*(nsteps+1)))
else:
    beta = np.loadtxt("beta_%04d.dat"%(start_opt))

Eqn = LaplaceEquation(num_points=npoints,\
                      n_steps=nsteps,\
                      boundary_target=10*np.cos(1.5*np.pi*np.linspace(-1.,1.,npoints)),\
                      heat_transfer_coefficient=2.0,\
                      fourier_number=0.2,\
                      beta=beta)

print("============================================================")
print("          FIELD INVERSION USING UNSTEADY ADJOINTS")
print("============================================================")
print("|      Iteration          |             Obj. Fn.           |")
print("------------------------------------------------------------")

for i in range(start_opt, end_opt):
    Eqn.DirectSolve()
    Eqn.Objective()
    Eqn.AdjointObjective()
    Eqn.AdjointSolve()
    print("|      %9d          |      %.15E     |"%(i, Eqn.obj))
    Eqn.beta = Eqn.beta - opt_step_length * Eqn.sens/np.linalg.norm(Eqn.sens)
    if (i+1)%20==0:
        np.savetxt("beta_%04d.dat"%(i+1), Eqn.beta)

print("------------------------------------------------------------")
