import numpy as np
from Heat2D_FE import *

npoints  = 65
nsteps   = 30000
opt_iter = 300

Eqn = LaplaceEquation(num_points=npoints,\
                      n_steps=nsteps,\
                      boundary_target=10*np.cos(1.5*np.pi*np.linspace(-1.,1.,npoints)),\
                      heat_transfer_coefficient=2.0,\
                      fourier_number=0.2,\
                      verbose=True)

Eqn.PlotSensorTarget()
Eqn.DirectSolve()
Eqn.PlotSensor(17,-1,'-b','Baseline')

Eqn.beta = np.loadtxt("beta_files/beta_%04d.dat"%(opt_iter))

Eqn.DirectSolve()
Eqn.PlotSensor(17,-1,'-r','Inversion')
Eqn.PlotShow()

Eqn.PlotSol(29999)
Eqn.PlotShow()
