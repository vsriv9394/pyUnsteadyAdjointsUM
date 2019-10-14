import numpy as np
import matplotlib as mpl
mpl.rcParams.update({'font.size':16, 'figure.figsize':[10,8]})
import matplotlib.pyplot as plt

npoints  = 65
nsteps   = 30000
opt_iter = 300

beta = np.loadtxt("beta_%04d.dat"%(opt_iter))

x,t = np.meshgrid(np.linspace(-1.,1.,65), np.linspace(0.,nsteps,nsteps+1))
plt.contourf(x, t, np.reshape(beta, [30001, 65], order='C'), 50)
plt.xlabel("X-direction")
plt.ylabel("Time")
plt.title("Field Inversion (Unsteady)")
plt.colorbar()
plt.show()
