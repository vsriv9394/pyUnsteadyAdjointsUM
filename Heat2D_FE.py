import numpy as np
import matplotlib as mpl
mpl.rcParams.update({'font.size':16, 'figure.figsize':(10,8)})
#mpl.rc('text', usetex=True)
import matplotlib.pyplot as plt
import ctypes as C
import sys

class LaplaceEquation:

    # INITIALIZATION FUNCTION =====================================================================================

    def __init__(self,\
                 num_points=33,\
                 fourier_number=0.25,\
                 conductivity=1.0,\
                 heat_transfer_coefficient=1.0,\
                 beta=None,\
                 init_states=None,\
                 boundary_target=None,\
                 n_steps=100,\
                 verbose=False):
    
        # INITIALIZE THE C LIBRARY FUNCTIONS ------------------------------------------------------------------
        
        self.lib = C.CDLL('./2DHeatEqn_ForwardEuler_src/lib2DHeatEquationForwardEuler.so')
        
        self._Update                    = self.lib.__getattr__('Update')
        self._AdjointUpdate             = self.lib.__getattr__('AdjointUpdate')
        self._Objective                 = self.lib.__getattr__('Objective')
        self._AdjointObjective          = self.lib.__getattr__('AdjointObjective')
        
        self._Update.restype            = None
        self._AdjointUpdate.restype     = None
        self._Objective.restype         = None
        self._AdjointObjective.restype  = None
        
        self._Update.argtypes           = [C.c_long,\
                                           C.c_long,\
                                           C.c_double,\
                                           C.c_double,\
                                           C.c_double,\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double)]
        
        self._AdjointUpdate.argtypes    = [C.c_long,\
                                           C.c_long,\
                                           C.c_double,\
                                           C.c_double,\
                                           C.c_double,\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double)]
        
        self._Objective.argtypes        = [C.c_long,\
                                           C.c_long,\
                                           C.c_double,\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double)]
        
        self._AdjointObjective.argtypes = [C.c_long,\
                                           C.c_long,\
                                           C.c_double,\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double),\
                                           C.POINTER(C.c_double)]
        
        # INITIALIZE THE PARAMETERS FOR THE PROBLEM -----------------------------------------------------------
        
        self.npoints = num_points
        self.F0      = fourier_number
        self.alpha   = conductivity
        self.alphab  = heat_transfer_coefficient
        self.nsteps  = n_steps
        self.h       = 2./(self.npoints-1)
        self.dt      = self.F0*self.h**2/self.alpha
        self.obj     = np.array([0.0])
        self.verbose = verbose
        
        if boundary_target is not None:
            self.ub  = boundary_target
        else:
            self.ub  = np.ones((self.npoints))
        
        # INTIALIZE THE COORDINATES HERE (FOR PLOTTING PURPOSES) ----------------------------------------------
        
        self.coord   = np.linspace(-1.,1.,self.npoints)
        
        self.x = np.zeros((self.npoints, self.npoints))
        self.y = np.zeros((self.npoints, self.npoints))
        
        for i in range(self.npoints):
            self.x[i,:] = self.coord
            self.y[i,:] = self.coord[i]
            self.x[:,i] = self.coord[i]
            self.y[:,i] = self.coord
        
        # INITIALIZE THE VARIABLES FOR DIRECT AND ADJOINT SOLUTION --------------------------------------------
        
        self.beta = np.ones((self.npoints * (self.nsteps+1)))
        
        if beta is not None:
            if np.shape(self.beta)[0]==np.shape(beta)[0] and len(np.shape(beta))==1:
                self.beta = beta
            else:
                print("Wrong dimensions for beta... Exiting\n")
                sys.exit()
        
        self.sens = np.zeros((self.npoints * (self.nsteps+1)))
        
        self.u    = np.zeros((self.npoints * self.npoints * (self.nsteps+1)))
        
        if init_states is not None:
            if np.shape(init_states)[0]==self.npoints*self.npoints and len(np.shape(init_states))==1:
                self.u[0:self.npoints*self.npoints] = init_states
            else:
                print("Wrong dimensions for initial states... Exiting\n")
                sys.exit()
        
        self.psi  = np.zeros((self.npoints * self.npoints * (self.nsteps+1)))
        
        self.history = []

    # GET FUNCTIONS ===============================================================================================    
    
    def GetFieldToPlotAtTime(self, istep):
        return np.reshape(self.u[self.npoints*self.npoints*istep:self.npoints*self.npoints*(istep+1)],
                          [self.npoints, self.npoints],
                          order='C')
    
    def GetAdjFieldToPlotAtTime(self, istep):
        return np.reshape(self.psi[self.npoints*self.npoints*istep:self.npoints*self.npoints*(istep+1)],
                          [self.npoints, self.npoints],
                          order='C')
    
    # EVALUATION FUNCTIONS ========================================================================================
    
    def DirectSolve(self):
        for i in range(1,self.nsteps+1):
            self._Update(self.npoints,
                         i,
                         self.F0,
                         self.h,
                         self.alphab,
                         self.ub.ctypes.data_as(C.POINTER(C.c_double)),
                         self.u.ctypes.data_as(C.POINTER(C.c_double)),
                         self.beta.ctypes.data_as(C.POINTER(C.c_double)))
            self.history.append(np.linalg.norm(self.GetFieldToPlotAtTime(i)))
            if self.verbose:
                sys.stdout.flush()
                sys.stdout.write("\r Iteration %5d with Field RMS value %.15E"%(i, self.history[-1]))
        if self.verbose:
            sys.stdout.write("\n")
    
    def Objective(self):
        self._Objective(self.npoints,
                        self.nsteps,
                        self.dt,
                        self.u.ctypes.data_as(C.POINTER(C.c_double)),
                        self.beta.ctypes.data_as(C.POINTER(C.c_double)),
                        self.obj.ctypes.data_as(C.POINTER(C.c_double)))
    
    def AdjointObjective(self):
        for i in range(self.nsteps,0,-1):
            self._AdjointObjective(self.npoints,
                                   i,
                                   self.dt,
                                   self.u.ctypes.data_as(C.POINTER(C.c_double)),
                                   self.beta.ctypes.data_as(C.POINTER(C.c_double)),
                                   self.psi.ctypes.data_as(C.POINTER(C.c_double)),
                                   self.sens.ctypes.data_as(C.POINTER(C.c_double)))
    
    def AdjointSolve(self):
        for i in range(self.nsteps,0,-1):
            self._AdjointUpdate(self.npoints,
                                i,
                                self.F0,
                                self.h,
                                self.alphab,
                                self.ub.ctypes.data_as(C.POINTER(C.c_double)),
                                self.u.ctypes.data_as(C.POINTER(C.c_double)),
                                self.beta.ctypes.data_as(C.POINTER(C.c_double)),
                                self.psi.ctypes.data_as(C.POINTER(C.c_double)),
                                self.sens.ctypes.data_as(C.POINTER(C.c_double)))
            if self.verbose:
                sys.stdout.flush()
                sys.stdout.write("\r Iteration %5d with Field RMS value %20.15E"%(i, np.linalg.norm(self.GetAdjFieldToPlotAtTime(i))))
        if self.verbose:
            sys.stdout.write("\n")
    
    # PLOTTING FUNCTIONS ==========================================================================================
    
    def PlotHist(self):
        
        plt.plot(self.history, '-b')
        plt.xlabel("Time step")
        plt.ylabel("RMS value over field")
        plt.title("Simulation History")
    
    def PlotAdjSol(self, istep):
        
        plt.contourf(self.x, self.y, self.GetAdjFieldToPlotAtTime(istep), 50)
        plt.title("Sensitivities at time step %d"%(istep))
        plt.colorbar()
    
    def PlotSol(self, istep):
        
        plt.contourf(self.x, self.y, self.GetFieldToPlotAtTime(istep), 50)
        plt.title("System state at time step %d"%(istep))
        plt.colorbar()

    def PlotSensorTarget(self):

        t = np.zeros((self.nsteps+1))
        for i in range(self.nsteps):
            t[i+1] = t[i] + self.dt
        u_match = 2.2*(np.exp(-0.1*t**2)-1.0)
        plt.plot(u_match, '-k', label='Target')
        
    
    def PlotSensor(self, x_id, y_id, style, label):
        
        plt.plot(self.u[(self.npoints-1)*self.npoints+int((self.npoints+1)/2)::self.npoints*self.npoints], style, label=label)
        plt.xlabel("Time step")
        plt.ylabel("Sensor measurement")
        plt.title("Sensor data")
        plt.legend()
    
    def PlotShow(self):
        
        plt.show()

# MAIN FUNCTION =======================================================================================================

if __name__=="__main__":

    npoints = 65

    Eqn = LaplaceEquation(num_points=npoints,\
                          n_steps=30000,\
                          boundary_target=10*np.cos(1.5*np.pi*np.linspace(-1.,1.,npoints)),\
                          heat_transfer_coefficient=2.0,\
                          fourier_number=0.2)

    print("============================================================")
    print(" Direct Solution")
    print("-----------------")

    Eqn.DirectSolve()
    Eqn.Objective()

    print(" Objective Function : %lf"%(Eqn.obj))
    print("========================================")
    print(" Adjoint Solution")
    print("------------------")

    Eqn.AdjointObjective()
    Eqn.AdjointSolve()

    print(" Sensitivities calculated using Adjoints")
    print("============================================================")
    
    Eqn.PlotHist()
    Eqn.PlotShow()
    
    Eqn.PlotSensorTarget()
    Eqn.PlotSensor(17,-1,'-b','Baseline')
    Eqn.PlotShow()
    
    Eqn.PlotSol(Eqn.nsteps-1)
    Eqn.PlotShow()

    Eqn.PlotAdjSol(Eqn.nsteps-19999)
    Eqn.PlotShow()
