import numpy as np
import os

class Boundary:
    def __init__(self, boundary_type, boundary_value):
        self.DefineBoundary(boundary_type, boundary_value)

    def DefineBoundary(self, boundary_type, boundary_value):
        self.type = boundary_type
        # either 'D' for Dirichlet or 'N' for Neumann
        if boundary_type == 'D' or boundary_type == 'N': pass
        else: pass
        
        self.value = boundary_value

class Space:
    def __init__(self):
        pass

    def CreateMesh(self, rowpts, colpts):
        #Setting number of gridpoints by rows and columns
        self.rowpts = rowpts
        self.colpts = colpts

        #Setting up velocity matrices
        self.u = np.ones((self.rowpts+2, self.colpts+2))
        self.v = np.ones((self.rowpts+2, self.colpts+2))
        self.u_star = np.zeros((self.rowpts+2, self.colpts+2))
        self.v_star = np.zeros((self.rowpts+2, self.colpts+2))
        self.u_next = np.zeros((self.rowpts+2, self.colpts+2))
        self.v_next = np.zeros((self.rowpts+2, self.colpts+2))
        self.u_c = np.zeros((self.rowpts,self.colpts))
        self.v_c = np.zeros((self.rowpts,self.colpts))

        #Setting up pressure matrices
        self.p = np.ones((self.rowpts+2, self.colpts+2))
        self.p_c = np.zeros((self.rowpts,self.colpts))

        #Setting up default source term
        self.SetSourceTerm()

    def SetDeltas(self, breadth, length):
        self.dx = length/(self.colpts-1)
        self.dy = breadth/(self.rowpts-1)

    def SetInitialU(self, U):
        self.u = U*self.u

    def SetInitialV(self, V):
        self.v = V*self.v

    def SetInitialP(self, P):
        self.p = P*self.p

    def SetSourceTerm(self, S_x = 0, S_y = 0):
        self.S_x = S_x
        self.S_y = S_y

class Fluid:
    def __init__(self, rho, mu):
        self.SetFluidProperties(rho,mu)

    def SetFluidProperties(self, rho, mu):
        self.rho = rho
        self.mu = mu






