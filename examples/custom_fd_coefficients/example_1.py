import numpy as np
from devito import Grid, Function, TimeFunction, Eq, Operator, solve
from devito import ConditionalDimension, Constant, first_derivative, second_derivative
from devito import left, right
from math import exp

import matplotlib
import matplotlib.pyplot as plt

from sympy import symbols, IndexedBase, Indexed, Idx

# Grid
Lx = 10
Nx = 11
dx = Lx/(Nx-1)

grid = Grid(shape=(Nx), extent=(Lx))
time = grid.time_dim
t = grid.stepping_dim
x = grid.dimensions

# time stepping parameters
t_end = 1.0
dt = 0.01
ns = int(t_end/dt)+1

# Devito computation
u = TimeFunction(name='u', grid=grid, time_order=4, space_order=4, save=ns)
u.data[:] = 0.0

# Main equations
eq = Eq(u.dt+u.dx)
#eq = Eq(u.dx, u)

print(eq)

e.match(p*q**r)
{p_: 4, q_: x, r_: 2}
>>> (p*q**r).xreplace(e.match(p*q**r))
4*x**2


#stencil = solve(eq, u.forward)

#print(stencil)

## bc's
#bc = [Eq(u[time+1,0], 0.0)]
#bc += [Eq(u[time+1,-1], 0.0)]

## IM mods:
#def emat(eta):
    #E = np.zeros((2,2))
    #E[0,0] = -(1.-eta)*(1.0-2.0*eta)/(1.0+eta)/(1.0+2.0*eta)
    #E[0,1] = -4.0*(1.0-eta)/(1.0+2.0*eta)
    #E[1,0] = -4.0*(2.0-eta)*(1.0-eta)/(1.0+eta)/(1.+2.0*eta)
    #E[1,1] = 3.0*(2.0-eta)*(1.0-2.0*eta)/eta/(1.0+2.0*eta)
    #return E
#E = emat(eta)
#bci = [Eq(u[time+1,201], E[0,0]*u[time+1,199]+E[0,1]*u[time+1,200])]
#bci += [Eq(u[time+1,202], E[1,0]*u[time+1,199]+E[1,1]*u[time+1,200])]
#for j in range(203,Nx+1):
    #bci += [Eq(u[time+1,j], 0.0)]

#x_t = x*x.spacing
#xg = np.linspace(0,Lx,Nx)

## Create the operators
#op = Operator([Eq(u.forward, stencil)]+bc+bci)

#op.apply(time_M=ns-1, dt=dt)

##fig, ax = plt.subplots()
##ax.plot(xg,u.data[0,:])
##ax.plot(xg,u.data[600,:])
##ax.set(xlabel='x', ylabel='u(x,2)')
##ax.grid()
##plt.show()

#fig, ax = plt.subplots()
#ax.plot(xg,u.data[0,:])
#ax.plot(xg,u.data[-1,:])
#ax.set(xlabel='x', ylabel='u(x,t)')
#ax.grid()
#plt.show()

#fig, ax = plt.subplots()
#ax.plot(xg[0:201],abs(u.data[-1,0:201]-u.data[0,0:201]))
#ax.set(xlabel='x', ylabel='|u(x,2)-u(x,0)|')
#ax.grid()
#plt.show()
