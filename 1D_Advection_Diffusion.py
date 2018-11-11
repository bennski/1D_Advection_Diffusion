# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:36:59 2016

@author: bennski

-----------------------------------------------------------------------------
 Backward time, centered space (BTCS) scheme paired with the Upwind method for 
 solving 1D advection-diffusion equation for a local nucleobase concentration 
 in a warm little pond.

 The equation solved is

       du     d  du     du
       -- = D -- -- - U --
       dt     dx dx     dx

 D:             diffusion coefficient
 U:             convective fluid velocity
 x = 0 = xmax:  bottom of pond
 x = xmax/2:    top of pond
-----------------------------------------------------------------------------
"""

from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np
import matplotlib.patches as mpatches
import numpy.matlib

styles = mpatches.ArrowStyle.get_styles()

def shift(l, n):
    return np.append(l[n:],l[:n])
    
# Initialize writer
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

# Initialization function for animation
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

# Set the diffusion of mass coefficient
D = 4e-10 # [m^2s^-1]
T = 273.15+65

xmax = 2*1. # in metres
tmax = 2000. # in seconds
level = 9

# The average fluid velocity resulting from a temperature gradient in the pond
delta_T = 0.01 # temperature differance from base of WLP to surface
beta = 564e-6 # volumetric thermal expansion coefficient of water at 65 degrees C
g = 9.81 # gravity [m s^-2]
U = np.sqrt(g*beta*delta_T*xmax/2.)

nx = (2**level) + 1
nt = (2**(level+1)) + 1 #Choosing nt to have twice as many grid points as nx

# Concentration array
u = np.matlib.zeros(shape=(nt, nx))

# Arrays for plotting
x = np.linspace(0,xmax,nx)
t = np.linspace(0,tmax,nt)

# Calculate delta_x, delta_t
delta_x = x[1] - x[0]
delta_t = t[2] - t[1]

#Initial condition
maxconc = 234.
sum_init = 0
for k in xrange(0,nx):
    u[0,k] = maxconc*np.e**(-(k-int(nx/2.))**2/(2*1**2)) # Gaussian
    sum_init = sum_init + u[0,k] # For mass frac & conservation of mass
    
# For plotting initial condition
r1 = np.zeros(shape=(nx))
for r in xrange(0,nx):
    r1[r] = u[0,r]*nx/sum_init

# Arrays for making system of equations
s = np.matlib.zeros(shape=(nx, 1))
c = np.matlib.zeros(shape=(nx, nx))

#Set up animation 
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(0, xmax), ylim=(0, 0.05))
line, = ax.plot([], [], '-', color='#550000', lw=2)
time_template = '%.1f minutes'
time_text = ax.text(0.01, 0.96, '', ha='left', va='center', transform=plt.gca().transAxes, fontsize=14, weight='bold', color='black')
    
ax.text(0.68, 0.955, r'u = %.1f cms$^{-1}$' % (U*100), ha='left', va='center', transform=plt.gca().transAxes, fontsize=14, weight='bold', color='black')
ax.arrow(1.55, 0.043, 0.2, 0.0, head_width=0.004, head_length=0.07, fc='k', ec='k')
    
# System of equations solver of finite difference approximation of schrodinger equation
def animate(n):	
    if n == -1:
        line.set_data(x,r1)
        time_text.set_text(time_template % (0))
    else:
        # Input Cyclic Boundary Conditions at x = 0 and x = nx-1
        c[0,0] = 1 + D*delta_t/delta_x**2
        s[0] = u[n,0]
        c[0,1] = U*delta_t/delta_x
        c[0,nx-1] = -D*delta_t/delta_x**2 - U*delta_t/delta_x
    
        c[nx-1,nx-1] = 1 + D*delta_t/delta_x**2 + U*delta_t/delta_x
        s[nx-1] = u[n,nx-1]
        c[nx-1,nx-2] = -U*delta_t/delta_x
        c[nx-1,0] = -D*delta_t/delta_x**2
        
        for j in xrange(1,nx-1):
            # Terms on RHS
            s[j] = u[n,j]
            # Coefficients on LHS with varying potential depending on position
            # 1st-order Upwind Method:
            c[j,j-1] = -D*delta_t/(delta_x**2) - U*delta_t/delta_x
            c[j,j] = (1 + 2*D*delta_t/(delta_x**2) + U*delta_t/delta_x)
            c[j,j+1] = -D*delta_t/(delta_x**2)

        # Compute the values of the solution using matrix division (i.e. cu = s -> u = c^-1*s)			
        sol = np.linalg.inv(c)*s
        # Input solutions into the solution array at time n+1 
        for k in xrange(0,nx):
            u[n+1, k] = sol[k]
            
        # Check conservation of mass
        Cons2 = 0
        for h in xrange(0,nx):
            Cons2 = Cons2 + sol[h]
        if np.abs(Cons2 - sum_init)/sum_init > 0.01:
            print 'ERROR: Conservation of mass was broken. %.3f %% of mass lost.'%(100*np.abs(Cons2 - sum_init)/sum_init)
            
        # Array to shift results into the convective fluid frame
        movingFrame = np.zeros(shape=(nx))
        for m in xrange(0,nx):
            movingFrame[m] = sol[m]
            
        line.set_data(x, np.transpose(shift(movingFrame,int((n+1)*U*delta_t/delta_x)%nx))/sum_init)
        time_text.set_text(time_template % (((n+1)*delta_t)/(60.))) # change units to minutes
    return line, time_text

        
ani = animation.FuncAnimation(fig, animate, np.arange(-1,nt-1),
                              interval=1, save_count=25, blit=True, init_func=init)
                           
                           
plt.xticks([0,0.5,1.0,1.5],[1.0,1.5,0.0,0.5],fontsize=14)
plt.yticks(fontsize=14)
plt.title(r'Mixing of Local Concentration in a Pond',fontsize=14)
plt.ylabel('Fraction of total initial concentration',fontsize=14)
plt.xlabel(r'Convection cell position, moving frame (m)',fontsize=14)
plt.tight_layout()

ani.save('1D_Advection_Diffusion.mp4', dpi = 300, writer=writer)
