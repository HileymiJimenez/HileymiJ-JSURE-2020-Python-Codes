# -*- coding:utf-8 -*-
"""
Spydervediction

this is a tempora script file.
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt
def model (ym, t):
   k = -1.1
   dydt = k * ym
   return dydt

# initial codition
y0 = 300000


# time points
t = np.linspace(0,4)

# solve ODE
ym = odeint(model,y0,t)

# plot data
td = [0,1,2,3,4]
yd = [300000, 73125, 28125, 632500, 588125]
plt.plot(td, yd, '*')

# plot results
plt.plot(t,ym)
plt.xlabel("time")
plt.ylabel("y(t)")
plt.show()
    
