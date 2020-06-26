# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 18:47:07 2020

@author: Jimenez
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt
def model(y, t, k):
    dydt = -k * y
    return dydt

# initial condition 
y0 = 4

# time points
t = np.linspace(0, 10)

# solve ODEs
k = 0.5
y1 = odeint (model, y0, t, args=(k,))
k = 0.8
y2 = odeint (model, y0, t, args=(k,))
k = 0.3
y3 = odeint (model, y0, t, args=(k,))

# plot results
plt.plot(t,y1,'r-',linewidth=2,label='k=0.1')
plt.plot(t,y2,'b--', linewidth=2, label='k=0.2')
plt.plot(t,y3,'g:', linewidth=2, label='k=0.5')
plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()
