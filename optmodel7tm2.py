# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 10:09:29 2020

@author: Jimenez
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
# from scipy.interpolate import interp1d

# Import CSV data file
# Column 1 = time (t)
# Column 2 = input (u)
# Column 3 = output (yp)
data = np.loadtxt('data5um.txt',delimiter=',')
# u0 = data[0,1]
# yp0 = data[0,1]
# t = data[:,0].T - data[0,0]
t = data[:,0]
# u = data[:,1].T
yp = data[:,1]

# specify number of steps
ns = len(t)
# delta_t = t[1]-t[0]
# create linear interpolation of the u data versus time
# uf = interp1d(t,u)

# define first-order plus dead-time approximation    
# def fopdt(y,t,uf,Km,taum,thetam)
def fopdt7(y,t,x,Km1,Km2,Km3,Km4,beta):
    # arguments
    #  y      = output
    #  t      = time
    #  uf     = input linear function (for time shift)
    #  Km     = model gain
    #  taum   = model time constant
    #  thetam = model time constant
    # time-shift u]
    Km4= x[3]
    beta = x[4]
  #  try:
    if t <= 2:
                Km1 = x[0]
                dydt = Km1*y**(2/3)-beta*y-Km4*y
    else: 
        if t <= 4:
                Km2 = x[1]
                dydt = Km2*y**(2/3)-beta*y-Km4*y
        else:
                Km3 = x[2]
                dydt = Km3*y**(2/3)-beta*y-Km4*y
  #  except:
    #    print('Error with time extrapolation: ' + str(t))
    #    um = 0
    # calculate derivative
    #    dydt = Km*y                # (-(y-yp0))*Km # + Km * (um-u0))/taum
    return dydt

# simulate FOPDT model with x=[Km,taum,thetam]
def sim_model(x):
    # input arguments
    Km1 = x[0]
    Km2 = x[1]
    Km3 = x[2]
    Km4 = x[3]
    beta = x[4]
    # taum = x[1]
    # thetam = x[2]
    # storage for model values
    ym = np.zeros(ns)  # model
    # initial condition
    ym[0] = 300000
    # loop through time steps    
    for i in range(0,ns-1):
        ts = [t[i],t[i+1]]
  #              y1 = odeint(fopdt,ym[i],ts,args=(uf,Km,taum,thetam))
        y1 = odeint(fopdt7,ym[i],ts,args=(x,Km1,Km2,Km3,Km4,beta))
        ym[i+1] = y1[-1]
    return ym

# define objective
def objective(x):
    # simulate model
    ym = sim_model(x)
    # calculate objective
    obj = 0.0
    for i in range(len(ym)):
        obj = obj + (ym[i]-yp[i])**2    
    # return result
    return obj

# initial guesses
x0 = np.zeros(5)
x0[0] = -20.5 # Km1
x0[1] = 20.5 # Km2 --> taum
x0[2] = 45.0 # Km3
x0[3] = 0.25
x0[4] = 0.999 # beta

# show initial objective
print('Initial SSE Objective: ' + str(objective(x0)))
print('Kp01: ' + str(x0[0]),', Kp02: ' + str(x0[1]), ' and Kp03: ' + str(x0[2]))
print('Kp04: ' + str(x0[3]))
print(' beta0: ' + str(x0[4]))
# optimize Km, taum, thetam
#solution = minimize(objective,x0)
solution = minimize(objective,x0,method='powell',options={'xtol': 1e-8, 'maxfev': 100000, 'disp': True})


# Another way to solve: with bounds on variables
#bnds = ((0.4, 0.6), (1.0, 10.0), (0.0, 30.0))
#solution = minimize(objective,x0,bounds=bnds,method='SLSQP')
x = solution.x

# show final objective
print('Final SSE Objective: ' + str(objective(x)))

print('Kp1: ' + str(x[0]),', Kp2: ' + str(x[1]),' and Kp3: ' + str(x[2]))
print('Kp4: ' + str(x[3]))
print(' beta: ' + str(x[4]))

# print('taup: ' + str(x[1]))
# print('thetap: ' + str(x[2]))

# calculate model with updated parameters
ym1 = sim_model(x0)
ym2 = sim_model(x)
# plot results
# plt.figure(1)
# plt.subplot(2,1,1)
plt.plot(t,yp,'ko-',linewidth=2,label='Experiment Data')
plt.plot(t,ym1,'b-',linewidth=2,label='Initial Guess')
plt.plot(t,ym2,'r--',linewidth=3,label='Optimized Model')
plt.xlabel('Days')
plt.ylabel('Number of Cells')
plt.xlabel('Day')
plt.title('Best Fit Model - Bertalanffy')
plt.legend(loc='best')
#plt.subplot(2,1,2)
#plt.plot(t,x[0],'bx-',linewidth=2)
#plt.plot(t,x[1],'r--',linewidth=3)
# plt.legend(['Measured','Interpolated'],loc='best')
# plt.ylabel('Input Data')
plt.savefig('outputmodel73Kmbeta.png')
plt.show()
