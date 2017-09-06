#!/usr/bin/env python

#Import necessary modules
from math import *
import numpy as np
from numpy.linalg import inv

def critical_conditions(x,sigma,alpha,N,chi12):
    phi_p = x[0]
    phi_s = x[1]

    # critical point from 2nd and third derivative
    df2 = -0.75*alpha*sigma**2.*(phi_p*sigma + phi_s)**(-0.5) + 1.0*chi12 - 1.0/(phi_p + phi_s - 1.0) + 1.0/(N*phi_p)
    df3 = 0.375*alpha*sigma**3.*(phi_p*sigma + phi_s)**(-1.5) + (phi_p + phi_s - 1.0)**(-2.) - 1.0/(N*phi_p**2.)
    return np.array([df2,df3])

def calculate_criticalpoint_jacobian(x,sigma,alpha,N,chi12):
    ''' df2/dphip,df2/dphis,df3/dphi,df3/dphis '''
    phi_p = x[0]
    phi_s = x[1]
    f1=0.375*alpha*sigma**3*(phi_p*sigma + phi_s)**(-1.5) + (phi_p + phi_s - 1)**(-2) - 1.0/(N*phi_p**2)
    f2=0.375*alpha*sigma**2*(phi_p*sigma + phi_s)**(-1.5) + (phi_p + phi_s - 1)**(-2)
    f3=-0.5625*alpha*sigma**4*(phi_p*sigma + phi_s)**(-2.5) - 2/(phi_p + phi_s - 1)**3 + 2.0/(N*phi_p**3)
    f4=-0.5625*alpha*sigma**3*(phi_p*sigma + phi_s)**(-2.5) - 2/(phi_p + phi_s - 1)**3
    return np.array([ [f1,f2],[f3,f4] ])


def criticalpoint_solve(alpha,N,sigma,chi12):
    """ Newton Raphson solver for the binary mixture"""
    R = [0,0] #phi_p, phi_s
    # initial guess
    new_R = [0.01,0.01] #phi_p, phi_s
    iter = 0
    max_iter = 2000
    # tolerance
    tolerance = 1e-5
    # step size
    step_size = 0.1
    while iter < max_iter :
        iter += 1
        R = new_R
        jacobian = calculate_criticalpoint_jacobian(R,sigma,alpha,N,chi12)
        invjac = inv(jacobian)
        f_r = critical_conditions(R,sigma,alpha,N,chi12)
        new_R = R - step_size*np.dot(invjac,f_r)
        if abs(new_R[0] - R[0]) < tolerance and abs(new_R[1]-R[1]) < tolerance:
            phi_p=R[0]
            phi_s=R[1]
            break

    return (phi_p,phi_s)

# check
#def calculate_critical_point(sigma,alpha,N,chi12):
#
#    x = np.zeros((2))
#    # initial guess
#    x[0] = 1e-5
#    x[1] = 1e-5
#
#    x_crit = fsolve(critical_conditions, x, args = (sigma, alpha, N, chi12))
#
#    return x_crit

def common_tangent_fun(x,phi_p1,sigma,alpha,N,chi12):
    "F1 = f'(phi_1a) - f'(phi_2a); F2 = (b-a)*f'(phi_1a) -[ f(phi_2a) - f(phi_1a) ]"

    F1 = -1.5*alpha*sigma*(phi_p1*sigma + x[1])**0.5 + 1.5*alpha*sigma*(x[0]*sigma + x[1])**0.5 + 1.0*chi12*phi_p1 - 1.0*chi12*x[0] - log(-phi_p1 - x[1] + 1.0) + log(-x[0] - x[1] + 1) + log(0.5*phi_p1)/N - log(0.5*x[0])/N

    F2 = -alpha*(phi_p1*sigma + x[1])**1.5 + alpha*(x[0]*sigma + x[1])**1.5 + 0.5*chi12*phi_p1**2 - 0.5*chi12*x[0]**2 + (-phi_p1 + x[0])*(-1.5*alpha*sigma*(phi_p1*sigma + x[1])**0.5 + 1.0*chi12*phi_p1 - log(-phi_p1 - x[1] + 1) - 1 + log(0.5*phi_p1)/N + 1.0/N) + (-phi_p1 - x[1] + 1)*log(-phi_p1 - x[1] + 1) - (-x[0] - x[1] + 1)*log(-x[0] - x[1] + 1) + phi_p1*log(0.5*phi_p1)/N - x[0]*log(0.5*x[0])/N

    return np.array([F1,F2])

def calculate_jacobian(x,phi_p1,sigma,alpha,N,chi12):
    "dF1/dphi_p2, dF1/dphi_s; dF2/dphi_p2, dF2/dphi_s"
    df1dphi = 0.75*alpha*sigma**2*(x[0]*sigma + x[1])**(-0.5) - 1.0*chi12 - 1./(-x[0] - x[1] + 1.) - 1.0/(N*x[0])

    df1dpsi = -0.75*alpha*sigma*(phi_p1*sigma + x[1])**(-0.5) + 0.75*alpha*sigma*(x[0]*sigma + x[1])**(-0.5) - 1.0/(-x[0] - x[1] + 1.) + 1.0/(-phi_p1 - x[1] + 1.)

    df2dphi = -1.5*alpha*sigma*(phi_p1*sigma + x[1])**0.5 + 1.5*alpha*sigma*(x[0]*sigma + x[1])**0.5 + 1.0*chi12*phi_p1 - 1.0*chi12*x[0] - log(-phi_p1 - x[1] + 1.) + log(-x[0] - x[1] + 1.) + log(0.5*phi_p1)/N - log(0.5*x[0])/N

    df2dpsi = -1.5*alpha*(phi_p1*sigma + x[1])**0.5 + 1.5*alpha*(x[0]*sigma + x[1])**0.5 + (-phi_p1 + x[0])*(-0.75*alpha*sigma*(phi_p1*sigma + x[1])**(-0.5) + 1./(-phi_p1 - x[1] + 1.0)) - log(-phi_p1 - x[1] + 1.) + log(-x[0] - x[1] + 1.)
    # Jacobian matrix
    return np.array([ [df1dphi,df1dpsi],[df2dphi,df2dpsi] ])

def NewtonRaphson_solve(alpha,N,sigma,chi12,guess):
   """ Newton Raphson solver for the binary mixture"""
   # get critical point
   critical_phi = criticalpoint_solve(alpha,N,sigma,chi12)
   # create points, phi_p1 in range up to critical point
   delta_phi = 0.0005 # point spacing
   phi_p_vals = np.arange(1e-5,critical_phi[0],delta_phi)
   phi_p_vals = phi_p_vals.tolist()
   # allocate roots
   R = [0,0] #phi_p2, phi_s
   # initial guess
   new_R = guess #phi_p2, phi_s
   iter = 0
   # allocate arrays
   phi_s = np.zeros((len(phi_p_vals),1))
   phi_p2 = np.zeros((len(phi_p_vals),1))
   phi_p1 = np.zeros((len(phi_p_vals),1))
   # maximum number of iterations
   max_iter = 2000
   # tolerance
   tolerance = 1e-5
   # step size
   step_size = 0.1
   #Loop to find the np.roots using Multivariate Newton-Rhapson
   for phi in phi_p_vals:
      iter = 0
      while iter < max_iter :
         iter += 1
         index = phi_p_vals.index(phi)
         R = new_R
         jacobian = calculate_jacobian(R,phi,sigma,alpha,N,chi12)
         invjac = inv(jacobian)
         f_r = common_tangent_fun(R,phi,sigma,alpha,N,chi12)
         new_R = R - step_size*np.dot(invjac,f_r)
         if abs(new_R[0] - R[0]) < tolerance and abs(new_R[1]-R[1]) < tolerance:
            phi_p1[index] = phi
            phi_p2[index] = new_R[0]
            phi_s[index] = new_R[1]
            break

   # build data list from phi_p1 and phi_p2 arrays
   phi_p1=phi_p1.tolist()
   phi_p2=phi_p2.tolist()
   phi_p2=phi_p2[::-1] #reverse the order of phi_p2
   phi_s=phi_s.tolist()
   phi_s_rev = phi_s[::-1]
   #put data together into one array
   phi_p = phi_p1 + phi_p2
   phi_s = phi_s + phi_s_rev

   return (phi_p,phi_s)
