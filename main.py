#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import VO
import argparse

parser = argparse.ArgumentParser(description='calculate VO phase diagram')
parser.add_argument('--outfile', type=argparse.FileType('w'), help="output file to be created",default='outfile.dat')
parser.add_argument('--N', type=int, help="polymer length",default=1000)
parser.add_argument('--sigma', type=float, help="polymer charge density",default=0.44)
parser.add_argument('--alpha', type=float, help="reduced temperature",default=3.655)
parser.add_argument('--chi12', type=float, help="chi parameter between polymers", default=0.0)
args = parser.parse_args()

N=args.N
sigma=args.sigma
alpha = args.alpha
chi12 = args.chi12

print 'Free Energy Function: \n'
print "phi_p*log(0.5*phi_p)/N + phi_s*log(0.5*phi_s) + (1 - phi_p - phi_s)*log(1 - phi_p - phi_s) -alpha*(phi_p*sigma + phi_s)**1.5 + 0.5*chi12*phi_p**2\n" 
print 'Polymer length: ',N
print 'Polymer charge density: ',sigma
print 'alpha: ',alpha
print 'chi12: ',chi12

initial_guess = [0.05,0.0]
print 'initial guess: ',initial_guess

z = VO.criticalpoint_solve(alpha,N,sigma,chi12)
print 'critical point: ', z

phi_p,phi_s= VO.NewtonRaphson_solve(alpha,N,sigma,chi12,initial_guess)

for i in range(len(phi_p)):
    args.outfile.write(str(phi_p[i][0])+' '+str(phi_s[i][0])+'\n')

args.outfile.close()

plt.plot(phi_p,phi_s,'o-')
plt.xlabel('phi_p')
plt.ylabel('phi_s')
plt.show()
