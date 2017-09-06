#!/usr/bin/env python

''' analytical derivatives check'''

from sympy import *

# define symbols
phi_p = Symbol('phi_p')
phi_s = Symbol('phi_s')
N = Symbol('N')
alpha = Symbol('alpha')
sigma = Symbol('sigma')
chi12 = Symbol('chi12')

# free energy
f = (phi_p/N)*log(phi_p/2.0) + phi_s*log(phi_s/2.0) + (1-phi_p-phi_s)*log(1-phi_p-phi_s) - alpha*(phi_p*sigma + phi_s)**1.5 + 0.5*chi12*phi_p*phi_p

print 'Free Energy function: \n'
print f
print '\n'

df_dphi_p = diff(f,phi_p)
#print 'df/dphi_p\n',df_dphi_p

# for critical conditions
print 'Critical point calculation: \n'
print 'Second derivative: \n'
df_dphi_p2 = diff(f,phi_p,2)
print df_dphi_p2
print '\n'
print 'Third derivative: \n'
df_dphi_p3 = diff(f,phi_p,3)
print df_dphi_p3
print '\n'
print 'For Jacobian \n'
print 'df2/dphi_p'
print diff(df_dphi_p2,phi_p)
print '\n'
print diff(df_dphi_p2,phi_s)
print '\n'
print diff(df_dphi_p3,phi_p)
print '\n'
print diff(df_dphi_p3,phi_s)

print 'For phase diagram: \n'
print "f'(phi_1a) - f'(phi_2a)\n"

phi_p1 = Symbol('phi_p1')
phi_p2 = Symbol('phi_p2')

f_phi_p1 = (phi_p1/N)*log(phi_p1/2.0) + phi_s*log(phi_s/2.0) + (1-phi_p1-phi_s)*log(1-phi_p1-phi_s) - alpha*(phi_p1*sigma + phi_s)**1.5 + 0.5*chi12*phi_p1*phi_p1
f_phi_p2 = (phi_p2/N)*log(phi_p2/2.0) + phi_s*log(phi_s/2.0) + (1-phi_p2-phi_s)*log(1-phi_p2-phi_s) - alpha*(phi_p2*sigma + phi_s)**1.5 + 0.5*chi12*phi_p2*phi_p2

F1=diff(f_phi_p1,phi_p1)-diff(f_phi_p2,phi_p2)
print F1
print '\n'

print "(b-a)*f'(phi_1a) -[ f(phi_2a) - f(phi_1a) ]\n"
F2 = (phi_p2-phi_p1)*diff(f_phi_p1,phi_p1)-(f_phi_p2-f_phi_p1)

print '\n'
print F2
print '\n'

# For jacobian
print 'For Jacobian \n'
#dF1/dphi_p2
print 'df1/dphi_p2\n'
print diff(F1,phi_p2)
print '\n'
#dF1/dphi_s
print 'df1/dphis\n'
print diff(F1,phi_s)
print '\n'
#dF2/dphi_p2
print 'df2/dphi_p2'
print '\n'
print diff(F2,phi_p2)
print '\n'
#dF2/dphi_s
print 'df2/dphi_s'
print '\n'
print diff(F2,phi_s)
