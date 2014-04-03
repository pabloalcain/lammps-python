#!/usr/bin/python
from numpy import *
from pylab import *
rc=5.4
rc_coul=60.0

N=5000

V0=373.118
u0=1.5

Va=2666.647
ua=1.6

Vr=3088.118
ur=1.7468

f=open('pandha_m.table','w')


##Medium for neutrons
A=linspace(0,rc,N+1)[1:]

V=V0*exp(-u0*A)/A
#V=V-V[N-1]
F=V0*exp(-u0*A)/A**2 *(u0*A+1)


print>>f,'# Pandha medium potential for same species\n'
print>>f,'NN'
print>>f,'N %i\n' % N

for i in range(N):
    print>>f , i+1,A[i],V[i],F[i]

##Medium for neutron-proton
A=linspace(0,rc,N+1)[1:]

V=Vr*exp(-ur*A)/A-Va*exp(-ua*A)/A
#V=V-V[N-1]
F=Vr*exp(-ur*A)/A**2 *(ur*A+1)-Va*exp(-ua*A)/A**2 *(ua*A+1)

print>>f,'# Pandha medium potential for different species\n'
print>>f,'NP'
print>>f,'N %i\n' % N

for i in range(N):
    print>>f , i+1,A[i],V[i],F[i]


#Medium for proton-proton    


ldas=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,40,50,100]
for i in ldas:
    print i
    rc_coul=max(3*i,rc)
    N_coul=int(N*rc_coul/rc)
    N_cou=5000
    Vc=1.44
    uc=1.0/i
    A=linspace(0,rc_coul,N_coul+1)[1:]
    
    V=V0*exp(-u0*A)/A+Vc*exp(-uc*A)/A
    #V=V-V[N_coul-1]
    F=V0*exp(-u0*A)/A**2 *(u0*A+1)+Vc*exp(-uc*A)/A**2 *(uc*A+1)

    print>>f,'# Pandha medium potential for protons (with coulomb, lambda=%g)\n'%i
    print>>f,'PP%i'%int(i)
    print>>f,'N %i\n' % N_coul
    
    for i in range(N_coul):
        print>>f , i+1,A[i],V[i],F[i]
