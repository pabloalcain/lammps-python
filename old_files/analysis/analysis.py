#!/usr/bin/python
import sys, time, fnmatch
sys.path.append('./src')
import minkowski, lindemann, gofr, mste
import matplotlib
matplotlib.use('Agg')
from numpy import *
from pylab import *
from os import chdir, listdir


def Calc():
    print "Calculando g(r)........."
    calc_gofr()
    print "Calculando lindemann...."
    lind=calc_lind()
    print "Calculando minkowski...."
    [V,DV,S,DS,C,DC,L,DL]=calc_mink()
    print "Calculando mste........."
    calc_mste()
    return lind,V,DV,S,DS,C,DC,L,DL

def calc_gofr():
    gofr.gofr("evol.lammpstrj",200)
    A=loadtxt("gr.dat")
    S=0*A
    r=A[:,0]
    S[:,0]=r
    q=array(range(len(r)))*2*pi/r[-1]
    for i in range(3):
        S[:,i+1]=1-imag(fft((A[:,i+1]-1)*r))*(r[1]-r[0])*4*pi*dens/q
    savetxt('sk.dat',S,header='r, skall, sknn, skpp, sknp')

    figure()
    plot(A[:,0],A[:,1],'.-')
    title("PCF (all)")
    xlabel("r")
    ylabel("g(r)")
    grid()
    savefig("gr_all.pdf")

    figure()
    plot(S[:,0],S[:,1],'.-')
    title("Structure (all)")
    xlabel("k")
    ylabel("S(k)")
    xlim((0,4))
    grid()
    savefig("sk_all.pdf")

    
    figure()
    plot(A[:,0],A[:,2],'.-')
    title("PCF (nn)")
    xlabel("r")
    ylabel("g(r)")
    grid()
    savefig("gr_nn.pdf")

    figure()
    plot(S[:,0],S[:,2],'.-')
    title("Structure (nn)")
    xlabel("k")
    ylabel("S(k)")
    xlim((0,4))
    grid()
    savefig("sk_nn.pdf")

    
    figure()
    plot(A[:,0],A[:,3],'.-')
    title("PCF (pp)")
    xlabel("r")
    ylabel("g(r)")
    grid()
    savefig("gr_pp.pdf")
    
    figure()
    plot(S[:,0],S[:,3],'.-')
    title("Structure (pp)")
    xlabel("k")
    ylabel("S(k)")
    xlim((0,4))
    grid()
    savefig("sk_pp.pdf")


    figure()
    plot(A[:,0],A[:,4],'.-')
    title("PCF (np)")
    xlabel("r")
    ylabel("g(r)")
    grid()
    savefig("gr_np.pdf")  

    figure()
    plot(S[:,0],S[:,4],'.-')
    title("Structure (np)")
    xlabel("k")
    ylabel("S(k)")
    xlim((0,4))
    grid()
    savefig("sk_np.pdf")
    return
    
def calc_lind():
    lind=lindemann.lindemann("evol.lammpstrj")
    return lind

def calc_mink():
    minkowski.minkowski("evol.lammpstrj",1.8,0.8,False)
    f=open("topologia.dat")
    n=0
    for line in f:
        n=n+1
        if n==1: continue
        a=line.split()
    [n, V,DV,S,DS,C,DC,L,DL] =  [float (x) for x in a]
    f.close()
    return V,DV,S,DS,C,DC,L,DL

def calc_mste():
    mste.mste("evol.lammpstrj")
    A=loadtxt("mste_dist.dat")
    figure()
    plot(A[:,0],A[:,1],'r.')
    title("Cluster Distribution")
    xlabel("Cluster mass")
    ylabel("Number of clusters")
    grid()
    savefig("clus.png")
    return
    
def Therm():
    print "Calculando Termo........"
    f=open("evol.log")

    read=False
    n=0
    E=0
    E2=0
    T=0
    T2=0
    P=0
    P2=0
    for line in f:
        if line[:4] == "Loop": 
            break
        if not read and line[:4] == "Step":
            read = True
            continue
        if not read:
            continue
        n+=1
        [aux, temp, energ, aux, pres]=[float(x) for x in line.split()]
        T+=temp
        T2+=temp**2
        E+=energ
        E2+=energ**2
        P+=pres
        P2+=pres**2
    T=T/n
    E=E/n
    P=P/n
    T2=T2/n
    E2=E2/n
    P2=P2/n

    DT= T2-T**2
    DE= E2-E**2
    DP= P2-P**2
    f.close()
    return T,DT,E,DE,P,DP


if len(sys.argv)==2:
    chdir(sys.argv[1])
th=open("therm_data.dat","w")
print >> th, "#lambda, x, N, dens, temp, Dtemp, energ, Denerg, pres, Dpres, lind, vol, Dvol, sur, Dsur, cur, Dcur, eul, Deul"
chdir('med')

for file in sorted(listdir('.')):
    if fnmatch.fnmatch(file,'cou*'):
        chdir(file)
        lda=file[3:]
        for file in sorted(listdir('.')):
            if fnmatch.fnmatch(file,'x*'):
                chdir(file)
                x=file[1:]
                for file in sorted(listdir('.')):
                    if fnmatch.fnmatch(file,'N*'):
                        chdir(file)
                        N=file[1:]
                        for file in sorted(listdir('.')):
                            if fnmatch.fnmatch(file,'d*'): 
                                chdir(file)
                                dens=file[1:]
                                for file in sorted(listdir('.')):
                                    if fnmatch.fnmatch(file,'T0.001'):
                                        chdir(file)
                                        temp=file[1:]
                                        print "lambda=%s, x=%s, N=%s, d=%s, T=%s" % (lda, x, N, dens, temp)
                                        [lind,V,DV,S,DS,C,DC,L,DL]=Calc()
                                        [T,DT,E,DE,P,DP]=Therm()
                                        print >> th, lda, x, N, dens, T,DT,E,DE,P,DP, lind,V,DV,S,DS,C,DC,L,DL
                                        chdir('..')
                                chdir('..')
                        chdir('..')
                chdir('..')
        chdir('..')
th.close()
