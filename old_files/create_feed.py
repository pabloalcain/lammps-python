#!/usr/bin/python
#*-* coding: utf-8 *-*
import os
from numpy import linspace
dens=linspace(0.014,0.015,2)
lda=[20]
pot=[1]
npart=[4000]
l=0
N=len(dens)*len(lda)*len(pot)*len(npart)

for i in pot:
    for j in lda:
        for mm in npart:
  	    for k in dens:
                l=l+1
                f=open('feed_temp','w')
                str="potential %i        #Tipo de potencial: 1:Med, 2:Stiff, 3:Horo\n\
coulomb %i          #Rango de Coulomb (lambda enteros entre 0 y 20)\n\
data None #Archivo de data (None para random)\n\
x 0.5              #x\n\
N %i               #Numero de particulas\n\
dens %g            #Densidad [1/fm^3]\n\
temp 1.600,0.1,60  #Temperatura inicial, final, nro pasos\n\
tempstyle gauss    #Estilo de temperatura\n\
Nsteps 50000       #Cantidad de pasos por temperatura \n\
Ndump 1000         #Cantidad de pasos entre dumps"%(i,j,mm,k)
                print>>f, str
                f.close()
                print("Corriendo condicion %i de %i, pot=%i, lda=%i, dens=%g, N=%i" %(l,N,i,j,k,mm))
                os.system("python simple.py feed_temp")
