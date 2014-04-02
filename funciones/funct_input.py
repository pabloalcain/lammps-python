#*-* coding: utf-8 *-*
import os
import random

def create_input(infile):
    lines=open(infile,'r').readlines()
    for line in lines:
        uncom=line.split("#")[0]
        ide=uncom.split(" ")[0]
        val=uncom.split(" ")[1]
        if ide=="potential":
            pot=get_potential(val)
            continue
        if ide=="coulomb":
            cou=get_coulomb(val)
            continue
        if ide=="data":
            fname=get_data(val)
            continue
        if ide=="x":
            x=get_x(val)
            continue
        if ide=="N":
            N=get_N(val)
            continue
        if ide=="dens":
            rho=get_dens(val)
            continue
        if ide=="temp":
            [Ti,Tf,NTemp]=get_temp(val)
            continue
        if ide=="tempstyle":
            sty=get_tempstyle(val)
            continue
        if ide=="Nsteps":
            Nsteps=get_nsteps(val)
            continue
        if ide=="Ndump":
            Ndump=get_ndump(val)
            continue
        raise ValueError ("Error, opción %s no válida" %ide)
    inp=write_input(pot,cou,fname,x,N,rho,Ti)
    fpath=set_filepath(pot,cou,x,N,rho)
    return inp,N,Ti,Tf,NTemp,sty,Nsteps,Ndump,fpath

def set_filepath(pot,cou,x,N,rho):
    fpath="data/"
    if pot==1:
	pot_path="med/"
    elif pot==2:
	pot_path="stiff/"
    elif pot==3:
	pot_path="horo/"

    cou_path="cou%g/"%cou
    x_path="x"+str(x)+"/"
    N_path="N"+str(N)+"/"
    rho_path="d"+str(rho)+"/"
    fpath=fpath+pot_path+cou_path+x_path+N_path+rho_path

    if not os.path.exists(fpath):
	os.makedirs(fpath)
    return fpath     
def get_potential(val):
    if val =="1" or val=="med":
        pot=1
    elif val=="2" or val=="stiff":
        pot=2
    elif val=="3" or val=="horo":
        pot=3
    else:
        raise NameError ("Error, opción %s no reconocida" % val)
    return pot

def get_coulomb(val):
    try:
	cou=int(val)
    except ValueError:
	raise ValueError ("Error, %s no se puede pasar a entero" % val)
    if (cou<0):# or cou > 20):
	raise ValueError ("Sólo soportados enteros entre 0 y 20 para lambda")
    return cou

def get_data(val):
    if val!="None" and val!="Lattice":
        if not (os.path.exists(val)):
            raise NameError ("Error, %s no existe" % val)
    return val

def get_x(val):
    try:
        x=float(val)
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a float" %val)
    if (x>1 or x<0):
        raise ValueError ("Error, x debe estar entre 0 y 1")
    return x


def get_N(val):
    try:
        N=int(val)
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a entero" %val)
    if (N<=0):
        raise ValueError ("Error, N debe ser mayor a 0")
    return N

def get_dens(val):
    try:
        rho=float(val)
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a float" %val)
    if (rho<=0):
        raise ValueError ("Error, rho debe ser mayor a 0")
    return rho

def get_temp(val):
    try:
        Ti=float(val.split(',')[0])
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a float" %val.split(',')[0])
    if (Ti<=0):
        raise ValueError ("Error, Ti debe ser mayor a 0")

    try:
        Tf=float(val.split(',')[1])
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a float" %val.split(',')[1])
    if (Tf<=0):
        raise ValueError ("Error, Tf debe ser mayor a 0")
    
    try:
        NTemp=int(val.split(',')[2])
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a int" %val.split(',')[2])
    if (NTemp<=0):
        raise ValueError ("Error, NTemp debe ser mayor a 0")
    return [Ti,Tf,NTemp]

def get_tempstyle(val):
    if val =="1" or val=="lin":
        sty=1
    elif val=="2" or val=="exp":
        sty=2
    elif val=="3" or val=="gauss":
        sty=3
    else:
        raise NameError ("Error, opción %s no reconocida" % val)
    return sty


def get_nsteps(val):
    try:
        Nsteps=int(val)
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a entero" %val)
    if (Nsteps<=0):
        raise ValueError ("Error, Nsteps debe ser mayor a 0")
    return Nsteps

def get_ndump(val):
    try:
        Ndump=int(val)
    except ValueError:
        raise ValueError ("Error, %s no se puede pasar a entero" %val)
    if (Ndump<=0):
        raise ValueError ("Error, Ndump debe ser mayor a 0")
    return Ndump

def write_input(pot,cou,fname,x,N,rho,Ti):
    if fname=="None":
        L=(float(N)/rho)**(1.0/3.0)
        Nprot=int(x*N)
        Nneut=N-Nprot
        data_i="region\t\tbox block 0 %f 0 %f 0 %f\n\n" % (L,L,L)+\
            "create_box\t2 box\n\n"+\
            "create_atoms\t1 random %i %i box\n" % (Nprot, random.randint(0,10000))+\
            "create_atoms\t2 random %i %i box\n\n" % (Nneut, random.randint(0,10000))+\
            "mass\t\t1 938.0\n"+\
            "mass\t\t2 938.0\n"+\
            "velocity\tall create %f %i\n\n" % (Ti, random.randint(0,10000))

    elif fname=="Lattice":
	L=int((N/2)**(1.0/3.0))+1
	N=L**3*2
	x=0.5
        data_i="lattice\tcustom %f a1 1.0 0.0 0.0 a2 0.0 1.0 0.0 a3 0.0 0.0 1.0 " % rho+\
        "basis 0.0 0.0 0.0 basis 0.5 0.5 0.5\n" +\
	"region\tbox block 0 %d 0 %d 0 %d\n" % (L, L, 2*L) +\
	"create_box\t2 box\n"+\
	"create_atoms\t2 box basis 1 1\n"+\
	"mass\t\t1 938.0\n"+\
	"mass\t\t2 938.0\n\n"+\
	"velocity\t all create %f %i\n\n" %(Ti, random.randint(0,10000))
	pass
    else:	
    	data_i="read_data %s\n\n" % fname+\
        "mass\t\t1 938.0\n"+\
        "mass\t\t2 938.0\n\n"
        
    potential="variable\ttabla string "
    if pot==1:
        potential=potential+"\"tablas/pandha_m.table\"\n\n"
    elif pot==2:
        potential=potential+"\"tablas/pandha_s.table\"\n\n"
    elif pot==3:
        potential=potential+"\"tablas/horo.table\"\n\n"

    if cou==0:
        potential=potential+"pair_coeff\t1 1 ${tabla} NN 5.4\n"
    else:
        potential=potential+"pair_coeff\t1 1 ${tabla} PP%i %g\n"%(cou,max(cou,5.4))
    potential=potential+"pair_coeff\t1 2 ${tabla} NP 5.4\n"+\
                    "pair_coeff\t2 2 ${tabla} NN 5.4\n\n"

    minimize="minimize\t0 1.0 1000 100000\n\n" if (fname=="None") else ""
    str="# Nuclear model\n\n"+\
        "package\tgpu force/neigh 0 1 -1\n\n"+\
        "units\t\tlj\n"+\
        "atom_style\tatomic\n\n"+\
	"timestep\t0.10\n"+\
        data_i+\
        "pair_style\ttable/gpu linear 5000\n\n"+\
        potential+\
        "kspace_style\tnone\n\n"+\
        "neighbor\t0.3 bin\n"+\
        "neigh_modify\tevery 1 delay 0 check yes one 8000 page 80000\n\n"+\
        "thermo_style custom step temp epair etotal press\n"+\
        "thermo\t\t1000\n\n"+\
	"min_style\thftn\n"+\
        minimize
    if cou!=0:
        str=str+"pair_coeff\t1 1 ${tabla} PP%i %g.0\n"%(cou, max(cou,5.4))
    return str
#
#    fout=open("in.nuclear", "w")
#    
#    print>>fout, str
