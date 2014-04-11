subroutine Init_Files
  use Parameters
  implicit none
  open(2,file=fname)

  open(24,file='topologia.dat')
  write(24,'(a60)') "#nconf, volm, varv, surm, vars, curm, varc, eulm, vare"

  if (neg) then
     open(25,file='toponeg.dat')
     write(25,'(a60)') "#nconf, volm, varv, surm, vars, curm, varc, eulm, vare"
  endif

  open(40,file='evol_topo.dat')
  write(40,*) "#iconf, vol, sur, cur, eul"
end subroutine Init_Files

subroutine Init_Data
  use Statistics
  volm=0
  surm=0
  curm=0
  eulm=0

  volmn=0
  surmn=0
  curmn=0
  eulmn=0

  volm2=0
  surm2=0
  curm2=0
  eulm2=0

  volm2n=0
  surm2n=0
  curm2n=0
  eulm2n=0
end subroutine Init_Data
 
subroutine Read_Data
  use Parameters
  use Particles
  implicit none
  integer i, j
  real aux
  real xi, xf, yi, yf, zi, zf
  character(20)::  as, a2

  nconf=1
  rclus=rmin
  fname=trim(fname)
  ipbc=.true.
  read(2,*,iostat=fstat)
  if (fstat.lt.0) return
  read(2,*)
  read(2,*)
  read(2,*) npart
  read(2,*)
  read(2,*) xi, xf
  read(2,*) yi, yf
  read(2,*) zi, zf
  read(2,*)
  if (.not.(xf-xi.eq.yf-yi.and.zf-zi.eq.yf-yi)) then
     write(06,'(a)') "Chequear dimensiones de la caja, no es cuadrada"
     stop
  else
     boxsz=xf-xi
  endif
  if(allocated(coords)) deallocate(coords)
  allocate(coords(3,27*npart))
  coords=0
  read_part: do i=1,npart
     read(2,*) aux, aux, (coords(j,i),j=1,3)
  enddo read_part
end subroutine Read_Data

function bool2int(bool)
  logical :: bool
  integer :: bool2int
  if(bool) then
     bool2int=1
  else
     bool2int=0
  endif
  return
end function bool2int

subroutine Init_Lattice
  use Parameters
  use Info_Lattice
  use Particles
  implicit none
  integer ipart
  real (kind = 8) :: xmin, xmax
  real (kind = 8) :: ymin, ymax
  real (kind = 8) :: zmin, zmax
  real (kind = 8) :: dx, dy, dz

  xi=0
  xf=boxsz

  pbc: if (.not.ipbc) then 
     
     xmin=minval(coords(1,:))
     xmax=maxval(coords(1,:))
     ymin=minval(coords(2,:))
     ymax=maxval(coords(2,:))
     zmin=minval(coords(3,:))
     zmax=maxval(coords(3,:))
     
     dx=xmax-xmin
     dy=ymax-ymin
     dz=zmax-zmin
     
     if (dx>=max(dy,dz)) then
        xf=xmax
        xi=xmin
     else if (dy>=max(dx,dz)) then
        xf=ymax
        xi=ymin
     else
        xf=zmax
        xi=zmin
     endif
  else
     ipextra=0
     do ipart=1,npart
        ! ipart cerca de x=0
        if (coords(1,ipart) <= rmin+rcell) then
           ipextra=ipextra+1
           coords(1,npart+ipextra)=coords(1,ipart)+boxsz !la copio
           coords(2,npart+ipextra)=coords(2,ipart) !la copio
           coords(3,npart+ipextra)=coords(3,ipart) !la copio

           !ipart cerca de x=boxsz
        else if (coords(1,ipart) >= boxsz-(rmin+rcell)) then
           ipextra=ipextra+1
           coords(1,npart+ipextra)=coords(1,ipart)-boxsz !la copio
           coords(2,npart+ipextra)=coords(2,ipart)
           coords(3,npart+ipextra)=coords(3,ipart)
        endif
        ! ipart cerca de y=0
        if (coords(2,ipart) <= rmin+rcell) then
           ipextra=ipextra+1
           coords(1,npart+ipextra)=coords(1,ipart)
           coords(2,npart+ipextra)=coords(2,ipart)+boxsz !la copio
           coords(3,npart+ipextra)=coords(3,ipart)
           !ipart cerca de y=boxsz
        else if (coords(2,ipart) >= boxsz-(rmin+rcell)) then
           ipextra=ipextra+1
           coords(1,npart+ipextra)=coords(1,ipart)
           coords(2,npart+ipextra)=coords(2,ipart)-boxsz !la copio
           coords(3,npart+ipextra)=coords(3,ipart)
        endif
        ! ipart cerca de z=0
        if (coords(3,ipart) <= rmin+rcell) then
           ipextra=ipextra+1
           coords(1,npart+ipextra)=coords(1,ipart)
           coords(2,npart+ipextra)=coords(2,ipart)
           coords(3,npart+ipextra)=coords(3,ipart)+boxsz !la copio
           !ipart cerca de z=boxsz
        else if (coords(3,ipart) >= boxsz-(rmin+rcell)) then
           ipextra=ipextra+1
           coords(1,npart+ipextra)=coords(1,ipart)
           coords(2,npart+ipextra)=coords(2,ipart)
           coords(3,npart+ipextra)=coords(3,ipart)-boxsz !la copio
        endif
     enddo
  endif pbc
  
  ! Dimension lineal de la caja digitalizada
  lado=int((xf-xi)/rcell)
  lx=lado
  ly=lado
  lz=lado
  
  if (allocated(iarray)) deallocate(iarray)
  allocate(iarray(0:lx-1,0:ly-1,0:lz-1))
  
  iarray=.false.
  
  rcellm=rcell*0.5
  dd=rmin*rmin
  
  ncells=lx*ly*lz
  if(allocated(lattice)) deallocate(lattice)
  allocate(lattice(0:ncells-1))

  lattice=.false.
  nbits=0
end subroutine Init_Lattice

subroutine Digitalize
  use Parameters
  use Particles
  use Info_Lattice
  implicit none
  integer i, j, k
  integer ilat, ipart
  real (kind=8) :: x, y, z
  real (kind=8) :: x2, y2, z2
  real (kind=8) :: d2
  logical :: flag
  do i=0,lz-1
     z=xi+i*rcell+rcellm  
     do j=0,ly-1
        y=xi+j*rcell+rcellm  
        do k=0,lx-1
           x=xi+k*rcell+rcellm  
           flag=.false.
           do ipart=1,npart+ipextra
              x2=(x-coords(1,ipart))**2
              if (x2.gt.dd) cycle
              y2=(y-coords(2,ipart))**2
              if (y2+x2.gt.dd) cycle
              z2=(z-coords(3,ipart))**2
              d2=x2+y2+z2
              if(d2.lt.dd) then
                 flag=.true.
                 exit
              endif
           enddo
           ilat=i+lx*(j+ly*k)
           if(flag) then
              nbits=nbits+1
              iarray(i,j,k)=.true.
              lattice(ilat)=.true.
           else
              nbitsn=nbitsn+1
           endif
        enddo
     enddo
  enddo
end subroutine Digitalize

subroutine Write_Data
  use Info_Lattice
  use Parameters
  use Statistics
  implicit none
  integer i, j, k
  real (kind=8) :: x, y, z
  
  open(10,file='digit.XYZ')
  write(10,*) nbits
  write(10,*) 1.

  if (neg) then
     open(11,file='digitneg.XYZ')
     write(11,*) nbitsn
     write(11,*) 1.
  endif
  
  
  do i=0,lz-1
     z=xi+i*rcell+rcellm  
     do j=0,ly-1
        y=xi+j*rcell+rcellm  
        do k=0,lx-1
           x=xi+k*rcell+rcellm
           if(iarray(i,j,k)) then
              write(10,*) 'H',  x, y, z
           else
              if (neg) write(11,*) 'He', x, y, z
           endif
        enddo
     enddo
  enddo
  write(24,'(i6, 8(f12.3))') &
       nconf, volm, volm2, surm, surm2, curm, curm2, eulm, eulm2
  close(10)
  close(24)
  
  if (neg) then
     write(25,'(i6, 8(f12.3))') &
          nconf, volmn, volm2n, surmn, surm2n, curmn, curm2n, eulmn, eulm2n
     close(11)
     close(25)
  endif
  
  close(40)
end subroutine Write_Data

subroutine Update_Statistics
  use Parameters
  use Functionals
  use Statistics
  implicit none

  volm=volm+voli*rcell**3
  surm=surm+suri*rcell**2
  curm=curm+curi*rcell
  eulm=eulm+euli
  
  volm2=volm2+voli*voli*rcell**6
  surm2=surm2+suri*suri*rcell**4
  curm2=curm2+curi*curi*rcell**2
  eulm2=eulm2+euli*euli
  
  write(40,*) iconf, vol, sur, cur, eul
  
end subroutine Update_Statistics

subroutine Update_Statistics_Neg
  use Parameters
  use Functionals
  use Particles
  use Statistics
  implicit none
  volmn=volmn+voli*(rcell/boxsz)**3
  surmn=surmn+suri*(rcell/boxsz)**2
  curmn=curmn+curi*(rcell/boxsz)
  eulmn=eulmn+euli
  
  volm2n=volm2n+voli*voli*(rcell/boxsz)**6
  surm2n=surm2n+suri*suri*(rcell/boxsz)**4
  curm2n=curm2n+curi*curi*(rcell/boxsz)**2
  eulm2n=eulm2n+euli*euli
end subroutine Update_Statistics_Neg

subroutine Get_Statistics
  use Statistics
  use Parameters
  implicit none
  volm=volm/nconf
  surm=surm/nconf
  curm=curm/nconf
  eulm=eulm/nconf

  volm2=volm2/nconf
  surm2=surm2/nconf
  curm2=curm2/nconf
  eulm2=eulm2/nconf

  volm2=sqrt((volm2-volm*volm))
  surm2=sqrt((surm2-surm*surm))
  curm2=sqrt((curm2-curm*curm))
  eulm2=sqrt((eulm2-eulm*eulm))

  if (neg) then
     volmn=volmn/nconf
     surmn=surmn/nconf
     curmn=curmn/nconf
     eulmn=eulmn/nconf
     
     volm2n=volm2n/nconf
     surm2n=surm2n/nconf
     curm2n=curm2n/nconf
     eulm2n=eulm2n/nconf

     volm2n=sqrt((volm2n-volmn*volmn))
     surm2n=sqrt((surm2n-surmn*surmn))
     curm2n=sqrt((curm2n-curmn*curmn))
     eulm2n=sqrt((eulm2n-eulmn*eulmn))
  endif

end subroutine Get_Statistics

subroutine Minkowski_Functionals
  use Info_Lattice
  use Statistics
  use Parameters
  use Functionals
  implicit none
  logical, allocatable :: tmp(:)
  integer :: idx,i,j,k
  integer :: ilat,itmp
  integer :: ncellstmp
  vol=0
  sur=0
  cur=0
  eul=0
  ! Workspace
  ncellstmp=(lx+2)*(ly+2)*(lz+2)
  if(allocated(tmp)) deallocate(tmp)
  allocate(tmp(0:ncellstmp-1))
  tmp=.false.
  
  pbc: if (ipbc) then
     do j=0,ly-1
        do i=0,lx-1
           ilat=( i )+( lx )*(( j )+( ly )*(lz-1))
           itmp=(i+1)+(lx+2)*((j+1)+(ly+2)*( 0  ))
           tmp(itmp)=lattice(ilat)
        enddo
     enddo

     do k=0,lz-1
        do i=0,lx-1
           ilat=( i )+( lx )*((ly-1)+( ly )*( k ))
           itmp=(i+1)+(lx+2)*(( 0  )+(ly+2)*(k+1))
           tmp(itmp)=lattice(ilat)
        enddo
     enddo

     do k=0,lz-1
        do j=0,ly-1
           ilat=(lx-1)+( lx )*(( j )+( ly )*( k ))
           itmp=( 0  )+(lx+2)*((j+1)+(ly+2)*(k+1))
           tmp(itmp)=lattice(ilat)
        enddo
     enddo
  endif pbc

  !  Ahora calcula
  do k=0,lz-1
     do j=0,ly-1 
        do i=0,lx-1
           idx=i+lx*(j+ly*k)
           if(lattice(idx)) then
              call Minko_3D(lx+2,ly+2,i+1,j+1,k+1,&
                   tmp,volume,surface,curvature,euler3D,&
                   ncellstmp)
              tmp(i+1+(lx+2)*(j+1+(ly+2)*(k+1)))=.true.
              vol=vol+volume
              sur=sur+surface
              cur=cur+curvature
              eul=eul+euler3D
           endif
        enddo
     enddo
  enddo
  voli=vol
  suri=sur
  curi=cur
  euli=eul
  deallocate(tmp)
  return
end subroutine Minkowski_Functionals
  
subroutine Minko_3D(lx,ly,jx,jy,jz,lattice,&
     volume,surface,curvature,euler3D,&
     ncellstmp)
  implicit none
  integer bool2int
  Integer, parameter  ::  volume_body=1,  &   !(a*a*a, where a is 
                                ! lattice displacement)
       surface_body=-6,&   !(!6*a*a, open body)
       surface_face=2, &   !(2*a*a,open face)
       curv_body=3,    &   !(3*a, open body)
       curv_face=-2,   &   !(!2*a, open face)
       curv_edge=1,    &   !(a, open line)
       euler3D_body=-1,&   !(open body)
       euler3D_face=1, &   !(open face)
       euler3D_edge=-1,&   !(open line)
       euler3D_vertex=1    !(vertices)
  !
  logical     ::  lattice(0:ncellstmp-1)
  Integer     ::  nfaces,nedges,nvert,ncellstmp
  Integer     ::  jx,jy,jz,lx,ly
  Integer     ::  volume,surface,curvature,euler3D
  Integer     ::  i0,jxi,jyi,kc1,kc2,kc3,j0,jyj,jzj,k4,k7,kc7,&
       kc1kc4kc5,k0,jzk,jzi,k9,k10

  nfaces=0
  nedges=0
  nvert=0
  !
  do i0=-1,1,2
     jxi=jx+i0
     jyi=jy+i0
     jzi=jz+i0
     kc1=1-bool2int(lattice(jxi+Lx*(jy+Ly*jz)))
     kc2=1-bool2int(lattice(jx+Lx*(jyi+Ly*jz)))
     kc3=1-bool2int(lattice(jx+Lx*(jy+Ly*jzi)))
     nfaces=nfaces+kc1+kc2+kc3
     do j0=-1,1,2
        jyj=jy+j0  
        jzj=jz+j0
        k4=Lx*(jyj+Ly*jz)
        k7=Lx*(jy+Ly*jzj)
        kc7=1-bool2int(lattice(jx+k7))
        kc1kc4kc5=kc1*(1-bool2int(lattice(jxi+k4)))*(1-bool2int(lattice(jx+k4)))
        nedges=nedges+kc1kc4kc5+&
             kc2*(1-bool2int(lattice(jx+Lx*(jyi+Ly*jzj))))*kc7+&
             kc1*(1-bool2int(lattice(jxi+k7)))*kc7
        if(kc1kc4kc5.ne.0) then
           do k0=-1,1,2
              jzk=jz+k0
              k9=Lx*(jy+Ly*jzk)
              k10=Lx*(jyj+Ly*jzk)
              nvert=nvert+(1-bool2int(lattice(jxi+k9)))*&
                   (1-bool2int(lattice(jxi+k10)))*&
                   (1-bool2int(lattice(jx+k10)))*(1-bool2int(lattice(jx+k9)))
           enddo ! k0
        endif ! kc1kc4kc5
     enddo ! j0
  enddo !i0
  !
  volume=volume_body
  surface=surface_body+surface_face*nfaces
  curvature=curv_body+curv_face*nfaces+curv_edge*nedges
  euler3D=euler3D_body+euler3D_face*nfaces+&
       euler3D_edge*nedges+euler3D_vertex*nvert
end subroutine Minko_3D

subroutine End_Data
  use Particles
  use Info_Lattice
  deallocate(coords)
  deallocate(iarray)
  deallocate(lattice)
end subroutine End_Data
