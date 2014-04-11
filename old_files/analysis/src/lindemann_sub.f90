subroutine Init_Data
  use Particles
  use Parameters
  use Statistics
  implicit none
  real xi, xf, yi, yf, zi, zf
  open(10, file=trim(fname))
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*) npart
  read(10,*)
  read(10,*) xi, xf
  read(10,*) yi, yf
  read(10,*) zi, zf
  read(10,*)
  if (.not.(xf-xi.eq.yf-yi.and.zf-zi.eq.yf-yi)) then
     write(06,'(a)') "Chequear dimensiones de la caja, no es cuadrada"
     stop
  else
     boxsz=xf-xi
  endif
  rewind(10)
  allocate (x(3,npart))
  allocate (r(3,npart))
  allocate (r2(3,npart))
  x=0
  r=0
  r2=0
end subroutine Init_Data

subroutine End_Data
  use Particles
  use Statistics
  deallocate (x)
  deallocate (r)
  deallocate (r2)
end subroutine End_Data

subroutine Read_Data
  use Parameters
  use Particles
  implicit none
  integer i, j
  integer iso
  real aux
  real xi, xf, yi, yf, zi, zf
  character(20)::  as, a2
  
  fname=trim(fname)
  read(10,*,iostat=fstat)
  if (fstat.lt.0) return
  read(10,*)
  read(10,*)
  read(10,*) npart
  read(10,*)
  read(10,*) xi, xf
  read(10,*) yi, yf
  read(10,*) zi, zf
  read(10,*)
  if (.not.(xf-xi.eq.yf-yi.and.zf-zi.eq.yf-yi)) then
     write(06,'(a)') "Chequear dimensiones de la caja, no es cuadrada"
     stop
  else
     boxsz=xf-xi
  endif
  read_part: do i=1,npart
     read(10,*) iso, aux, (x(j,i),j=1,3)
  enddo read_part
end subroutine Read_Data

subroutine Get_Statistics
  use Parameters
  use Particles
  use Statistics
  implicit none
  integer i, j
  
  do i=1,npart
     do j=1,3
        r (j,i) = r (j,i) + x(j,i)   
        r2(j,i) = r2(j,i) + x(j,i)**2
     enddo
  enddo
end subroutine Get_Statistics
  
  

subroutine Get_Lind
  use Particles
  use Parameters
  use Statistics
  implicit none
  integer i, j
  
  lind=0
  do i=1,npart
     do j=1,3
        lind=lind+r2(j,i)/nconf-(r(j,i)/nconf)**2
     enddo
  enddo
  lind=sqrt(lind/npart)/boxsz
end subroutine Get_Lind
