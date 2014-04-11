subroutine Init_Files
  use Parameters
  use Particles
  implicit none

  open(11,file='gr.dat')
  write(11,*) '#r, grall, grnn, grpp, grnp'
end subroutine Init_Files

subroutine Init_Data
  use Particles
  use Parameters
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
  allocate (isos(npart))
  allocate (gr(4,nbins))
  gr=0
  rmax=boxsz/2
  deltar=rmax/nbins
end subroutine Init_Data

subroutine End_Data
  use Particles
  deallocate(x)
  deallocate(isos)
  deallocate(gr)
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
  ipbc=.true.
  read(10,*,iostat=fstat)
  if (fstat.lt.0) return
  nprot=0
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
     isos(i)=(iso.eq.2)
     if (isos(i)) nprot=nprot+1
  enddo read_part
  nneut=npart-nprot
end subroutine Read_Data

subroutine Calc_gr
  use Parameters
  use Particles
  implicit none
  integer i, j
  integer px, py, pz
  integer, dimension(3) :: image
  real, dimension(3) :: r
  real dist
  integer, parameter :: all=1, nn=2, pp=3, np=4
  integer bin
  real rdig, acum

  part: do i=1,npart
     do j=i+1,npart
        r=x(:,i)-x(:,j)
        where(r>boxsz/2) r=r-boxsz
        where(r<-boxsz/2) r=r+boxsz
        dist=dot_product(r,r)
        dist=sqrt(dist)
        if (dist.ge.rmax) cycle
        bin=ceiling(dist*nbins/rmax)
        rdig=((bin-0.5)*deltar)
        acum=1.0/rdig**2
        same: if (isos(i).eqv.isos(j)) then
           prot: if (isos(i)) then
              gr(pp,bin)=gr(pp,bin)+acum
           else
              gr(nn,bin)=gr(nn,bin)+acum
           end if prot
        else
           gr(np,bin)=gr(np,bin)+acum
        end if same
        gr(all,bin)=gr(np,bin)+acum
     enddo
  enddo part
end subroutine Calc_gr

subroutine Write_Data
  use Particles
  use Parameters
  implicit none
  integer i, j
  integer, parameter :: all=1, nn=2, pp=3, np=4
  real, parameter :: pi=2*acos(0.)

  gr=gr*boxsz**3/(4*pi*deltar)

  gr(all,:)=4*gr(all,:)/(nconf*npart**2)
  gr(pp ,:)=2*gr(pp ,:)/(nconf*nprot**2)
  gr(nn ,:)=2*gr(nn ,:)/(nconf*nneut**2)
  gr(np ,:)=1*gr(np ,:)/(nconf*nprot*nneut)
  
  do i=1,nbins
     write(11,*) (i-0.5)*deltar, (gr(j,i),j=1,4)
  enddo
  close(11)
  close(10)
end subroutine Write_Data
