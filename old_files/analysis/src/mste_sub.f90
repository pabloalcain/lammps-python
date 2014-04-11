subroutine Init_Files
  use Parameters
  use Particles
  implicit none

  open(11,file='mste_clus.dat')
  write(11,*) fname, npart, nprot
  
  open(13,file='mste_dist.dat')
  write(13,*) '#mass, total n. of fragments, means'
end subroutine Init_Files

subroutine Init_Data
  use Particles
  use mod_mst
  implicit none
  real xi, xf, yi, yf, zi, zf
  integer iso,i
  
  nprot=0
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
  endif
  read_part: do i=1,npart
     read(10,*) iso
     if (iso.eq.2) nprot=nprot+1
  enddo read_part
  rewind(10)
  nprot=npart/2
  allocate (x(3,npart))
  allocate (p(3,npart))
  allocate (x_shifted(3,npart))
  allocate (isos(npart))
  allocate (lista_part_ord(npart))
  allocate (lista_part_aux(npart))
  allocate (lista_part_rem(npart))
  allocate (index_part_clu(npart))
  allocate (bordes_clusters(0:npart))
  allocate (mdistt(npart))
  allocate (mass_dist(npart))
  mdistt=0
  mass_dist=0
end subroutine Init_Data

subroutine End_Data
  use Particles
  use mod_mst
  deallocate (x)
  deallocate (p)
  deallocate (x_shifted)
  deallocate (isos)
  deallocate (lista_part_ord)
  deallocate (lista_part_aux)
  deallocate (lista_part_rem)
  deallocate (index_part_clu)
  deallocate (bordes_clusters)
  deallocate (mdistt)
  deallocate (mass_dist)
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
  
  nconf=1
  ipbc=.true.
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
     read(10,*) iso, aux, (x(j,i),j=1,3), (p(j,i),j=1,3)
     isos(i)=(iso.eq.2)
  enddo read_part
  p=p*mass
end subroutine Read_Data

subroutine Get_Statistics
  use mod_mst
  use Particles
  use Statistics
  implicit none
  mdistt=mdistt+mass_dist
end subroutine Get_Statistics

subroutine Write_Files
  use mod_mst
  use Particles
  use Parameters
  use Statistics
  implicit none
  integer i
  integer nmass

  nmass=0
  do i = 1, npart
     if(mass_dist(i)>0) then
        nmass = nmass + 1
     endif
  enddo


200 format(214(i5,1x))
997 format(2(i6,1x),2(f12.6,1x))
998 format(200(i5,1x))
  
  write(11,997) iconf, nmass, 1.0, 1.0
  write(11,200)(lista_part_ord(i),i=1,npart)     !particle indexes
  write(11,200)(bordes_clusters(i),i=1,npart)    !cluster limits
  
  do i = 1,npart
     if(mass_dist(i) > 0 ) then
        write(11,998,advance='no') i
     endif
  enddo
  write(11,*)
  do i = 1,npart
     if(mass_dist(i) > 0 ) then
        write(11,998,advance='no') mass_dist(i)
     endif
  enddo
end subroutine Write_Files

subroutine Write_Data
  use Parameters
  use Particles
  use mod_mst
  implicit none
  integer i

  do i=1,npart
     write(13,*) i, mdistt(i),mdistt(i)*(1.0/(nconf*npart))
  enddo
  close(13)
  close(11)
end subroutine Write_Data
