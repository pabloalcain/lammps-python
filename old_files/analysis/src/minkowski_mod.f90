module Functionals
  real (kind=8)  :: cur=0,vol=0,sur=0,eul=0
  real (kind=8)  :: voli=0,suri=0,curi=0,euli=0
  integer        :: volume,surface,curvature,euler3D
end module Functionals

module Statistics
  real (kind=8)  :: volm=0,surm=0,curm=0,eulm=0
  real (kind=8)  :: volmn=0,surmn=0,curmn=0,eulmn=0
  real (kind=8)  :: volm2=0,surm2=0,curm2=0,eulm2=0
  real (kind=8)  :: volm2n=0,surm2n=0,curm2n=0,eulm2n=0
end module Statistics

module Parameters
  character(20) :: fname
  real (kind=8) :: rmin, rcell, rclus
  integer :: iconf=0, nconf=0
  logical :: ipbc
  integer :: fstat
  logical :: neg
end module Parameters

module Particles
  integer       :: npart, ipextra
  real (kind=8) :: dens, boxsz
  real (kind=8), dimension(:,:), allocatable :: coords
end module Particles

module Info_Lattice
  real (kind=8)        :: xf, xi
  integer              :: lx, ly, lz, lado
  logical, allocatable :: iarray(:,:,:)
  logical, allocatable :: lattice(:)
  integer              :: ncells
  real (kind=8)        :: rcellm, dd
  integer              :: nbits, nbitsn
end module Info_Lattice
