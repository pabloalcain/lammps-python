module Parameters
  character(20) :: fname
  integer :: nbins
  real :: deltar, rmax
  integer :: iconf=0, nconf=0
  logical :: ipbc
  integer :: fstat
end module Parameters

module Particles
  integer       :: npart, nprot, nneut
  real :: boxsz
  real,  dimension(:,:), allocatable :: x
  logical, dimension(:), allocatable :: isos
  double precision,  dimension(:,:), allocatable :: gr
end module Particles

module Statistics
  integer :: pasada
end module Statistics
