module Parameters
  character(20) :: fname
  real :: lind
  integer :: fstat
  integer :: nconf, iconf
end module Parameters


module Particles
  integer       :: npart
  real :: boxsz
  real,  dimension(:,:), allocatable :: x
end module Particles

module Statistics
  integer :: pasada
  real (kind=8), dimension(:,:), allocatable :: r, r2
end module Statistics
