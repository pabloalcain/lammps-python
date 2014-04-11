program main
  implicit none
  real (kind=8) a,b
  character(20) c
  a=1.5
  b=0.5
  c="evol.lammpstrj"
  print*,c,a,b
  call minkowski(c,a,b)
end program main
