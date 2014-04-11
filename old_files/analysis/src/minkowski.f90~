subroutine minkowski(i_fname, i_rmin, i_rcell,i_neg)
  use Functionals
  use Parameters
  use Info_Lattice
  implicit none
  character(20) :: i_fname
  real (kind=8) :: i_rmin, i_rcell
  logical :: i_neg

  neg=i_neg
  fname=i_fname
  rmin=i_rmin
  rcell=i_rcell

  call Init_Files
  call Init_Data
  
  iconf=0
  do while (.true.)
     call Read_Data
     if (fstat.lt.0) exit
     iconf=iconf+1
!     write(06,*), "Configuracion", iconf
     call Init_Lattice
     call Digitalize
     call Minkowski_Functionals
     call Update_Statistics
     if (neg) then
        lattice = .not.lattice
        call Minkowski_Functionals
        call Update_Statistics_Neg
     endif
  enddo
  close(2)
  nconf=iconf
  call Get_Statistics
  call Write_Data
  call End_Data
  return
end subroutine minkowski
