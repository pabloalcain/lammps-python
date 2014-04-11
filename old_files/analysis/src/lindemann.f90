subroutine lindemann(i_fname, o_lind)
  use Parameters
  implicit none
  real, intent (out) :: o_lind
  character(20), intent(in) :: i_fname
  
  fname=i_fname
  call Init_Data
  iconf=0
  do while (.true.)
     call Read_Data
     if (fstat.lt.0) exit
     iconf=iconf+1
!     write(06,*), "Configuracion", iconf
     call Get_Statistics
  enddo
  nconf=iconf
  call Get_Lind
  close(10)
  o_lind=lind
  call End_Data
  return
end subroutine lindemann
