subroutine mste(i_fname)
  use mod_mst
  use Parameters
  implicit none
  integer     :: nclus_max
  integer     :: iaux1, iaux2, isort1, isort2, ipart, kpart, iclus
  integer     :: i
  character(20), intent(in) :: i_fname
  fname=i_fname
  call Init_Data
  call Init_Files

  iconf=0
  do while (.true.)
     call Read_Data
     if (fstat.lt.0) exit
     iconf=iconf+1
!     write(06,*), "Configuracion", iconf
     call mst
     call Write_Files
     call Get_Statistics
  enddo
  nconf=iconf
  call Write_Data
  call End_Data
  close(10)
  return
end subroutine mste
