subroutine gofr(i_fname, i_nbins)
  use Parameters
  implicit none
  integer, intent(in) :: i_nbins
  character(20), intent(in) :: i_fname
  
  nbins=i_nbins
  fname=i_fname
  call Init_Data
  call Init_Files
  iconf=0
  do while (.true.)
     call Read_Data
     if (fstat.lt.0) exit
     iconf=iconf+1
     !write(06,*), "Configuracion", iconf
     call Calc_gr
  enddo
  nconf=iconf
  call Write_Data
  call End_Data
end subroutine gofr
