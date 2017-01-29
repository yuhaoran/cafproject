subroutine finalize
use variables
use cubefft
use penfft
implicit none
save
!include 'fftw3.f'

!call sfftw_destroy_plan(plan_fft_fine)
!call sfftw_destroy_plan(plan_ifft_fine)
call destroy_cubefft_plan

!call sfftw_destroy_plan(planx)
!call sfftw_destroy_plan(iplanx)
!call sfftw_destroy_plan(plany)
!call sfftw_destroy_plan(iplany)
!call sfftw_destroy_plan(planz)
!call sfftw_destroy_plan(iplanz)

call destroy_penfft_plan

endsubroutine
