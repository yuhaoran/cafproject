module cubefft
  contains
  subroutine create_cubefft_plan
    use variables
    implicit none
    save
    include 'fftw3.f'
    call sfftw_plan_dft_r2c_3d(plan_fft_fine,nfe,nfe,nfe,rho_f,rho_f,FFTW_MEASURE)
    call sfftw_plan_dft_c2r_3d(plan_ifft_fine,nfe,nfe,nfe,rho_f,rho_f,FFTW_MEASURE)
  endsubroutine create_cubefft_plan

  subroutine destroy_cubefft_plan
    use variables
    implicit none
    save
    include 'fftw3.f'
    call sfftw_destroy_plan(plan_fft_fine)
    call sfftw_destroy_plan(plan_ifft_fine)
  endsubroutine destroy_cubefft_plan

endmodule
