implicit none
include 'fftw3.f'
integer,parameter :: ng=512 ! number of grid per dimension
integer(8) plan_fft,plan_ifft ! fft plans
real den(ng+2,ng,ng) ! fft array, additional "+2" for storing Nyquist frequency
call sfftw_plan_dft_r2c_3d(plan_fft,ng,ng,ng,den,den,FFTW_MEASURE)
call sfftw_plan_dft_c2r_3d(plan_ifft,ng,ng,ng,den,den,FFTW_MEASURE)
den(:ng,:,:)=1.23
call sfftw_execute(plan_fft) ! forward transform
call sfftw_execute(plan_ifft) ! backward transform
den=den/ng/ng/ng ! normalize
call sfftw_destroy_plan(plan_fft)
call sfftw_destroy_plan(plan_ifft)
end
