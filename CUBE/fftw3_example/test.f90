implicit none
include 'fftw3.f'
integer,parameter :: ng=2048 ! number of grid per dimension
integer(8) plan_fft,plan_ifft ! fft plans
real den(ng+2,ng,ng) ! fft array, additional "+2" for storing Nyquist frequency

integer i
integer(4) i4

print*, 'ilp64 test: default integer'
i=2147483647
i=i+1
print*,i


print*, 'call forward plan'
call sfftw_plan_dft_r2c_3d(plan_fft,ng,ng,ng,den,den,FFTW_MEASURE)
print*, 'call backward plan'
call sfftw_plan_dft_c2r_3d(plan_ifft,ng,ng,ng,den,den,FFTW_MEASURE)
print*, 'initialize den'
den(:ng,:,:)=1.23
print*, 'call forward fft'
call sfftw_execute(plan_fft) ! forward transform
print*, 'call backward plan'
call sfftw_execute(plan_ifft) ! backward transform
den=den/ng/ng/ng ! normalize
print*, 'destroy'
call sfftw_destroy_plan(plan_fft)
call sfftw_destroy_plan(plan_ifft)
end
