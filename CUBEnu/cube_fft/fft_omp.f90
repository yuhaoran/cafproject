program test
  use omp_lib
  use iso_fortran_env, only : int64
  implicit none
  include 'fftw3.f'
  integer,parameter :: ng=512 ! number of grid per dimension
  real,parameter :: pi=3.14159
  integer(8) plan_fft,plan_ifft ! fft plans
  real den(ng+2,ng,ng) ! fft array, additional "+2" for storing Nyquist frequency
  real temp_r,temp_theta

  integer i,j,k,l,seedsize
  integer(int64) time64
  integer(4) i4,t1,t2,tt1,tt2,ttt1,ttt2,t_rate,stat

  integer,allocatable :: iseed(:)




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print*, 'ilp64 test: default integer'
  i=2147483647
  i=i+1
  print*,i

  ! fft plan
  call sfftw_init_threads(stat)
  print*, stat
  stop
  call system_clock(t1,t_rate)
  print*, 'call forward plan'
  call sfftw_plan_dft_r2c_3d(plan_fft,ng,ng,ng,den,den,FFTW_MEASURE)
  print*, 'call backward plan'
  call sfftw_plan_dft_c2r_3d(plan_ifft,ng,ng,ng,den,den,FFTW_MEASURE)
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  call system_clock(t1,t_rate)
  print*, 'initialize den'
  call random_seed(size=seedsize)
  !print*,'min seedsize =', seedsize
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
  call system_clock(time64)
  do i=1,seedsize
    iseed(i)=lcg(time64)
    !print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
  enddo
  call random_seed(put=iseed)
  call random_number(den(:ng,:,:))
  deallocate(iseed)
  print*,'Box-Muller transform'
  do k=1,ng
  do j=1,ng
  do i=1,ng,2
    temp_theta=2*pi*den(i,j,k)
    temp_r=sqrt(-2*log(1-den(i+1,j,k)))
    den(i,j,k)=temp_r*cos(temp_theta)
    den(i+1,j,k)=temp_r*sin(temp_theta)
  enddo
  enddo
  enddo
  print*, den(1:10,1,1)
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  call system_clock(t1,t_rate)
  print*, 'call forward fft'
  call sfftw_execute(plan_fft) ! forward transform
  print*, 'call backward plan'
  call sfftw_execute(plan_ifft) ! backward transform
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  print*,'check results'
  den=den/ng/ng/ng ! normalize
  print*, den(1:10,1,1)
  print*, 'destroy'
  call sfftw_destroy_plan(plan_fft)
  call sfftw_destroy_plan(plan_ifft)

  !!!!!!!!!!!!!!!!!!!!!!!!

contains
  function lcg(s) !// Linear congruential generator
    implicit none
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  endfunction lcg
end
