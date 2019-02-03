program test
  !use iso_fortran_env, only : int64
  use,intrinsic :: ISO_C_BINDING
  use omp_lib
  implicit none
  include 'fftw3.f03'
  integer,parameter :: ng=304 ! number of grid per dimension
  real,parameter :: pi=3.14159
  integer(8) plan_fft,plan_ifft ! fft plans
  real den(ng+2,ng,ng) ! fft array, additional "+2" for storing Nyquist frequency
  real temp_r,temp_theta
  real total_start_cpu, total_finish_cpu, totaltime_cpu
  double precision :: start, finish, total_start, total_finish
  double precision :: createtime, initialtime, setuptime, executime, totaltime
  integer i,j,k,l,seedsize,stat
  integer(c_int64_t) time64
  integer(4) i4,t1,t2,t_rate,ncore

  integer,allocatable :: iseed(:)




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cpu_time(total_start_cpu)
  total_start=omp_get_wtime() 
  start = omp_get_wtime()
  print*, 'ilp64 test: default integer'
  i=2147483647
  i=i+1
  print*,i

  call sfftw_init_threads(stat)
  write(*,*) "Status:", stat
  write(*,*) "ng=",ng
  ncore=omp_get_max_threads()
  print*, 'max threads =',ncore
  !ncore=4
  call sfftw_plan_with_nthreads(ncore)
  write(*,*) "Using ", ncore, "threads."

  ! fft plan
  call system_clock(t1,t_rate) 
  print*, 'call forward plan'
  call sfftw_plan_dft_r2c_3d(plan_fft,ng,ng,ng,den,den,FFTW_MEASURE)
  print*, 'call backward plan'
  call sfftw_plan_dft_c2r_3d(plan_ifft,ng,ng,ng,den,den,FFTW_MEASURE)
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  ! initialize matrix
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
 
  ! transforms
do i=1,1
  call system_clock(t1,t_rate)
  print*, 'FFT'
  call sfftw_execute(plan_fft) ! forward transform
  print*, 'iFFT'
  call sfftw_execute(plan_ifft) ! backward transform
  den=den/ng/ng/ng ! normalize
  print*, den(1:10,1,1)
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
enddo  

  ! destroy plan
  call system_clock(t1,t_rate)
  print*, 'destroy'
  call sfftw_destroy_plan(plan_fft)
  call sfftw_destroy_plan(plan_ifft)
  call fftw_cleanup_threads()
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';


contains
  function lcg(s) !// Linear congruential generator
    implicit none
    integer :: lcg
    integer(c_int64_t) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_c_int64_t)
    end if
    s = mod(s * 279470273_c_int64_t, 4294967291_c_int64_t)
    lcg = int(mod(s, int(huge(0), c_int64_t)), kind(0))
  endfunction lcg
end
