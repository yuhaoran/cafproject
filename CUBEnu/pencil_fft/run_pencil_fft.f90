program run_pencil_fft
  use iso_fortran_env, only : int64
  use parameters
  use pencil_fft
  implicit none

  real temp_r,temp_theta
  integer i,j,k,seedsize
  integer(int64) time64
  integer,allocatable :: iseed(:)

  call geometry

  if (this_image()==1) print*, 'call fft plan'
  call create_penfft_plan

  if (this_image()==1) print*, 'initialize random number'
  call random_seed(size=seedsize)
  !print*,'min seedsize =', seedsize
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
  call system_clock(time64)
  do i=1,seedsize
    time64=time64*this_image()
    iseed(i)=lcg(time64)
    !print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
  enddo
  call random_seed(put=iseed)
  call random_number(r3)

  deallocate(iseed)

  if (this_image()==1) print*,'Box-Muller transform'
  do k=1,ng
  do j=1,ng
  do i=1,ng,2
    temp_theta=2*pi*r3(i,j,k)
    temp_r=sqrt(-2*log(1-r3(i+1,j,k)))
    r3(i,j,k)=temp_r*cos(temp_theta)
    r3(i+1,j,k)=temp_r*sin(temp_theta)
  enddo
  enddo
  enddo
  !r3=this_image()
  r0=r3(:,:,ng)
  print*, r3(1,1,1),r3(ng,ng,ng)

  sync all

  if (head) print*, 'call ftran'
  call pencil_fft_forward
  if (head) print*, 'call btran'
  call pencil_fft_backward
  sync all
  print*, 'precision on image',this_image(),maxval(abs(r3(:,:,ng)-r0))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
