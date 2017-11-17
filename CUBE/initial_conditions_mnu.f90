program mnu_ic
  implicit none

!!- These should be/are in the parameters module -!!
  !Numbers
  logical, parameter :: head = .true.
  real, parameter :: pi = 4.0*atan(1.d0)

  !Strings
  character(len=*), parameter :: output = '/scratch/dbi208/caf/image1'

  !Cosmological parameters
  real, parameter :: a_i_nu = 1./(1.0+10.0) !Starting scalefactor
  real, parameter :: mnu = 0.05 !Neutrino mass eV
  real, parameter :: box = 100.0 !Box Size Mpc/h
  real, parameter :: omega_m = 0.32 !Energy fraction of matter

  !Numerical parameters
  integer, parameter :: nn = 1 !Number of nodes
  integer, parameter :: nc = 256 !Number of cells
  integer, parameter :: nc2p = 2 !Ratio of cells to particles
  !e.g. nc2p=2 -> (nc/2)**3 total particles

  !ID Format
  integer, parameter :: izipi = 4
!!- -!!
  
  !Fermi-Dirac CDF 
  integer, parameter :: ncdf = 10000
  real, dimension(2,ncdf) :: cdf

  !Units
  real, parameter :: vp2s = 1.0/(300.*sqrt(omega_m)*box/a_i_nu/2./nc)
  real, parameter :: fd = 50.2476/mnu/a_i_nu !kBcTnu/mnu
  
  !Seed
  integer(4) seedsize
  integer(4), allocatable :: iseed(:)
  real, allocatable :: rseed_all(:,:)

  !Useful variables and small arrays
  integer :: i,j,k,n

  !Particle information
  integer(izipi) :: iq
  real(8), dimension(3) :: rng
  real, dimension(3) :: xq,vq

  !Check
  logical :: ok

  !Setup
  if (head) then
     write(*,*) ''
     write(*,*) 'Homogeneous Initial Conditions for Massive Neutrinos'
     write(*,*) ''
     write(*,*) 'vp2s/(km/s)=',vp2s
     write(*,*) 'fd/(km/s)=',fd
  end if

  !Read seed
  if (head) write(*,*) 'Reading seeds'
  call random_seed(size=seedsize)
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
  allocate(rseed_all(seedsize,nn**3))
  open(11,file=output//'/seed_1.bin',status='old',access='stream')
  read(11) iseed
  close(11)
  call random_seed(put=iseed)
     
  !Read cdf table
  if (head) write(*,*) 'Reading cdf'
  open(11,file='./CDFTable.txt')
  read(11,*) cdf
  close(11)

  !Create particles
  if (head) write(*,*) 'Computing particle positions and velocities'
  open(unit=81,file=output//'/xv_nu.bin',status='replace',access='stream')
  open(unit=82,file=output//'/id_nu.bin',status='replace',access='stream')
  do k=1,nc,nc2p
     do j=1,nc,nc2p
        do i=1,nc,nc2p

           !Compute Lagrangian positions and Thermal velocities
           
           !Positions stored in xq
           xq=(/i,j,k/)-0.5
           
           !Velocities stored in vq
           call random_number(rng)
           !!Convert from uniform to FD
           ok=.false.
           do n=1,ncdf-1
              if (cdf(2,n+1).gt.rng(1)) then
                 !Linear interpolate
                 rng(1)=(cdf(1,n)*(cdf(2,n+1)-rng(1))+cdf(1,n+1)*(rng(1)-cdf(2,n)))/(cdf(2,n+1)-cdf(2,n))
                 ok=.true.
                 exit
              end if
           end do
           if (.not.ok) ERROR STOP 13

           !!Store fraction of max velocity
           !!!Fermi-Dirac CDF Approximated by Gaussian
           iq=nint(approxCDF(rng(1))*int(2,8)**(8*izipi)-int(2,8)**(8*izipi-1),kind=izipi)

           !!Amplitude and Angle
           rng(1)=rng(1)*fd*vp2s !!Holds velocity amplitude
           rng(2)=2.*rng(2)-1. !cosTheta in (-1,1)
           rng(3)=rng(3)*2.*pi !Phi in 0 to 2*pi
           !!Direction
           vq(1)=rng(1)*sqrt(1.-rng(2)**2.)*cos(rng(3))
           vq(2)=rng(1)*sqrt(1.-rng(2)**2.)*sin(rng(3))
           vq(3)=rng(1)*rng(2)

           !Write out
           write(81) xq
           write(81) vq
           write(82) iq

           !Diagnostics
           if (i+(j-1)*nc+(k-1)*nc**2 .lt. 10 .and. head) then
              write(*,*) ''
              write(*,*) i,j,k
              write(*,*) 'xq',xq
              write(*,*) 'vq',vq,sum(vq**2)**0.5
              write(*,*) 'iq',iq
              write(*,*) 'viq',vp2s*fd*invertCDF(1.d0*(iq+int(2,8)**(8*izipi-1))/int(2,8)**(8*izipi))
           end if

        end do
     end do
  end do
  close(81)
  close(82)

  if (head) write(*,*) 'Finished neutrino ic'

contains

  function approxCDF(v) result(c)
    implicit none
    real(8), intent(in) :: v
    real(8) :: c
    real(8), parameter :: s=3.5
    c=1.-exp(-(v/s)**2.)
  end function approxCDF

  function invertCDF(c) result(v)
    implicit none
    real(8), intent(in) :: c
    real(8) :: v
    real(8), parameter :: s=3.5
    v=s*sqrt(log(1./(1.-c)))
  end function invertCDF

end program mnu_ic
