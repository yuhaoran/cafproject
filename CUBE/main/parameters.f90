module parameters
  implicit none
  save

  ! output directory for both IC and snapshots
  character(*),parameter :: opath='../output/universe1/'

  ! simulation parameters
  integer,parameter :: izipx=2 ! 1 or 2, integer*? for particle location
  integer,parameter :: izipv=2 ! 1 or 2, integer*? for particle velocity
  integer(izipx),parameter :: ishift=-(2**(izipx*8-1))
  real(8),parameter :: rshift=0.5-ishift

  ! (hereafter 'number of fine cells' = 'nf')
  ! (hereafter 'number of coarse cells' = 'nc')
  ! (hereafter 'per dimension' = '/dim')
  integer,parameter :: nn=2 ! number of imgages (nodes) /dim
  integer,parameter :: ncell=4 ! number of nf in each nc, /dim
  integer,parameter :: nnt=2 ! number of tiles /image/dim
  integer,parameter :: nc=24 ! nc/image/dim, in physical volume, >=24
  integer,parameter :: nt=nc/nnt ! nc/tile/dim, in physical volume, >=12

  integer,parameter :: nf=nc*ncell ! >=96
  integer,parameter :: nft=nt*ncell ! >=48

  ! ngrid /image/dim for pencil-fft
# ifdef FFTFINE
    integer,parameter :: ng=nf ! fine grid fft, for IC, dsp, convert_zip_to_xv
# else
    integer,parameter :: ng=nc ! coarse grid fft, for N-body main code
# endif
  integer,parameter :: npen=ng/nn ! ng /dim in shorter side of the pencil, for pencil decomposition

  integer,parameter :: ncore=1 ! number of cores per image
  integer,parameter :: n_nest=4 ! number of nested threads

  integer,parameter :: ncb=6 ! nc in buffer /dim, single side; 6 by default
  integer,parameter :: nce=nc+2*ncb ! extended nc
  integer,parameter :: nte=nt+2*ncb ! extended nt

  integer,parameter :: nfb=ncb*ncell ! 24
  integer,parameter :: nf_cutoff=16 ! beyond this length, fine force is zero

  integer,parameter :: nfe=nft+2*nfb ! 96

  integer,parameter :: np_nc=ncell ! number of particles / coarse cell / dim

  real,parameter :: rsoft=0.1 ! PP softening length
  logical,parameter :: np_2n3=.false. ! if there are 2*N**3 particles

  ! cosmological parameters
  real,parameter :: z_i=20.0   ! initial redshift
  real,parameter :: z_i_nu=z_i ! initial redshift for neutrinos
  real,parameter :: a_i=1/(1+z_i) ! initial scale factor

  real,parameter :: box=100.0  ! simulation scale /dim, in unit of Mpc/h
  real,parameter :: h0=67.74    ! Hubble constant
  real,parameter :: s8=0.8276   ! \sigma_8
  real,parameter :: ratio_nudm_dim=2 ! ratio of number of particles for neutrino/CDM, /dim
  real,parameter :: m_neu=0.05 ! neutrino mass
  real,parameter :: omega_nu=m_neu/93.14/(h0/100.)**2
  real,parameter :: pi=4*atan(1.)


  !!real,parameter :: omega_c=0.27
  !!real,parameter :: omega_b=0.05
  !!real,parameter :: omega_m=omega_c+omega_b+omega_nu
  !!real,parameter :: omega_l=1-omega_m
  !!real,parameter :: omega_l=0.73
  !!real,parameter :: omega_m=1-omega_l

  real,parameter :: omega_c=0.2588964265232328
  real,parameter :: omega_b=0.048597561544344206
  real,parameter :: omega_m=omega_c+omega_b !+omega_nu
  real,parameter :: omega_l=1-omega_m

  real,parameter :: omega_ch=0.
  real,parameter :: bias=1
  real,parameter :: power_index=2
  real,parameter :: wde=-1

  real,parameter :: f_nl=0
  real,parameter :: g_nl=0
  real,parameter :: n_s=0.96
  real,parameter :: scalar_amp=2.142e-9

  integer,parameter :: nts=400 ! maximum number of timesteps
  real,parameter :: ra_max=0.1
  real(8),parameter :: v_resolution=2.1/(2**(izipv*8))
  real(8),parameter :: x_resolution=1.0/2**(izipx*8)

  ! transfer function
  integer, parameter      :: nk_tf=2000

  !! MPI image variables !!
  integer(4) :: rank,icx,icy,icz,inx,iny,inz,ipx,ipy,ipz
  logical head
  ! checkpoint variables
  integer,parameter :: nmax_redshift=100
  integer cur_checkpoint, n_checkpoint[*]
  real z_checkpoint(nmax_redshift)[*]
  logical checkpoint_step[*], final_step[*]

  ! 128 byte (equivalent 32 4-byte variables) header in zip2
  type sim_header
    ! standard cubep3m 18 variables
    integer nplocal
    real a, t, tau
    integer nts
    real dt_f_acc, dt_pp_acc, dt_c_acc
    integer cur_checkpoint,cur_proj,cur_halo
    real mass_p
    real v_i2r(3)
    real shake_offset(3) ! -18
    ! more simulation config info
    real box
    integer(4) rank
    integer(2) nn,nnt,nt,ncell,ncb
    integer(1) izipx,izipv

    ! cosmology
    real h0
    real omega_m
    real omega_l
    real s8
    real m_neu(3)
    real vsim2phys
    real z_i
  endtype

  type(sim_header) sim
  integer m1,m2,m3,m

  contains
    subroutine geometry
      rank=this_image()-1            ! MPI_rank
      icz=rank/(nn**2)+1             ! image_z
      icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
      icx=mod(rank,nn)+1             ! image_x

      m1=icx ! pencil_fft convension
      m2=icy
      m3=icz
      m=num_images()

      ! adjacent images
      inx=modulo(icx-2,nn)+1
      iny=modulo(icy-2,nn)+1
      inz=modulo(icz-2,nn)+1
      ipx=modulo(icx,nn)+1
      ipy=modulo(icy,nn)+1
      ipz=modulo(icz,nn)+1

      head=(this_image()==1)

      sync all
    endsubroutine

    subroutine print_header(s)
      type(sim_header),intent(in) :: s
      if (this_image()==1) then
      print*,'-------------------------------- CUBE info --------------------------------'
      print*,'| nplocal      =',s%nplocal
      print*,'| a,t,tau      =',s%a,s%t,s%tau
      print*,'| nts          =',s%nts
      print*,'| dt f,pp,c    =',s%dt_f_acc,s%dt_pp_acc,s%dt_c_acc
      print*,'| cur_steps    =',int(s%cur_checkpoint,2),int(s%cur_proj,2),int(s%cur_halo,2)
      print*,'| mass_p       =',s%mass_p
      print*,'| v_i2r        =',s%v_i2r, '(1.0)/1'
      print*,'| shake_offset =',s%shake_offset
      print*,'| '
      print*,'| box/(Mpc/h)=',s%box
      print*,'| rank       =',s%rank
      print*,'| nn         =',s%nn
      print*,'| nnt        =',s%nnt
      print*,'| nt         =',s%nt, ' ( nf_tile=',int(ncell*(nt+2*ncb),2),')'
      print*,'| ncell      =',s%ncell
      print*,'| ncb        =',s%ncb
      print*,'| izip x,v =',s%izipx,s%izipv
      print*,'| '
      print*,'| h0         =',s%h0
      print*,'| omega_m    =',s%omega_m
      print*,'| omega_l    =',s%omega_l
      print*,'| sigma_8    =',s%s8
      print*,'| m_neu/eV   =',s%m_neu
      print*,'| vsim2phys  =',s%vsim2phys, '(km/s)/(1.0)'
      print*,'| z_i        =',s%z_i
      print*,'------------------------------------------------------------------------------'
      endif
      sync all
    endsubroutine

    function image1d(cx,cy,cz)
      integer image1d,cx,cy,cz
      image1d=cx+nn*(cy-1)+nn**2*(cz-1)
    endfunction

    pure function image2str(nimage)
      character(:),allocatable :: image2str
      character(20) :: str
      integer,intent(in) :: nimage
      write(str,'(i6)') nimage
      image2str=trim(adjustl(str))
    endfunction

    pure function z2str(z)
      character(:),allocatable :: z2str
      character(20) :: str
      real,intent(in) :: z
      write(str,'(f7.3)') z
      z2str=trim(adjustl(str))
    endfunction

    function output_dir()
      character(:),allocatable :: output_dir
      character(20) :: str_z,str_i
      write(str_i,'(i6)') rank
      write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)
      output_dir=opath//'node'//trim(adjustl(str_i))//'/'
    endfunction

    function output_prefix()
      character(:),allocatable :: output_prefix
      character(20) :: str_z,str_i
      write(str_i,'(i6)') rank
      write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)
      output_prefix=opath//'node'//trim(adjustl(str_i))//'/'//trim(adjustl(str_z))
    endfunction

    function output_suffix()
      character(:),allocatable :: output_suffix
      character(20) :: str_i
      write(str_i,'(i6)') rank
      output_suffix='_'//trim(adjustl(str_i))//'.bin'
    endfunction

    function output_name(zipname)
      character(*) ::  zipname
      character(:),allocatable :: output_name
      output_name=output_prefix()//zipname//output_suffix()
    endfunction

    function ic_name(zipname)
      character(*) ::  zipname
      character(:),allocatable :: ic_name
      character(20) :: str_z,str_i
      write(str_i,'(i6)') rank
      write(str_z,'(f7.3)') z_i
      ic_name=opath//'node'//trim(adjustl(str_i))//'/'//trim(adjustl(str_z))//zipname//output_suffix()
    endfunction
endmodule
