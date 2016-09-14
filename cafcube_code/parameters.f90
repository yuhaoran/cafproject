module parameters
implicit none
save
!parameters

character(*),parameter :: ipath='./init/2lpt/'
character(*),parameter :: opath='./output/'

integer,parameter :: izipx=2
integer,parameter :: izipv=2
integer(izipx),parameter :: ishift=-(2**(izipx*8-1))
real,parameter :: rshift=0.5-ishift

integer,parameter :: nn=1
integer,parameter :: ncell=4
integer,parameter :: rsoft=0.1

integer,parameter :: nnt=2 ! n_tile/node/dim
integer,parameter :: nc=64 ! nc_physical/node/dim
integer,parameter :: nt=nc/nnt ! nc_physical/tile/dim
integer,parameter :: npen=nc/nn ! pencil

integer,parameter :: nf=nc*ncell ! 96
integer,parameter :: nft=nt*ncell ! 48
integer,parameter :: nfpen=npen*ncell

integer,parameter :: ncore=1
integer,parameter :: n_nest=4

integer,parameter :: ncb=6 ! 6 by default
integer,parameter :: nce=nc+2*ncb
integer,parameter :: nte=nt+2*ncb

integer,parameter :: nfb=ncb*ncell ! 24
integer,parameter :: nf_cutoff=16

integer,parameter :: nfe=nft+2*nfb ! 96

real,parameter :: z_i=50
real,parameter :: z_i_nu=5
real,parameter :: a_i=1/(1+z_i)

real,parameter :: box=300 ! Mpc/h
real,parameter :: h0=68
real,parameter :: ratio_nudm_dim=2
real,parameter :: m_neu=0.05
real,parameter :: omega_nu=m_neu/93.14/(h0/100.)**2

!real,parameter :: omega_c=0.27
!real,parameter :: omega_b=0.05
!real,parameter :: omega_m=omega_c+omega_b+omega_nu
!real,parameter :: omega_l=1-omega_m
real,parameter :: omega_l=0.73
real,parameter :: omega_m=1-omega_l

real,parameter :: omega_ch=0.7
real,parameter :: bias=1
real,parameter :: power_index=2
real,parameter :: wde=-1

real,parameter :: f_nl=0
real,parameter :: g_nl=0
real,parameter :: n_s=0.96
real,parameter :: scalar_amp=2.46e-9

integer,parameter :: nts=400 ! maximum number of timesteps
real,parameter :: ra_max=0.05
real,parameter :: v_resolution=2.1/(2**(izipv*8))
real,parameter :: x_resolution=1.0/2**(izipx*8)

!! MPI images !!
integer,parameter :: rank=this_image()-1            ! MPI_rank
integer,parameter :: icz=rank/(nn**2)+1             ! image_z
integer,parameter :: icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
integer,parameter :: icx=mod(rank,nn)+1             ! image_x
! adjacent images
integer,parameter :: inx=modulo(icx-2,nn)+1
integer,parameter :: iny=modulo(icy-2,nn)+1
integer,parameter :: inz=modulo(icz-2,nn)+1
integer,parameter :: ipx=modulo(icx,nn)+1
integer,parameter :: ipy=modulo(icy,nn)+1
integer,parameter :: ipz=modulo(icz,nn)+1

logical,parameter :: head=(rank==0)

! 128 byte (equivalent 32 4-byte variables) header in zip2
type sim_header
  ! standard cubep3m 18 variables
  integer nplocal
  real a, t, tau
  integer nts
  real dt_f_acc, dt_pp_acc, dt_c_acc
  integer cur_checkpoint,cur_proj,cur_halo
  real mass_p
  real v_r2i(3)
  real shake_offset(3) ! -18
  ! more simulation config info
  real box
  integer(4) rank
  integer(1) nn,nnt,nt,ncell,ncb
  integer(1) izipx,izipv,izip2
  logical(4) tile_struc ! 22
  real h0
  real omega_m
  real omega_l
  real s8
  real m_neu(3)
  real vsim2phys
  real z_i
endtype

type(sim_header) sim

contains

subroutine print_header(s)
  type(sim_header),intent(in) :: s
  print*,'-------------------------------- CUBEP3M info --------------------------------'
  print*,'| nplocal      =',s%nplocal
  print*,'| a,t,tau      =',s%a,s%t,s%tau
  print*,'| nts          =',s%nts
  print*,'| dt f,pp,c    =',s%dt_f_acc,s%dt_pp_acc,s%dt_c_acc
  print*,'| cur_steps    =',int(s%cur_checkpoint,2),int(s%cur_proj,2),int(s%cur_halo,2)
  print*,'| mass_p       =',s%mass_p
  print*,'| v_r2i        =',s%v_r2i
  print*,'| shake_offset =',s%shake_offset
  print*,'| '
  print*,'| box/(Mpc/h)=',s%box
  print*,'| rank       =',int(s%rank,2)
  print*,'| nn         =',s%nn
  print*,'| nnt        =',s%nnt
  print*,'| nt         =',s%nt, ' ( nf_tile=',int(ncell*(nt+2*ncb),2),')'
  print*,'| ncell      =',s%ncell
  print*,'| ncb        =',s%ncb
  print*,'| izip x,v,2 =',s%izipx,s%izipv,s%izip2
  print*,'| tile_struc =',s%tile_struc
  print*,'| h0         =',s%h0
  print*,'| omega_m    =',s%omega_m
  print*,'| omega_l    =',s%omega_l
  print*,'| sigma_8    =',s%s8
  print*,'| m_neu/eV   =',s%m_neu
  print*,'| vsim2phys  =',s%vsim2phys
  print*,'| z_i        =',s%z_i
  print*,'------------------------------------------------------------------------------'
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

function cumsum6(input)
implicit none
integer input(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
integer cumsum6(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
integer nsum,ihx,ihy,ihz,igx,igy,igz
nsum=0
do ihz=1,nnt
do ihy=1,nnt
do ihx=1,nnt
do igz=1-ncb,nt+ncb
do igy=1-ncb,nt+ncb
do igx=1-ncb,nt+ncb
  nsum=nsum+input(igx,igy,igz,ihx,ihy,ihz)
  cumsum6(igx,igy,igz,ihx,ihy,ihz)=nsum
enddo
enddo
enddo
enddo
enddo
enddo
endfunction

endmodule
