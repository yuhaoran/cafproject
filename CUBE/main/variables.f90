!#define zipconvert
module variables
use parameters
implicit none
save

! parameters
integer,parameter :: npnode=nf**3
real,parameter :: density_buffer=2.0
integer,parameter :: npmax=npnode*(nte*1./nt)**3*density_buffer
integer,parameter ::  nseedmax=200
real,parameter :: vbuf=0.9
real,parameter :: dt_max=1
real,parameter :: dt_scale=1
integer(8),parameter :: unit8=1
integer,parameter :: NULL=0
real,parameter :: GG=1.0/6.0/pi

! variables
integer its[*]
real dt[*],dt_old[*],dt_mid[*]
real dt_fine(nn**3),dt_pp(nn**3),dt_coarse(nn**3),dt_vmax(nn**3)
real a[*],da[*],a_mid[*],tau[*],t[*] ! time step
real f2_max_fine(nnt,nnt,nnt)[*],f2_max_pp(nnt,nnt,nnt)[*],f2_max_coarse[*]

integer iseed(nseedmax), iseedsize
integer itx,ity,itz,ix,iy,iz,i_dim
integer i,j,k,l,ip,ipp,pp
integer nplocal[*], nptile(nnt,nnt,nnt)
integer(8) nptotal, npcheck

real mass_p

! FFT plans
integer(8) plan_fft_fine,plan_ifft_fine

real v_i2r(3)[*],v_i2r_new(3)[*]
real vmax(3)[*],vmax_new(3)[*]
! n^3
#ifdef zipconvert
  integer(1) xic_new(3,npmax)
  integer(2) vic_new(3,npmax)
#endif
integer(izipx) x(3,npmax)[*], x_new(3,npmax/nnt**3)
integer(izipv) v(3,npmax)[*], v_new(3,npmax/nnt**3)
#ifdef PID
  integer(2) pid(4,npmax)[*], pid_new(4,npmax/nnt**3)
#endif
integer(1) rhoc_i1(nt,nt,nt,nnt,nnt,nnt)
integer(4) rhoc_i4(nc**2)

real rho_f(nfe+2,nfe,nfe,ncore)
real crho_f(nfe+2,nfe,nfe,ncore)
real kern_f(3,nfe/2+1,nfe,nfe)
real force_f(3,nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1,ncore)
!integer rhoce1d(nce**3), rhoequiv(nce,nce,nce) ! rhoce is a coarray
!equivalence(rhoequiv,rhoce1d)

! rho in physical tiles and 6 buffers of tiles
! n^3
!integer rhotile(nt,nt,nt,nnt,nnt,nnt)[*]
integer rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
integer cum(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]


! coarse kernel arrays
real ck(3,nc,nc,nc)
!real kern_c(3,nc*nn/2,nc,npen+2)
real kern_c(3,nc*nn/2+1,nc,npen)
!real tmp_kern_c(3,nc*nn,nc,npen+2)
real tmp_kern_c(3,nc*nn+2,nc,npen)

real force_c(3,0:nc+1,0:nc+1,0:nc+1)[*]

character (10) :: img_s, z_s

!equivalence(rhoce,rhoce1d)

contains

  function cumsum3(input)
    implicit none
    integer input(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer cumsum3(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer nsum,igx,igy,igz
    nsum=0
    do igz=1-2*ncb,nt+2*ncb
    do igy=1-2*ncb,nt+2*ncb
    do igx=1-2*ncb,nt+2*ncb
    	nsum=nsum+input(igx,igy,igz)
    	cumsum3(igx,igy,igz)=nsum
    enddo
    enddo
    enddo
  endfunction cumsum3

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
