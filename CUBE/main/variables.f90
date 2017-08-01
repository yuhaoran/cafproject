module variables
  use parameters
  implicit none
  save

  ! parameters
  integer(8),parameter :: np_image=(nc*np_nc)**3*merge(2,1,np_2n3) ! average number of particles per image
  integer(8),parameter :: np_image_max=np_image*(nte*1./nt)**3*image_buffer
  integer(8),parameter :: np_tile_max=np_image/nnt**3*(nte*1./nt)**3*tile_buffer
  integer(8),parameter ::  nseedmax=200
  integer(8),parameter :: unit8=1
  real,parameter :: vbuf=0.9
  real,parameter :: dt_max=1
  real,parameter :: dt_scale=1
  real,parameter :: GG=1.0/6.0/pi

  ! variables
  integer(8) istep
  real dt[*],dt_old[*],dt_mid[*]
  real dt_fine[*],dt_pp[*],dt_coarse[*],dt_vmax[*]
  real a[*],da[*],a_mid[*],tau[*],t[*] ! time step
  real f2_max_fine(nnt,nnt,nnt)[*],f2_max_pp(nnt,nnt,nnt)[*],f2_max_coarse[*]

  integer(4) iseed(nseedmax), iseedsize
  integer(8) itx,ity,itz,ix,iy,iz,i_dim
  integer(8) i,j,k,l,ip,ipp,pp
  integer(8) nplocal[*], nptile(nnt,nnt,nnt)
  integer(8) npglobal, npcheck
  real(8) xq(3),deltax(3),deltav(3),vreal(3)

  real mass_p

  ! FFT plans
  integer(8) plan_fft_fine,plan_ifft_fine
  real vmax,overhead_tile[*],overhead_image[*],sigma_vi,sigma_vi_new
  !real vdisp(506,2),sigma_vi_old,sigma_vi
  real(4) svz(500,2),svr(100,2)
  real(8) sigma_vci,sigma_vfi,sigma_vres,sigma_vci_old,sigma_vfi_old,sigma_vres_old
  real(8) std_vsim_c[*],std_vsim_res[*],std_vsim[*]
  ! n^3
  integer(izipx) xp(3,np_image_max)[*], xp_new(3,np_tile_max)
  integer(izipv) vp(3,np_image_max)[*], vp_new(3,np_tile_max)
#ifdef PID
    integer(8) pid(np_image_max)[*], pid_new(np_tile_max)
#endif

  real rho_f(nfe+2,nfe,nfe)
  real crho_f(nfe+2,nfe,nfe)
  real kern_f(nfe/2+1,nfe,nfe,3)
  real force_f(3,nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1)

  integer(4) rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  real(4) vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt) ! cannot have >7 dims
  integer(8) cum(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]

  ! the following variables are introduced because
  ! gcc only allows <= 7 ranks in arrays
  real(4) vtransx(3,ncb,nt+2*ncb,nt+2*ncb,nnt,nnt)[*]
  real(4) vtransy(3,nt+2*ncb,ncb,nt+2*ncb,nnt,nnt)[*]
  real(4) vtransz(3,nt+2*ncb,nt+2*ncb,ncb,nnt,nnt)[*]

  ! coarse kernel arrays
  real ck(3,nc,nc,nc)
  real kern_c(nc*nn/2+1,nc,npen,3)
  real crho_c(nc*nn+2,nc,npen) !!! temp
  real force_c(3,0:nc+1,0:nc+1,0:nc+1)[*]

  character (10) :: img_s, z_s

  !equivalence(rhoce,rhoce1d)

contains

  function cumsum3(rho_input)
    implicit none
    integer(4) rho_input(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8) cumsum3(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8) nsum,igx,igy,igz
    nsum=0
    do igz=1-2*ncb,nt+2*ncb
    do igy=1-2*ncb,nt+2*ncb
    do igx=1-2*ncb,nt+2*ncb
      nsum=nsum+rho_input(igx,igy,igz)
      cumsum3(igx,igy,igz)=nsum
    enddo
    enddo
    enddo
  endfunction cumsum3

  function cumsum6(rho_input)
    implicit none
    integer(4) rho_input(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8) cumsum6(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8) nsum,ihx,ihy,ihz,igx,igy,igz
    nsum=0
    do ihz=1,nnt
    do ihy=1,nnt
    do ihx=1,nnt
    do igz=1-ncb,nt+ncb
    do igy=1-ncb,nt+ncb
    do igx=1-ncb,nt+ncb
      nsum=nsum+rho_input(igx,igy,igz,ihx,ihy,ihz)
      cumsum6(igx,igy,igz,ihx,ihy,ihz)=nsum
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  endfunction

  real function interp_sigmav(aa,rr)
    implicit none
    integer(8) ii,i1,i2
    real aa,rr,term_z,term_r
    i1=1
    i2=500
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (aa>svz(ii,1)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    term_z=svz(i1,2)+(svz(i2,2)-svz(i1,2))*(aa-svz(i1,1))/(svz(i2,1)-svz(i1,1))
    i1=1
    i2=100
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (rr>svz(ii,1)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    term_r=svr(i1,2)+(svr(i2,2)-svr(i1,2))*(rr-svr(i1,1))/(svr(i2,1)-svr(i1,1))
    interp_sigmav=term_z*term_r
    print*,term_z,term_r
  endfunction

endmodule
