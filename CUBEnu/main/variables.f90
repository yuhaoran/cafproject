module variables
  use omp_lib
  use parameters
  implicit none
  save

  ! parameters
  integer(8),parameter :: np_image=(nc*np_nc)**3*merge(2,1,body_centered_cubic) ! average number of particles per image
  integer(8),parameter :: np_image_max=np_image*(nte*1./nt)**3*image_buffer
  integer(8),parameter :: np_tile_max=np_image/nnt**3*(nte*1./nt)**3*tile_buffer
  integer(8),parameter ::  nseedmax=200
  integer(8),parameter :: unit8=1
  real,parameter :: dt_max=1
  real,parameter :: dt_scale=1
  real,parameter :: GG=1.0/6.0/pi

  logical neutrino_flag
  integer n_checkpoint_neu ! the checkpoint to switch on neutrinos

  ! variables
  integer(8) istep
  real dt[*],dt_old[*],dt_mid[*],dt_e
  real a[*],da[*],a_mid[*],tau[*],t[*] ! time step
  real f2_max_fine(nnt,nnt,nnt)[*],f2_max_pp[*],f2_max_coarse[*]

  integer(4) iseed(nseedmax), iseedsize, nth,ith
  integer(8) itx,ity,itz,ix,iy,iz,i_dim
  integer(8) i,j,k,l,ip,ipp,pp,nzero
  integer(8) nptile(nnt,nnt,nnt), npcheck
  real(8) xq(3),deltax(3),deltav(3),vreal(3)

  ! FFT plans
  integer(8) plan_fft_fine,plan_ifft_fine
  real vmax,vmax_nu,overhead_tile[*],overhead_image[*]
  real sigma_vi,sigma_vi_new,sigma_vi_nu,sigma_vi_new_nu
  !real vdisp(506,2),sigma_vi_old,sigma_vi
  real(4) svz(500,2),svr(100,2)
  real(8) sigma_vci,sigma_vfi,sigma_vres,sigma_vci_old,sigma_vfi_old,sigma_vres_old
  real(8) std_vsim_c[*],std_vsim_res[*],std_vsim[*]
  real(8) std_vsim_c_nu[*],std_vsim_res_nu[*],std_vsim_nu[*]
  ! n^3
  integer(izipx) xp(3,np_image_max)[*], xp_new(3,np_tile_max)
  integer(izipv) vp(3,np_image_max)[*], vp_new(3,np_tile_max)
#ifdef PID
    integer(8) pid(np_image_max)[*], pid_new(np_tile_max)
#endif
#ifdef AFIELD
    real(4) afield(3,np_image_max)
#endif

  real rho_f(nfe+2,nfe,nfe)
  real crho_f(nfe+2,nfe,nfe)
  real kern_f(nfe/2+1,nfe,nfe,3)
  real force_f(3,nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1)

  integer(4) rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  real(4) vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt) ! cannot have >7 dims
  !integer(8) cum(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  !integer(8) cum5(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]

  integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),codimension[*] :: idx_b_l,idx_b_r
  integer(8),dimension(nt,nt,nnt,nnt,nnt),codimension[*] :: ppl0,pplr,pprl,ppr0,ppl,ppr

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
  integer(8),parameter :: np_pp_max=np_image/nnt**3*(1+2./nt)**3*tile_buffer
  integer hoc(1-ncell:nft+ncell,1-ncell:nft+ncell,1-ncell:nft+ncell)
  integer ll(np_pp_max)
  integer(8) npairs,itest1,nlast
  real(8) xvec1(3),xvec2(3),xvec21(3),rmag,force_pp(3),rcut,pcut,f_tot(3)
  integer(4) ivec1(3),np,ii,jj,kk,np1,np2,l1,l2
  integer igx,igy,igz,lp,ip1,ip2
  integer t1,t2,tt1,tt2,ttt1,ttt2,t_rate

contains

  !function cumsum3(rho_input)
  !  implicit none
  !  integer(4) rho_input(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  !  integer(8) cumsum3(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  !  integer(8) nsum,igx,igy,igz
  !  nsum=0
  !  do igz=1-2*ncb,nt+2*ncb
  !  do igy=1-2*ncb,nt+2*ncb
  !  do igx=1-2*ncb,nt+2*ncb
  !    nsum=nsum+rho_input(igx,igy,igz)
  !    cumsum3(igx,igy,igz)=nsum
  !  enddo
  !  enddo
  !  enddo
  !endfunction cumsum3

  !pure function spine3(rho) ! cumulative sum record at the yz-plane
  !  implicit none
  !  integer(4),intent(in) :: rho(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  !  integer(8) spine3(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  !  integer(8) nsum,igz,igy
  !  nsum=0
  !  do igz=1-2*ncb,nt+2*ncb
  !  do igy=1-2*ncb,nt+2*ncb
  !    nsum=nsum+sum(rho(:,igy,igz))
  !    spine3(igy,igz)=nsum
  !  enddo
  !  enddo
  !endfunction

  pure subroutine spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)
    implicit none
    integer(4),intent(in) :: rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8),dimension(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb),intent(out) :: idx_ex_r
    integer(8),dimension(nt,nt),intent(out) :: pp_l,pp_r,ppe_l,ppe_r
    integer(8) nsum,np_phy
    integer igz,igy
    ! spine := yz-plane to record cumulative sum
    nsum=0
    do igz=1-2*ncb,nt+2*ncb
    do igy=1-2*ncb,nt+2*ncb
      nsum=nsum+sum(rhoce(:,igy,igz))
      idx_ex_r(igy,igz)=nsum ! right index
    enddo
    enddo
    do igz=1,nt
    do igy=1,nt
      ppe_r(igy,igz)=idx_ex_r(igy,igz)-sum(rhoce(nt+1:,igy,igz))
    enddo
    enddo
    nsum=0
    do igz=1,nt
    do igy=1,nt
      pp_l(igy,igz)=nsum+1
      np_phy=sum(rhoce(1:nt,igy,igz))
      nsum=nsum+np_phy
      pp_r(igy,igz)=nsum
      ppe_l(igy,igz)=ppe_r(igy,igz)-np_phy+1
    enddo
    enddo
  endsubroutine

  pure subroutine spine_image(rhoc,idx_b_l,idx_b_r,ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr)
    implicit none
    integer(4),intent(in) :: rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),intent(out) :: idx_b_l,idx_b_r
    integer(8),dimension(nt,nt,nnt,nnt,nnt),intent(out) :: ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr
    integer(8) nsum,nsum_p,np_phy,ihz,ihy,ihx,igz,igy,ctile_mass(nnt,nnt,nnt),ctile_mass_p(nnt,nnt,nnt)
    ! spine := yz-plane to record cumulative sum
    nsum=0;nsum_p=0
    do ihz=1,nnt ! sum cumulative tile mass first
    do ihy=1,nnt
    do ihx=1,nnt
      ctile_mass(ihx,ihy,ihz)=nsum
      ctile_mass_p(ihx,ihy,ihz)=nsum_p
      nsum=nsum+sum(rhoc(:,:,:,ihx,ihy,ihz))
      nsum_p=nsum_p+sum(rhoc(1:nt,1:nt,1:nt,ihx,ihy,ihz))
    enddo
    enddo
    enddo
    do ihz=1,nnt ! calculate extended spine cumulative index on both sides
    do ihy=1,nnt
    do ihx=1,nnt
      nsum=ctile_mass(ihx,ihy,ihz)
      do igz=1-ncb,nt+ncb
      do igy=1-ncb,nt+ncb
        idx_b_l(igy,igz,ihx,ihy,ihz)=nsum
        nsum=nsum+sum(rhoc(:,igy,igz,ihx,ihy,ihz))
        idx_b_r(igy,igz,ihx,ihy,ihz)=nsum
      enddo
      enddo
    enddo
    enddo
    enddo
    do ihz=1,nnt ! calculate physical spine
    do ihy=1,nnt
    do ihx=1,nnt
      nsum_p=ctile_mass_p(ihx,ihy,ihz)
      do igz=1,nt
      do igy=1,nt
        ppl(igy,igz,ihx,ihy,ihz)=nsum_p
        np_phy=sum(rhoc(1:nt,igy,igz,ihx,ihy,ihz))
        nsum_p=nsum_p+np_phy
        ppr(igy,igz,ihx,ihy,ihz)=nsum_p

        ppe_r0(igy,igz,ihx,ihy,ihz)=idx_b_r(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt+1:,igy,igz,ihx,ihy,ihz))
        ppe_rl(igy,igz,ihx,ihy,ihz)= ppe_r0(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt-ncb+1:nt,igy,igz,ihx,ihy,ihz))

        ppe_l0(igy,igz,ihx,ihy,ihz)=idx_b_l(igy,igz,ihx,ihy,ihz)+sum(rhoc(:0,igy,igz,ihx,ihy,ihz))
        ppe_lr(igy,igz,ihx,ihy,ihz)=ppe_l0(igy,igz,ihx,ihy,ihz) +sum(rhoc(1:ncb,igy,igz,ihx,ihy,ihz))
      enddo
      enddo
    enddo
    enddo
    enddo
  endsubroutine

#ifdef OLD
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

  pure function spine6(rho) ! cumulative sum record at the (yz,itx,ity,itz)-plane
    implicit none
    integer(4),intent(in) :: rho(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8) spine6(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8) nsum,ihz,ihy,ihx,igz,igy
    nsum=0
    do ihz=1,nnt
    do ihy=1,nnt
    do ihx=1,nnt
      do igz=1-ncb,nt+ncb
      do igy=1-ncb,nt+ncb
        nsum=nsum+sum(rho(:,igy,igz,ihx,ihy,ihz))
        spine6(igy,igz,ihx,ihy,ihz)=nsum
      enddo
      enddo
    enddo
    enddo
    enddo
  endfunction
#endif

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
