subroutine particle_mesh
  use omp_lib
  use variables
  use cubefft
  use pencil_fft
  implicit none
  save

  ! force settings
  !logical,parameter :: fine_force=.true.
  !logical,parameter :: coarse_force=.true.
  !logical,parameter :: pp_force=.false.
  !logical,parameter :: ext_pp_force=.false.

  integer(8) nlast,nlen
  integer(4) ithread, nthread
  integer(8) idxf(3),np
  integer(8) idx1(3), idx2(3)
  real tempx(3), dx1(3), dx2(3)
  real r3t(-1:nt+2,-1:nt+2,-1:nt+2) ! coarse density on tile, with buffer=2

  if (head) then
    print*, ''
    print*, 'particle mesh'
  endif

  cum=cumsum6(rhoc)
  !nthread=1
  !ithread=omp_get_thread_num()+1
  !print*,'ithread'

  vmax=0
  f2_max_fine(1:nnt,1:nnt,1:nnt)=0
  f2_max_pp=0
  f2_max_coarse=0

  if (head) print*, '  pm fine over',nnt**3,'tiles'

  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    !if (head) print*,'    tile',int(itx,1),int(ity,1),int(itz,1)
    ! fine_cic_mass ------------------------------------------------------------
    rho_f=0
    crho_f=0
    !if (head) print*,'      fine_cic_mass'
    do k=2-ncb,nt+ncb-1
    do j=2-ncb,nt+ncb-1
    do i=2-ncb,nt+ncb-1
      nlast=cum(i-1,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np ! loop over particle
        ip=nlast+l
        tempx=4.*((/i,j,k/)-1)+4*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution !-0.5
        idx1 = floor(tempx) + 1
        idx2 = idx1 + 1
        dx1 = idx1 - tempx
        dx2 = 1 - dx1
        idx1=idx1+nfb
        idx2=idx2+nfb
        rho_f(idx1(1),idx1(2),idx1(3))=rho_f(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
        rho_f(idx2(1),idx1(2),idx1(3))=rho_f(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
        rho_f(idx1(1),idx2(2),idx1(3))=rho_f(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
        rho_f(idx1(1),idx1(2),idx2(3))=rho_f(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
        rho_f(idx1(1),idx2(2),idx2(3))=rho_f(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
        rho_f(idx2(1),idx1(2),idx2(3))=rho_f(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
        rho_f(idx2(1),idx2(2),idx1(3))=rho_f(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
        rho_f(idx2(1),idx2(2),idx2(3))=rho_f(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p
      enddo
    enddo
    enddo
    enddo
    ! fine force ---------------------------------------------------------------
    !if (head) print*,'      fine_fft'
    call sfftw_execute(plan_fft_fine)
    crho_f(:,:,:)=rho_f(:,:,:) ! back up
    do i_dim=1,3
      !if (head) print*,'      fine_ifft dim',int(i_dim,1)
      rho_f(::2,:,:)=-crho_f(2::2,:,:)*kern_f(:,:,:,i_dim)
      rho_f(2::2,:,:)=crho_f(::2,:,:)*kern_f(:,:,:,i_dim)
      call sfftw_execute(plan_ifft_fine)
      rho_f=rho_f/real(nfe)/real(nfe)/real(nfe)
      force_f(i_dim,:,:,:)=rho_f(nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1)
    enddo
    f2_max_fine(itx,ity,itz)=maxval(sum(force_f(:,:,:,:)**2,1))
    ! fine velocity ------------------------------------------------------------
    !if (head) print*,'      fine velocity'
    do k=1,nt
    do j=1,nt
    do i=1,nt ! loop over coarse cell
      nlast=cum(i-1,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np ! loop over particle
        ip=nlast+l
        tempx=4.*((/i,j,k/)-1)+4*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution !-0.5
        idx1 = floor(tempx) + 1
        idx2 = idx1 + 1
        dx1 = idx1 - tempx
        dx2 = 1 - dx1
        idx1=idx1+nfb
        idx2=idx2+nfb
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))

        vreal=vreal+force_f(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
        vreal=vreal+force_f(:,idx2(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
        vreal=vreal+force_f(:,idx1(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
        vreal=vreal+force_f(:,idx1(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
        vreal=vreal+force_f(:,idx1(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
        vreal=vreal+force_f(:,idx2(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
        vreal=vreal+force_f(:,idx2(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
        vreal=vreal+force_f(:,idx2(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)

        vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi_new*vrel_boost)*vreal)/pi,kind=izipv)
      enddo

    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  sigma_vi=sigma_vi_new
  sync all
  !-----------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (head) print*, '  pm coarse'
  ! coarse_cic_mass ------------------------------------------------------------
  if (head) print*, '    coarse cic mass'
  r3=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt ! loop over tile
    r3t=0
    do k=0,nt+1
    do j=0,nt+1
    do i=0,nt+1
      nlast=cum(i-1,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np ! loop over particle
        ip=nlast+l
        tempx=((/i,j,k/)-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution-0.5
        idx1(:)=floor(tempx(:))+1
        idx2(:)=idx1(:)+1
        dx1(:)=idx1(:)-tempx(:) ! CIC contribution to idx1
        dx2(:)=1-dx1(:) ! CIC contribution to idx2
        r3t(idx1(1),idx1(2),idx1(3))=r3t(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
        r3t(idx2(1),idx1(2),idx1(3))=r3t(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
        r3t(idx1(1),idx2(2),idx1(3))=r3t(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
        r3t(idx1(1),idx1(2),idx2(3))=r3t(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
        r3t(idx1(1),idx2(2),idx2(3))=r3t(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
        r3t(idx2(1),idx1(2),idx2(3))=r3t(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
        r3t(idx2(1),idx2(2),idx1(3))=r3t(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
        r3t(idx2(1),idx2(2),idx2(3))=r3t(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p
      enddo
    enddo
    enddo
    enddo
    ! put center part of r3t into subset of r3
    r3((itx-1)*nt+1:itx*nt,(ity-1)*nt+1:ity*nt,(itz-1)*nt+1:itz*nt)=r3t(1:nt,1:nt,1:nt)
  enddo
  enddo
  enddo
  sync all
  ! coarse force ---------------------------------------------------------------
  if (head) print*, '    coarse cic force'
  if (head) print*,'      coarse_fft'
  call pencil_fft_forward
  ! save complex rho_c into crho_c
  crho_c(::2,:,:)=real(cxyz)
  crho_c(2::2,:,:)=imag(cxyz)
  do i_dim=1,3
    if (head) print*,'      coarse_ifft dim',int(i_dim,1)
    rxyz(::2,:,:)=-crho_c(2::2,:,:)*kern_c(:,:,:,i_dim)
    rxyz(2::2,:,:)=crho_c(::2,:,:)*kern_c(:,:,:,i_dim)
    call pencil_fft_backward
    force_c(i_dim,1:nc,1:nc,1:nc)=r3
  enddo
  sync all
  ! sync force_c buffer for CIC force
  if (head) print*, '      sync force_c buffer'
  force_c(:,0,:,:)=force_c(:,nc,:,:)[image1d(inx,icy,icz)]
  force_c(:,nc+1,:,:)=force_c(:,1,:,:)[image1d(ipx,icy,icz)]
  sync all
  force_c(:,:,0,:)=force_c(:,:,nc,:)[image1d(icx,iny,icz)]
  force_c(:,:,nc+1,:)=force_c(:,:,1,:)[image1d(icx,ipy,icz)]
  sync all
  force_c(:,:,:,0)=force_c(:,:,:,nc)[image1d(icx,icy,inz)]
  force_c(:,:,:,nc+1)=force_c(:,:,:,1)[image1d(icx,icy,ipz)]
  sync all
  ! coarse_max_dt
  f2_max_coarse=maxval(sum(force_c**2,1))
  sync all
  ! coarse velocity ------------------------------------------------------------
  if (head) print*, '    coarse cic velocity'
  do itz=1,nnt ! loop again
  do ity=1,nnt
  do itx=1,nnt ! loop over tiles
    do k=1,nt
    do j=1,nt
    do i=1,nt
      nlast=cum(i-1,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np ! loop over particle
        ip=nlast+l
        tempx=((/itx,ity,itz/)-1)*nt+((/i,j,k/)-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution-0.5
        idx1(:)=floor(tempx(:))+1
        idx2(:)=idx1(:)+1
        dx1(:)=idx1(:)-tempx(:)
        dx2(:)=1-dx1(:)
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        vreal=vreal+force_c(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
        vreal=vreal+force_c(:,idx2(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
        vreal=vreal+force_c(:,idx1(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
        vreal=vreal+force_c(:,idx1(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
        vreal=vreal+force_c(:,idx1(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
        vreal=vreal+force_c(:,idx2(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
        vreal=vreal+force_c(:,idx2(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
        vreal=vreal+force_c(:,idx2(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)
        vmax=max(vmax,maxval(vreal+vfield(:,i,j,k,itx,ity,itz)))
        vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
      enddo
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  sim%vsim2phys=(1.5/a)*box*h0*sqrt(omega_m)/nf_global
  sync all

  if (head) print*, '  constrain dt'
  dt_fine=sqrt( 1.0 / (sqrt(maxval(f2_max_fine))*a_mid*GG) )
  dt_coarse=sqrt( real(ncell) / (sqrt(f2_max_coarse)*a_mid*GG) )
  dt_pp=1000
  dt_vmax=vbuf*20/vmax
  sync all

  do i=1,nn**3
    dt_fine=min(dt_fine,dt_fine[i])
    dt_coarse=min(dt_coarse,dt_coarse[i])
    dt_pp=min(dt_pp,dt_pp[i])
    dt_vmax=min(dt_vmax,dt_vmax[i])
  enddo
  sync all

endsubroutine
