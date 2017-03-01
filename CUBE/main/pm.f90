!#define debug_force
!#define write_file
!#define finengp
!#define coarsengp
module pm
contains

subroutine particle_mesh
use omp_lib
use variables
use cubefft
use pencil_fft
implicit none
save

! force settings
logical,parameter :: fine_force=.true.
logical,parameter :: coarse_force=.true.
logical,parameter :: pp_force=.false.
logical,parameter :: ext_pp_force=.false.

integer idx1min(3),idx1max(3),idx2min(3),idx2max(3) !!! debug
integer nlast,nlen
integer ithread, nthread
integer idxf(3),np
integer idx1(3), idx2(3)
real vtemp(3), tempx(3), dx1(3), dx2(3)
real r3t(-1:nt+2,-1:nt+2,-1:nt+2) ! coarse density on tile, with buffer=2

if (head) then
  print*, ''
  print*, 'particle mesh'
endif

cum=cumsum6(rhoc)

nthread=1
!ithread=omp_get_thread_num()+1
ithread=1
!print*,'ithread',ithread

vmax_new=0
f2_max_fine=0
f2_max_pp=0
f2_max_coarse=0

if (fine_force) then
if (head) print*, '  pm fine'
if (head) print*, '    first loop'
do itz=1,nnt
do ity=1,nnt
do itx=1,nnt
  ! fine_cic_mass ---------------------------------------------------
  rho_f(:,:,:,ithread)=0
  crho_f(:,:,:,ithread)=0
  do k=2-ncb,nt+ncb-1
  do j=2-ncb,nt+ncb-1
  do i=2-ncb,nt+ncb-1
    nlast=cum(i-1,j,k,itx,ity,itz)
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np ! loop over particle
      ip=nlast+l
#ifndef finengp
      ! CIC
      !tempx = 4.*((/i,j,k/)-1) + 4*((x(:,ip)+int(-128,1))+128.5)/256 - 0.5
      tempx=4.*((/i,j,k/)-1)+4*(x(:,ip)+ishift+rshift)*x_resolution !-0.5
      idx1 = floor(tempx) + 1
      idx2 = idx1 + 1
      dx1 = idx1 - tempx
      dx2 = 1 - dx1
      idx1=idx1+nfb
      idx2=idx2+nfb ! positive integer indeces
      !idx1min=min(idx1,idx1min)
      !idx1max=max(idx1,idx1max)
      !idx2min=min(idx2,idx2min)
      !idx2max=max(idx2,idx2max)

      rho_f(idx1(1),idx1(2),idx1(3),ithread)=rho_f(idx1(1),idx1(2),idx1(3),ithread)+dx1(1)*dx1(2)*dx1(3)*mass_p
      rho_f(idx2(1),idx1(2),idx1(3),ithread)=rho_f(idx2(1),idx1(2),idx1(3),ithread)+dx2(1)*dx1(2)*dx1(3)*mass_p
      rho_f(idx1(1),idx2(2),idx1(3),ithread)=rho_f(idx1(1),idx2(2),idx1(3),ithread)+dx1(1)*dx2(2)*dx1(3)*mass_p
      rho_f(idx1(1),idx1(2),idx2(3),ithread)=rho_f(idx1(1),idx1(2),idx2(3),ithread)+dx1(1)*dx1(2)*dx2(3)*mass_p
      rho_f(idx1(1),idx2(2),idx2(3),ithread)=rho_f(idx1(1),idx2(2),idx2(3),ithread)+dx1(1)*dx2(2)*dx2(3)*mass_p
      rho_f(idx2(1),idx1(2),idx2(3),ithread)=rho_f(idx2(1),idx1(2),idx2(3),ithread)+dx2(1)*dx1(2)*dx2(3)*mass_p
      rho_f(idx2(1),idx2(2),idx1(3),ithread)=rho_f(idx2(1),idx2(2),idx1(3),ithread)+dx2(1)*dx2(2)*dx1(3)*mass_p
      rho_f(idx2(1),idx2(2),idx2(3),ithread)=rho_f(idx2(1),idx2(2),idx2(3),ithread)+dx2(1)*dx2(2)*dx2(3)*mass_p
#else
      ! NGP
      idxf=ceiling( (4.*((/i,j,k/)-1) + 4*((x(:,ip)+int(-128,1))+128.5)/256) )
      idxf=idxf+nfb ! positive integer indeces
      rho_f(idxf(1),idxf(2),idxf(3),ithread)=rho_f(idxf(1),idxf(2),idxf(3),ithread)+mass_p
#endif
    enddo
  enddo
  enddo
  enddo

  !!! print*, 'real-space rho_f =', sum(rho_f*1.d0), maxval(rho_f)
#ifdef write_file
  open(10,file='testrhof.dat',status='replace',access='stream')
  write(10) rho_f(:nfe,:,:,1)
  close(10)
#endif

  ! fine force --------------------------------------------------------
  call sfftw_execute(plan_fft_fine)
!  print*, 'sum of Fourier-space rho_f =', sum(rho_f*1d0), maxval(rho_f)
  crho_f(:,:,:,ithread)=rho_f(:,:,:,ithread) ! back up
  do i_dim=1,3
    rho_f(::2,:,:,ithread)=-crho_f(2::2,:,:,ithread)*kern_f(i_dim,:,:,:)
    rho_f(2::2,:,:,ithread)=crho_f(::2,:,:,ithread)*kern_f(i_dim,:,:,:)
    call sfftw_execute(plan_ifft_fine)
    rho_f=rho_f/real(nfe)/real(nfe)/real(nfe)
    force_f(i_dim,:,:,:,ithread)=rho_f(nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1,ithread)
  enddo
  ! if (head) print*, '      max force_f', maxval(abs(force_f(1,:,:,:,1))),maxval(abs(force_f(2,:,:,:,1))),maxval(abs(force_f(3,:,:,:,1)))

  ! max force
  f2_max_fine(itx,ity,itz)=maxval(sum(force_f(:,:,:,:,ithread)**2,1))
  !print*, 'f2_max_fine(itx,ity,itz)', f2_max_fine(itx,ity,itz)
  ! fine max velocity

  !print*, sqrt(maxval(sum(force_f(:,:,:,:,ithread)**2,1)))

  ! fine velocity ---------------------------------------------------
  do k=1,nt
  do j=1,nt
  do i=1,nt ! loop over coarse cell
    nlast=cum(i-1,j,k,itx,ity,itz)
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np ! loop over particle
      ip=nlast+l
#ifndef finengp
      ! CIC
      tempx=ncell*((/i,j,k/)-1)+ncell*(x(:,ip)+ishift+rshift)*x_resolution !-0.5
      idx1 = floor(tempx) + 1
      idx2 = idx1 + 1
      dx1 = idx1 - tempx
      dx2 = 1 - dx1
      idx1=idx1+nfb
      idx2=idx2+nfb ! positive integer indeces

      vtemp=v(:,ip)*v_i2r
      vtemp=vtemp+force_f(:,idx1(1),idx1(2),idx1(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx2(1),idx1(2),idx1(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx1(1),idx2(2),idx1(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx1(1),idx1(2),idx2(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_f(:,idx1(1),idx2(2),idx2(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
      vtemp=vtemp+force_f(:,idx2(1),idx1(2),idx2(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_f(:,idx2(1),idx2(2),idx1(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx2(1),idx2(2),idx2(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)
#else
      ! NGP
      idxf=ceiling( (ncell*((/i,j,k/)-1) + ncell*((x(:,ip)+int(-128,1))+128.5)/256) )
      idxf=idxf+nfb ! positive integer indeces
      vtemp=v(:,ip)*v_i2r
      vtemp=vtemp+force_f(:,idxf(1),idxf(2),idxf(3),ithread)*a_mid*dt/6/pi
#endif
      vmax_new(1)=max(vmax_new(1),abs(vtemp(1)))
      vmax_new(2)=max(vmax_new(2),abs(vtemp(2)))
      vmax_new(3)=max(vmax_new(3),abs(vtemp(3)))
    enddo
  enddo
  enddo
  enddo

enddo
enddo
enddo
sync all

! sync vmax over images
if (head) print*, '    co_max vmax_new over images'
do i=1,nn**3
  vmax_new=max(vmax_new,vmax_new(:)[i])
enddo
v_i2r_new=vmax_new*v_resolution
if (head) print*, '    vmax_new =', vmax_new
sync all


if (head) print*, '    second loop'
do itz=1,nnt
do ity=1,nnt
do itx=1,nnt

  ! fine_cic_mass ------------------------------------------
  rho_f(:,:,:,ithread)=0
  crho_f(:,:,:,ithread)=0
  do k=2-ncb,nt+ncb-1
  do j=2-ncb,nt+ncb-1
  do i=2-ncb,nt+ncb-1
    nlast=cum(i-1,j,k,itx,ity,itz)
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np ! loop over particle
      ip=nlast+l
#ifndef finengp
      ! CIC
      tempx=4.*((/i,j,k/)-1)+4*(x(:,ip)+ishift+rshift)*x_resolution !-0.5
      idx1 = floor(tempx) + 1
      idx2 = idx1 + 1
      dx1 = idx1 - tempx
      dx2 = 1 - dx1
      idx1=idx1+nfb
      idx2=idx2+nfb ! positive integer indeces
      !idx1min=min(idx1,idx1min)
      !idx1max=max(idx1,idx1max)
      !idx2min=min(idx2,idx2min)
      !idx2max=max(idx2,idx2max)

      rho_f(idx1(1),idx1(2),idx1(3),ithread)=rho_f(idx1(1),idx1(2),idx1(3),ithread)+dx1(1)*dx1(2)*dx1(3)*mass_p
      rho_f(idx2(1),idx1(2),idx1(3),ithread)=rho_f(idx2(1),idx1(2),idx1(3),ithread)+dx2(1)*dx1(2)*dx1(3)*mass_p
      rho_f(idx1(1),idx2(2),idx1(3),ithread)=rho_f(idx1(1),idx2(2),idx1(3),ithread)+dx1(1)*dx2(2)*dx1(3)*mass_p
      rho_f(idx1(1),idx1(2),idx2(3),ithread)=rho_f(idx1(1),idx1(2),idx2(3),ithread)+dx1(1)*dx1(2)*dx2(3)*mass_p
      rho_f(idx1(1),idx2(2),idx2(3),ithread)=rho_f(idx1(1),idx2(2),idx2(3),ithread)+dx1(1)*dx2(2)*dx2(3)*mass_p
      rho_f(idx2(1),idx1(2),idx2(3),ithread)=rho_f(idx2(1),idx1(2),idx2(3),ithread)+dx2(1)*dx1(2)*dx2(3)*mass_p
      rho_f(idx2(1),idx2(2),idx1(3),ithread)=rho_f(idx2(1),idx2(2),idx1(3),ithread)+dx2(1)*dx2(2)*dx1(3)*mass_p
      rho_f(idx2(1),idx2(2),idx2(3),ithread)=rho_f(idx2(1),idx2(2),idx2(3),ithread)+dx2(1)*dx2(2)*dx2(3)*mass_p
#else
      ! NGP
      idxf=ceiling( (4.*((/i,j,k/)-1) + 4*((x(:,ip)+int(-128,1))+128.5)/256) )
      idxf=idxf+nfb ! positive integer indeces
      rho_f(idxf(1),idxf(2),idxf(3),ithread)=rho_f(idxf(1),idxf(2),idxf(3),ithread)+mass_p
#endif
    enddo
  enddo
  enddo
  enddo

  ! fine force ---------------------------------------------------------
  call sfftw_execute(plan_fft_fine)
!  print*, 'sum of Fourier-space rho_f =', sum(rho_f*1d0), maxval(rho_f)
  crho_f(:,:,:,ithread)=rho_f(:,:,:,ithread) ! back up
  do i_dim=1,3
    rho_f(::2,:,:,ithread)=-crho_f(2::2,:,:,ithread)*kern_f(i_dim,:,:,:)
    rho_f(2::2,:,:,ithread)=crho_f(::2,:,:,ithread)*kern_f(i_dim,:,:,:)
    call sfftw_execute(plan_ifft_fine)
    rho_f=rho_f/real(nfe)**3
    force_f(i_dim,:,:,:,ithread)=rho_f(nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1,ithread)
  enddo

  !print*, 'v_i2r_new'
  !print*, v_i2r_new
  !print*, vmax_new
  !print*, 'force2'
  !print*, maxval(abs(force_f(1,:,:,:,1))),maxval(abs(force_f(2,:,:,:,1))),maxval(abs(force_f(3,:,:,:,1)))

  ! fine velocity ------------------------------------------------
  do k=1,nt
  do j=1,nt
  do i=1,nt ! loop over coarse cell
    nlast=cum(i-1,j,k,itx,ity,itz)
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np ! loop over particle
      ip=nlast+l
#ifndef finengp
      tempx=4.*((/i,j,k/)-1)+4*(x(:,ip)+ishift+rshift)*x_resolution !-0.5
      idx1 = floor(tempx) + 1
      idx2 = idx1 + 1
      dx1 = idx1 - tempx
      dx2 = 1 - dx1
      idx1=idx1+nfb
      idx2=idx2+nfb ! positive integer indeces
      vtemp=v(:,ip)*v_i2r

      !print*, force_f(:,idx1(1),idx1(2),idx1(3),ithread)
      !print*, force_f(:,idx2(1),idx1(2),idx1(3),ithread)
      !print*, force_f(:,idx1(1),idx2(2),idx1(3),ithread)
      !print*, force_f(:,idx1(1),idx1(2),idx2(3),ithread)
      !print*, force_f(:,idx1(1),idx2(2),idx2(3),ithread)
      !print*, force_f(:,idx2(1),idx1(2),idx2(3),ithread)
      !print*, force_f(:,idx2(1),idx2(2),idx1(3),ithread)
      !print*, force_f(:,idx2(1),idx2(2),idx2(3),ithread)
      !print*, 'vtemp', vtemp
      vtemp=vtemp+force_f(:,idx1(1),idx1(2),idx1(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx2(1),idx1(2),idx1(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx1(1),idx2(2),idx1(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx1(1),idx1(2),idx2(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_f(:,idx1(1),idx2(2),idx2(3),ithread)*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
      vtemp=vtemp+force_f(:,idx2(1),idx1(2),idx2(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_f(:,idx2(1),idx2(2),idx1(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_f(:,idx2(1),idx2(2),idx2(3),ithread)*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)
      !print*, force_f(:,idx1(1),idx1(2),idx1(3),ithread)*a_mid*dt/6/pi
      !print*, 'vtemp', vtemp
#else
      ! NGP
      idxf=ceiling( (4.*((/i,j,k/)-1) + 4*((x(:,ip)+int(-128,1))+128.5)/256) )
      idxf=idxf+nfb ! positive integer indeces
      vtemp=v(:,ip)*v_i2r
      !print*, 'particle', ip
      !print*, 'v', vtemp
      vtemp=vtemp+force_f(:,idxf(1),idxf(2),idxf(3),ithread)*a_mid*dt/6/pi
      !print*, 'v', vtemp
#endif
      v(:,ip)=nint(vtemp/v_i2r_new,kind=izipv)
    enddo
  enddo
  enddo
  enddo

enddo
enddo
enddo

! update v_i2r by v_i2r_new
if (head) print*, '    update v_i2r by v_i2r_new'
v_i2r=v_i2r_new
vmax=vmax_new

!print*, 'fine mesh done'
!print*,'vmax', vmax

endif

sync all
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! call coarse_mesh
 ! call coarse_mass

if (coarse_force) then
if (head) print*, '  pm coarse'

! coarse_cic_mass ------------------------------------------
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
#ifndef coarsengp
      ! CIC
      !tempx(:) = ((/i,j,k/)-1) + (x(:,ip)+int(-128,1)+128.5)/256 - 0.5
      tempx=((/i,j,k/)-1)+(x(:,ip)+ishift+rshift)*x_resolution-0.5
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
#else
      ! NGP
      r3t(i,j,k)=r3t(i,j,k)+mass_p
#endif
    enddo
  enddo
  enddo
  enddo
  r3((itx-1)*nt+1:itx*nt,(ity-1)*nt+1:ity*nt,(itz-1)*nt+1:itz*nt)=r3t(1:nt,1:nt,1:nt)
  ! put center part of r3t into subset of r3
enddo
enddo
enddo
sync all
!print*, 'sum of r3', sum(r3*1.d0), mass_p, ip

!#ifdef debug_force
!r3=0
!r3(12,12,12)=1000
!#endif

! coarse force ----------------------------------------
if (head) print*, '    coarse cic force'
call pencil_fft_forward

!call fft_cube2pencil ! r3 -> cz
!call trans_zxy2xyz ! cz -> cx

! save complex rho_c into crho_c
crho_c(::2,:,:)=real(cxyz)
crho_c(2::2,:,:)=imag(cxyz)

do i_dim=1,3
  cxyz=cmplx(-crho_c(2::2,:,:)*kern_c(i_dim,:,:,:),crho_c(::2,:,:)*kern_c(i_dim,:,:,:))
  !cxyz=cmplx(-imag(cxyz)*kern_c(i_dim,:,:,:), real(cxyz)*kern_c(i_dim,:,:,:))
  call pencil_fft_backward
  !cxyz=crho*(0,1)*kern_c(i_dim,:,:,:)
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
!print*, 'f2_max_coarse = ',f2_max_coarse
!print*, 'min', minval(sum(force_c**2,1))

! coarse max velocity ---------------------------------------------
if (head) print*, '    coarse cic velocity'
do itz=1,nnt
do ity=1,nnt
do itx=1,nnt ! loop over tiles
  do k=1,nt
  do j=1,nt
  do i=1,nt
    nlast=cum(i-1,j,k,itx,ity,itz)
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np ! loop over particle
      ip=nlast+l
#ifndef coarsengp
      ! CIC
      !tempx(:) = ((/itx,ity,itz/)-1)*nt + ((/i,j,k/)-1) + (x(:,ip)+int(-128,1)+128.5)/256 - 0.5
      tempx=((/itx,ity,itz/)-1)*nt+((/i,j,k/)-1)+(x(:,ip)+ishift+rshift)*x_resolution-0.5
      idx1(:)=floor(tempx(:))+1
      idx2(:)=idx1(:)+1
      dx1(:)=idx1(:)-tempx(:)
      dx2(:)=1-dx1(:)
      vtemp=v(:,ip)*v_i2r
      vtemp=vtemp+force_c(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx2(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx1(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx1(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_c(:,idx1(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
      vtemp=vtemp+force_c(:,idx2(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_c(:,idx2(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx2(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)
#else
      ! NGP
      idxf=((/itx,ity,itz/)-1)*nt+(/i,j,k/)
      vtemp=v(:,ip)*v_i2r
      vtemp=vtemp+force_c(:,idxf(1),idxf(2),idxf(3))*a_mid*dt/6/pi
#endif

      vmax_new(1)=max(vmax_new(1),abs(vtemp(1)))
      vmax_new(2)=max(vmax_new(2),abs(vtemp(2)))
      vmax_new(3)=max(vmax_new(3),abs(vtemp(3)))

    enddo
  enddo
  enddo
  enddo
enddo
enddo
enddo
sync all

! sync vmax over images
if (head) print*, '    co_max vmax_new over images'
do i=1,nn**3
  vmax_new=max(vmax_new,vmax_new(:)[i])
enddo
v_i2r_new=vmax_new*v_resolution
if (head) print*, '    vmax_new =', vmax_new
sync all

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
#ifndef coarsengp
      ! CIC
      tempx=((/itx,ity,itz/)-1)*nt+((/i,j,k/)-1)+(x(:,ip)+ishift+rshift)*x_resolution-0.5
      idx1(:)=floor(tempx(:))+1
      idx2(:)=idx1(:)+1
      dx1(:)=idx1(:)-tempx(:)
      dx2(:)=1-dx1(:)
      vtemp=v(:,ip)*v_i2r
      vtemp=vtemp+force_c(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx2(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx1(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx1(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_c(:,idx1(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
      vtemp=vtemp+force_c(:,idx2(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
      vtemp=vtemp+force_c(:,idx2(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
      vtemp=vtemp+force_c(:,idx2(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)
#else
      ! NGP
      idxf=((/itx,ity,itz/)-1)*nt+(/i,j,k/)
      vtemp=v(:,ip)*v_i2r
      !print*, 'particle', ip
      !print*, 'v', vtemp
      vtemp=vtemp+force_c(:,idxf(1),idxf(2),idxf(3))*a_mid*dt/6/pi
      !print*, 'v',vtemp
#endif

      v(:,ip)=nint(vtemp/v_i2r_new,kind=izipv)

    enddo
  enddo
  enddo
  enddo
enddo
enddo
enddo
sync all

! update v_i2r by v_i2r_new
if (head) print*, '    update v_i2r by v_i2r_new'
v_i2r=v_i2r_new
vmax=vmax_new

endif

f2_max_fine(1,1,1)=maxval(f2_max_fine)
f2_max_pp(1,1,1)=maxval(f2_max_pp)
vmax(1)=maxval(vmax) ! max over x,y,z
sync all

! gather max forces and velocities
if (head) then
  do i=1,nn**3
    dt_fine(i)=f2_max_fine(1,1,1)[i]
    dt_pp(i)=f2_max_pp(1,1,1)[i]
    dt_coarse(i)=f2_max_coarse[i]
    dt_vmax(i)=vmax(1)[i]
  enddo
  dt_fine(1)=sqrt( 1.0 / (sqrt(maxval(dt_fine))*a_mid*GG) )
  dt_coarse(1)=sqrt( real(ncell) / (sqrt(maxval(dt_coarse))*a_mid*GG) )
  !dt_pp(1)=sqrt( real(ncell) / (sqrt(maxval(dt_coarse))*a_mid*GG) )
  dt_pp(1)=1000
  dt_vmax=vbuf*20/maxval(dt_vmax)
endif

if (head) then 
  !print*, 'particle mesh done'
  !print*, ''
endif
sync all

endsubroutine



endmodule
