#define ZA
#define mkdir
#define PID
!#define FIXED_IC

program initial_conditions
use penfft_fine
use iso_fortran_env , only : int64
implicit none
save

! nc: coarse grid per node per dim
! nf: fine grid per node per dim
real,parameter :: a=1/(1+z_i)
real,parameter :: Vphys2sim=1.0/(300.*sqrt(omega_m)*box/a/2/nf)
integer,parameter :: nk=1000
integer i,j,k,l,seedsize
real kmax,temp_r,temp_theta,pow,phi8
real(8) v8
integer(int64) :: t

integer nplocal
! power spectrum arrays
real, dimension(7,nk) :: tf    !CAMB
real, dimension(2,nc) :: pkm,pkn

integer,allocatable :: iseed(:)
real,allocatable :: rseed_all(:,:)

integer ind,dx,dxy,kg,mg,jg,ig,i0,j0,k0,itx,ity,itz,idx,imove,nf_shake
integer ileft,iright,nlen,nlast,g(3)
real kr,kx,ky,kz
real xi(10,nbin)

complex cx_temp(nf*nn/2+1,nf,nfpen)
real phi(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]
#ifndef ZA
  real phi1(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]
  real,dimension(nf,nf,nf) :: phixx,phiyy,phizz,phixy,phiyz,phizx
#endif
logical,parameter :: correct_kernel=.true.


! zip arrays
integer,parameter :: npt=nt*np_nc ! np / tile / dim
integer,parameter :: npb=ncb*np_nc
integer,parameter :: npmax=2*(npt+2*npb)**3
!integer rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
integer rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
integer rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
integer cume(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
integer(izipx) x(3,npmax)
integer(izipv) v(3,npmax)
#ifdef PID
  integer(2) pid(4,npmax)
#endif
real grad_max(3),vmax(3),v_i2r(3),vf

character (10) :: img_s, z_s
character (200) :: fn0,fn1,fn2,fn3,fn4

!equivalence(phixx,phixy)
!equivalence(phiyy,phiyz)
!equivalence(phizz,phizx)

if (head) then
  print*, 'Initial conditions on resolution', nf
  print*, 'Number of particles per side', np_nc*nc
  print*, 'np_2n3 =',np_2n3
endif


#ifdef mkdir
  call system('mkdir -p '//'.'//opath//'node'//image2str(this_image()-1))
#endif
fn0='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip0_'//image2str(this_image()-1)//'.dat'
fn1='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip1_'//image2str(this_image()-1)//'.dat'
fn2='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip2_'//image2str(this_image()-1)//'.dat'
fn3='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip3_'//image2str(this_image()-1)//'.dat'
#ifdef PID
  fn4='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zipid_'//image2str(this_image()-1)//'.dat'
#endif

sim%nplocal=1 ! will be overwritten
sim%a=1./(1+z_i)
sim%t=0
sim%tau=0

sim%nts=0

sim%dt_f_acc=1000
sim%dt_pp_acc=1000
sim%dt_c_acc=1000

sim%cur_checkpoint=0
sim%cur_proj=0
sim%cur_halo=0

sim%mass_p=real((nf)**3)/sim%nplocal ! will be overwritten
sim%v_r2i=1 ! will be overwritten
sim%shake_offset=0

sim%box=box
sim%rank=this_image()-1
sim%nn=int(nn,2)
sim%nnt=int(nnt,2)
sim%nt=int(nt,2)
sim%ncell=int(ncell,2)
sim%ncb=int(ncb,2)
sim%izipx=int(izipx,1)
sim%izipv=int(izipv,1)

sim%h0=h0
sim%omega_m=omega_m
sim%omega_l=omega_l
sim%s8=s8

sim%m_neu(1:3)=0
sim%vsim2phys=1.0/(300.*sqrt(omega_m)*box/a/2./ nf/nn)
sim%z_i=z_i


! initialize variables
if (np_2n3) then
  nf_shake=(ncell/np_nc)/2
else
  nf_shake=0
endif

! initvar
phi=0
tf=0
pkm=0
pkn=0


print*,'Creating FFT plans'
call create_penfft_fine_plan

! transferfnc
open(11,file='../tf/ith2_mnu0p05_z5_tk.dat',form='formatted')
read(11,*) tf
close(11)
!Delta^2
tf(2,:)=tf(2,:)**2 * tf(1,:)**(3+n_s) / (2*pi**2)
tf(3,:)=tf(3,:)**2 * tf(1,:)**(3+n_s) / (2*pi**2)
tf(6,:)=tf(6,:)**2 * tf(1,:)**(3+n_s) / (2*pi**2)
!dk
tf(4,1)=tf(1,2)/2
do k=2,nk-1
  tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
enddo
tf(4,nk)=tf(1,nk)-tf(1,nk-1)

v8=0
kmax=2*pi*sqrt(3.)*nf/2/box
do k=1,nk
  if (tf(1,k)>kmax) exit
  v8=v8+tf(2,k)*tophat(tf(1,k)*8)**2*tf(4,k)/tf(1,k)
enddo
tf(2:3,:)=tf(2:3,:)*(s8**2/v8)*Dgrow(a)**2

! noisemap
print*,'Generating random noise'

call random_seed()
call random_seed(size=seedsize)
allocate(iseed(seedsize))
allocate(rseed_all(seedsize,nn**3))
call system_clock(t)
  do i = 1, seedsize
      iseed(i) = lcg(t)
      !print*,t,seed(i)
  end do

call random_seed(put=iseed) ! generate seed using system time

!call random_seed(get=iseed)
print*, iseed

!call random_seed()
call random_number(cube)

!#ifdef IFYOUWANTBUG
  deallocate(iseed)
  deallocate(rseed_all)
!#endif

!! Box-Muller transform
print*,'Box-Muller transform'
do k=1,nf
do j=1,nf
do i=1,nf,2
  temp_theta=2*pi*cube(i,j,k)
  temp_r=sqrt(-2*log(1-cube(i+1,j,k)))
  cube(i,j,k)=temp_r*cos(temp_theta)
  cube(i+1,j,k)=temp_r*sin(temp_theta)
enddo
enddo
enddo

#ifdef FIXED_IC
  !!!!! test only
  open(11,file='initnoise.dat',access='stream')
  read(11) cube
  close(11)
  print*, 'READ IN NOISE MAP:', cube(1,1,1), cube(nf,nf,nf)
  !!!!! test only
#endif

!call cross_power(xi,cube,cube)
!open(11,file='initpower.dat',access='stream')
!write(11) xi
!close(11)

print*, 'Start ftran'
call fft_cube2pencil_fine
print*, 'Start transpose'
call trans_zxy2xyz_fine


! delta field
print*, 'Wiener filter noise'
do k=1,nfpen
do j=1,nf
do i=1,nf*nn/2+1
  ! global grid in Fourier space for i,j,k
  kg=(nn*(icz-1)+icy-1)*nfpen+k
  jg=(icx-1)*nf+j
  ig=i
  kz=mod(kg+nf/2-1,nf)-nf/2
  ky=mod(jg+nf/2-1,nf)-nf/2
  kx=ig-1
  kr=sqrt(kx**2+ky**2+kz**2)
  kr=max(kr,1.0)
  pow=interp_tf(2*pi*kr/box,1,2)/(4*pi*kr**3)
  cx(i,j,k)=cx(i,j,k)*sqrt(pow*nf**3)
enddo
enddo
enddo
if (head) cx(1,1,1)=0 ! DC frequency

sync all

cx_temp=cx ! backup delta

print*,'Start transpose'
call trans_xyz2zxy_fine
print*,'Start btran'
call ifft_pencil2cube_fine

! write initial overdensity
print*,'Write delta_L into file'
open(11,file='delta_L.dat',status='replace',access='stream')
write(11) cube/Dgrow(a)
close(11)
! potential field

do k=1,nfpen
do j=1,nf
do i=1,nf*nn+1,2
  kg=(nn*(icz-1)+icy-1)*nfpen+k
  jg=(icx-1)*nf+j
  ig=i
  kz=mod(kg+nf/2-1,nf)-nf/2
  ky=mod(jg+nf/2-1,nf)-nf/2
  kx=(ig-1)/2
  kz=2*sin(pi*kz/nf)
  ky=2*sin(pi*ky/nf)
  kx=2*sin(pi*kx/nf)
  kr=kx**2+ky**2+kz**2
  kr=max(kr,1.0/nf**2)
  cx((i+1)/2,j,k)=-4*pi/kr
enddo
enddo
enddo
if (head) cx(1,1,1)=0 ! DC frequency

!print*, cx(1:2,1,1)
!print*, cx_temp(1:2,1,1) 
!stop

sync all

if (correct_kernel) then

  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine
  
  open(11,file='laplace.dat',status='replace',access='stream')
  write(11) cube
  close(11)

  sync all

  phi8=0.0

  phi8=cube(9,1,1)[image1d(1,1,1)]+cube(1,9,1)[image1d(1,1,1)]+cube(1,1,9)[image1d(1,1,1)]
  phi8=phi8+cube(nf-7,1,1)[image1d(nn,1,1)]
  phi8=phi8+cube(1,nf-7,1)[image1d(1,nn,1)]
  phi8=phi8+cube(1,1,nf-7)[image1d(1,1,nn)]
  phi8=phi8/6

!print*, 'phi8='
!print*, cube(9,1,1)[image1d(1,1,1)]
!print*, cube(1,9,1)[image1d(1,1,1)]
!print*, cube(1,1,9)[image1d(1,1,1)]
!print*, cube(nf-7,1,1)[image1d(nn,1,1)]
!print*, cube(1,nf-7,1)[image1d(1,nn,1)]
!print*, cube(1,1,nf-7)[image1d(1,1,nn)]

  do k=1,nf
  do j=1,nf
  do i=1,nf
    kg=k+nf*(icx-1)
    jg=j+nf*(icy-1)
    ig=i+nf*(icz-1)
    kx=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kz=mod(ig+nf/2-1,nf)-nf/2
    kr=sqrt(kx**2+ky**2+kz**2)
    if (kr>8) then
      cube(i,j,k)=cube(i,j,k)-(phi8+1/8.)
    elseif (kr>0) then
      cube(i,j,k)=-1/kr
    elseif (kr==0) then
      cube(i,j,k)=-2.5
    endif
  enddo
  enddo
  enddo
  sync all

  call fft_cube2pencil_fine
  call trans_zxy2xyz_fine
endif

! Complex multiply density field with potential kernel
cx=real(cx)*cx_temp
cx_temp=cx  ! backup phi(k)

!!!!!! output phi in real space
call trans_xyz2zxy_fine
call ifft_pencil2cube_fine
phi=0
phi(1:nf,1:nf,1:nf)=cube ! phi1
open(11,file='phi1.dat',status='replace',access='stream')
write(11) cube
close(11)
call fft_cube2pencil_fine
call trans_zxy2xyz_fine
!!!!!!

#ifndef ZA
  ! 2LPT

  ! phi_xx
  do k=1,nfpen
  do j=1,nf
  do i=1,nf*nn+1,2
    kg=(nn*(icz-1)+icy-1)*nfpen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=(ig-1)/2
    kz=2*sin(pi*kz/nf) ! ?
    ky=2*sin(pi*ky/nf) ! ?
    kx=2*sin(pi*kx/nf) ! ?
    kr=kx**2+ky**2+kz**2
    cx((i+1)/2,j,k)=(-kx**2)/(4*pi)*cx_temp((i+1)/2,j,k)
  enddo
  enddo
  enddo
  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine
  phixx=cube

  ! phi_yy
  do k=1,nfpen
  do j=1,nf
  do i=1,nf*nn+1,2
    kg=(nn*(icz-1)+icy-1)*nfpen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=(ig-1)/2
    kz=2*sin(pi*kz/nf) ! ?
    ky=2*sin(pi*ky/nf) ! ?
    kx=2*sin(pi*kx/nf) ! ?
    kr=kx**2+ky**2+kz**2
    cx((i+1)/2,j,k)=(-ky**2)/(4*pi)*cx_temp((i+1)/2,j,k)
  enddo
  enddo
  enddo
  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine
  phiyy=cube

  ! phi_zz
  do k=1,nfpen
  do j=1,nf
  do i=1,nf*nn+1,2
    kg=(nn*(icz-1)+icy-1)*nfpen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=(ig-1)/2
    kz=2*sin(pi*kz/nf) ! ?
    ky=2*sin(pi*ky/nf) ! ?
    kx=2*sin(pi*kx/nf) ! ?
    kr=kx**2+ky**2+kz**2
    cx((i+1)/2,j,k)=(-kz**2)/(4*pi)*cx_temp((i+1)/2,j,k)
  enddo
  enddo
  enddo
  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine
  phizz=cube

  ! phi_xy
  do k=1,nfpen
  do j=1,nf
  do i=1,nf*nn+1,2
    kg=(nn*(icz-1)+icy-1)*nfpen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=(ig-1)/2
    kz=2*sin(pi*kz/nf) ! ?
    ky=2*sin(pi*ky/nf) ! ?
    kx=2*sin(pi*kx/nf) ! ?
    kr=kx**2+ky**2+kz**2
    cx((i+1)/2,j,k)=(-kx*ky)/(4*pi)*cx_temp((i+1)/2,j,k)
  enddo
  enddo
  enddo
  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine
  phixy=cube

  ! phi_yz
  do k=1,nfpen
  do j=1,nf
  do i=1,nf*nn+1,2
    kg=(nn*(icz-1)+icy-1)*nfpen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=(ig-1)/2
    kz=2*sin(pi*kz/nf) ! ?
    ky=2*sin(pi*ky/nf) ! ?
    kx=2*sin(pi*kx/nf) ! ?
    kr=kx**2+ky**2+kz**2
    cx((i+1)/2,j,k)=(-ky*kz)/(4*pi)*cx_temp((i+1)/2,j,k)
  enddo
  enddo
  enddo
  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine
  phiyz=cube

  ! phi_zx
  do k=1,nfpen
  do j=1,nf
  do i=1,nf*nn+1,2
    kg=(nn*(icz-1)+icy-1)*nfpen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=(ig-1)/2
    kz=2*sin(pi*kz/nf) ! ?
    ky=2*sin(pi*ky/nf) ! ?
    kx=2*sin(pi*kx/nf) ! ?
    kr=kx**2+ky**2+kz**2
    cx((i+1)/2,j,k)=(-kz*kx)/(4*pi)*cx_temp((i+1)/2,j,k)
  enddo
  enddo
  enddo
  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine
  phizx=cube

  ! phi2_source in real space
  cube=phixx*phiyy+phiyy*phizz+phizz*phixx-phixy**2-phiyz**2-phizx**2

  !! diff in real space
  !phixx=phi(0:nf-1,1:nf,1:nf)-2*phi(1:nf,1:nf,1:nf)+phi(2:nf+1,1:nf,1:nf)
  !phiyy=phi(1:nf,0:nf-1,1:nf)-2*phi(1:nf,1:nf,1:nf)+phi(1:nf,2:nf+1,1:nf)
  !phizz=phi(1:nf,1:nf,0:nf-1)-2*phi(1:nf,1:nf,1:nf)+phi(1:nf,1:nf,2:nf+1)

  !cube=phixx*phiyy+phiyy*phizz+phizz*phixx

  !phixy=(phi(2:nf+1,2:nf+1,1:nf)-phi(0:nf-1,2:nf+1,1:nf)+phi(0:nf-1,0:nf-1,1:nf)-phi(2:nf+1,0:nf-1,1:nf))/4
  !phiyz=(phi(1:nf,2:nf+1,2:nf+1)-phi(1:nf,0:nf-1,2:nf+1)+phi(1:nf,0:nf-1,0:nf-1)-phi(1:nf,2:nf+1,0:nf-1))/4
  !phizx=(phi(2:nf+1,1:nf,2:nf+1)-phi(0:nf-1,1:nf,2:nf+1)+phi(0:nf-1,1:nf,0:nf-1)-phi(2:nf+1,1:nf,0:nf-1))/4

  !cube=cube-phixy**2-phiyz**2-phizx**2 ! source term of phi2

  open(11,file='phi2_source.dat',status='replace',access='stream')
  write(11) cube
  close(11)

  sync all

  call fft_cube2pencil_fine
  call trans_zxy2xyz_fine

  ! solve phi2
  do k=1,nfpen
  do j=1,nf
  do i=1,nf*nn+1,2
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*nfpen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=(ig-1)/2
    kz=2*sin(pi*kz/nf) ! ?
    ky=2*sin(pi*ky/nf) ! ?
    kx=2*sin(pi*kx/nf) ! ?
    kr=kx**2+ky**2+kz**2
    kr=max(kr,1.0/nf**2)
    cx((i+1)/2,j,k)=(-4*pi/kr)*cx((i+1)/2,j,k)
  enddo
  enddo
  enddo
  if (head) cx(1,1,1)=0 ! DC frequency

  call trans_xyz2zxy_fine
  call ifft_pencil2cube_fine

  open(11,file='phi2.dat',status='replace',access='stream')
  write(11) cube
  close(11)

  phixx=phi(1:nf,1:nf,1:nf) ! backup phi1
  phiyy=cube ! backup phi2

  ! phi for delta_x
  print*,'Dgrow(a)', Dgrow(a)
  print*, a, Dgrow(a)
  print*, 'a=1',Dgrow(1.)

  phi(1:nf,1:nf,1:nf)=phixx-Dgrow(a)*phiyy ! corrected phi for positions

  open(11,file='phi12.dat',status='replace',access='stream')
  write(11) phi(1:nf,1:nf,1:nf)
  close(11)
#endif


!#ifdef ZA
!  phi(1:nf,1:nf,1:nf)=phixx
!#endif

! buffer phi

phi(:0,:,:)=phi(nf-nfb:nf,:,:)[image1d(inx,icy,icz)]
phi(nf+1:,:,:)=phi(1:nfb+1,:,:)[image1d(ipx,icy,icz)]
phi(:,:0,:)=phi(:,nf-nfb:nf,:)[image1d(icx,iny,icz)]
phi(:,nf+1:,:)=phi(:,1:nfb+1,:)[image1d(icx,ipy,icz)]
phi(:,:,:0)=phi(:,:,nf-nfb:nf)[image1d(icx,icy,inz)]
phi(:,:,nf+1:)=phi(:,:,1:nfb+1)[image1d(icx,icy,ipz)]

#ifndef ZA
  ! buffer phi1
  phi1(1:nf,1:nf,1:nf)=phixx
  phi1(:0,:,:)=phi1(nf-nfb:nf,:,:)[image1d(inx,icy,icz)]
  phi1(nf+1:,:,:)=phi1(1:nfb+1,:,:)[image1d(ipx,icy,icz)]
  phi1(:,:0,:)=phi1(:,nf-nfb:nf,:)[image1d(icx,iny,icz)]
  phi1(:,nf+1:,:)=phi1(:,1:nfb+1,:)[image1d(icx,ipy,icz)]
  phi1(:,:,:0)=phi1(:,:,nf-nfb:nf)[image1d(icx,icy,inz)]
  phi1(:,:,nf+1:)=phi1(:,:,1:nfb+1)[image1d(icx,icy,ipz)]
#endif
call destroy_penfft_fine_plan

vf=vfactor(a)
print*, 'vf',vf
grad_max(1)=maxval(abs(phi(0:nf-1,1:nf,1:nf)-phi(2:nf+1,1:nf,1:nf)))
grad_max(2)=maxval(abs(phi(1:nf,0:nf-1,1:nf)-phi(1:nf,2:nf+1,1:nf)))
grad_max(3)=maxval(abs(phi(1:nf,1:nf,0:nf-1)-phi(1:nf,1:nf,2:nf+1)))

vmax=grad_max/2/(4*pi)*vf

v_i2r=vmax*v_resolution

print*, 'grad_max',grad_max
print*, 'vmax',vmax
! create particles (no communication)

print*, 'phi', phi(1:4,1,1)

open(10,file=fn0,status='replace',access='stream')
open(11,file=fn1,status='replace',access='stream')
open(12,file=fn2,status='replace',access='stream')
#ifdef PID
open(14,file=fn4,status='replace',access='stream')
#endif

write(12) sim ! occupying 128 bytes

nplocal=0
do itz=1,nnt
do ity=1,nnt
do itx=1,nnt
  iright=0
  rhoce=0
  rholocal=0
  do k=1-npb,npt+npb ! calculate coarse mesh density
  do j=1-npb,npt+npb
  do i=1-npb,npt+npb
  do imove=0,nf_shake
    k0=nft*(itz-1)+(ncell/np_nc)*(k-1)+1+imove
    j0=nft*(ity-1)+(ncell/np_nc)*(j-1)+1+imove
    i0=nft*(itx-1)+(ncell/np_nc)*(i-1)+1+imove
    g(1)=ceiling( real(i-1)/np_nc+(0.5+imove)/ncell + (phi(i0-1,j0,k0)-phi(i0+1,j0,k0))/2/(4*pi)/ncell )
    g(2)=ceiling( real(j-1)/np_nc+(0.5+imove)/ncell + (phi(i0,j0-1,k0)-phi(i0,j0+1,k0))/2/(4*pi)/ncell )
    g(3)=ceiling( real(k-1)/np_nc+(0.5+imove)/ncell + (phi(i0,j0,k0-1)-phi(i0,j0,k0+1))/2/(4*pi)/ncell )
    rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1
    !if (i==1 .and. j==1 .and. k==1) then
    !  print*, 4*( real(i-1)/np_nc+0.5/ncell + (phi(i0-1,j0,k0)-phi(i0+1,j0,k0))/2/(4*pi)/ncell )
    !  print*, 4*( real(j-1)/np_nc+0.5/ncell + (phi(i0,j0-1,k0)-phi(i0,j0+1,k0))/2/(4*pi)/ncell )
    !  print*, 4*( real(k-1)/np_nc+0.5/ncell + (phi(i0,j0,k0-1)-phi(i0,j0,k0+1))/2/(4*pi)/ncell )
    !endif
  enddo
  enddo
  enddo
  enddo

  cume=cumsum3(rhoce)

  do k=1-npb,npt+npb ! create particles in extended mesh
  do j=1-npb,npt+npb
  do i=1-npb,npt+npb
  do imove=0,nf_shake
    k0=nft*(itz-1)+(ncell/np_nc)*(k-1)+1+imove
    j0=nft*(ity-1)+(ncell/np_nc)*(j-1)+1+imove
    i0=nft*(itx-1)+(ncell/np_nc)*(i-1)+1+imove
    g(1)=ceiling( real(i-1)/np_nc+(0.5+imove)/ncell + (phi(i0-1,j0,k0)-phi(i0+1,j0,k0))/2/(4*pi)/ncell )
    g(2)=ceiling( real(j-1)/np_nc+(0.5+imove)/ncell + (phi(i0,j0-1,k0)-phi(i0,j0+1,k0))/2/(4*pi)/ncell )
    g(3)=ceiling( real(k-1)/np_nc+(0.5+imove)/ncell + (phi(i0,j0,k0-1)-phi(i0,j0,k0+1))/2/(4*pi)/ncell )
    rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
    idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3))
    x(1,idx)=floor(( real(i-1)/np_nc+(0.5+imove)/ncell + (phi(i0-1,j0,k0)-phi(i0+1,j0,k0))/2/(4*pi)/ncell )/x_resolution)
    x(2,idx)=floor(( real(j-1)/np_nc+(0.5+imove)/ncell + (phi(i0,j0-1,k0)-phi(i0,j0+1,k0))/2/(4*pi)/ncell )/x_resolution)
    x(3,idx)=floor(( real(k-1)/np_nc+(0.5+imove)/ncell + (phi(i0,j0,k0-1)-phi(i0,j0,k0+1))/2/(4*pi)/ncell )/x_resolution)
    v(1,idx)=nint( (phi(i0-1,j0,k0)-phi(i0+1,j0,k0))/2/(4*pi)*vf/v_i2r(1), kind=izipv)
    v(2,idx)=nint( (phi(i0,j0-1,k0)-phi(i0,j0+1,k0))/2/(4*pi)*vf/v_i2r(2), kind=izipv)
    v(3,idx)=nint( (phi(i0,j0,k0-1)-phi(i0,j0,k0+1))/2/(4*pi)*vf/v_i2r(3), kind=izipv)
#ifdef PID
    pid(1,idx)=this_image()-1
    pid(2:4,idx)=floor(( ((/itx,ity,itz/)-1)*nft+(ncell/np_nc)*((/i,j,k/)-1)+0.5+imove )/nf*2**16-2**15)
#endif
    !if (itx==1 .and. ity==1 .and. itz==1 .and. i==1 .and. j==1 .and. k==1) then
    !  print*, 'idx',idx
    !  !print*, x(:,idx)
    !  print*, (x(:,idx)+ishift+rshift)*x_resolution*ncell + (g-1)*ncell
    !  print*, v(:,idx)*v_i2r
    !endif
    !if (itx==nnt .and. ity==nnt .and. itz==nnt .and. i==npt .and. j==npt .and. k==npt) then
    !  print*, 'idx',idx
    !  !print*, x(:,idx)
    !  print*, (x(:,idx)+ishift+rshift)*x_resolution*ncell + (g-1)*ncell
    !  print*, v(:,idx)*v_i2r
    !endif
  enddo
  enddo
  enddo
  enddo

  do k=1,nt ! compact particle
  do j=1,nt
    ileft=iright+1
    nlast=cume(nt,j,k)
    nlen=nlast-cume(0,j,k)
    iright=ileft+nlen-1
    x(:,ileft:iright)=x(:,nlast-nlen+1:nlast)
    v(:,ileft:iright)=v(:,nlast-nlen+1:nlast)
#ifdef PID
    pid(:,ileft:iright)=pid(:,nlast-nlen+1:nlast)
#endif
  enddo
  enddo

  write(10) x(:,1:iright)
  write(11) v(:,1:iright)
#ifdef PID
  write(14) pid(:,1:iright)
#endif
  write(12) rhoce(1:nt,1:nt,1:nt)

  nplocal=nplocal+iright

enddo
enddo
enddo

close(10)
close(11)
#ifdef PID
close(14)
#endif

sim%nplocal=nplocal
sim%mass_p=real((nf)**3)/nplocal
sim%v_r2i=1/v_i2r
rewind(12)
write(12) sim
close(12)

call print_header(sim)

if (head) print*, 'initial condition done'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

function interp_tf(kr,ix,iy)
implicit none
integer ix,iy,ii,i1,i2
real kr,xx,yy,x1,x2,y1,y2,interp_tf
i1=1
i2=nk
do while (i2-i1>1)
  ii=(i1+i2)/2
  if (kr>tf(ix,ii)) then
    i1=ii
  else
    i2=ii
  endif
enddo
x1=log(tf(ix,i1))
y1=log(tf(iy,i1))
x2=log(tf(ix,i2))
y2=log(tf(iy,i2))
xx=log(kr)
yy=y1+(y2-y1)*(xx-x1)/(x2-x1)
interp_tf=exp(yy)
endfunction interp_tf


function tophat(x)
  implicit none
  real :: x,tophat
  if (x/=0) then
    tophat=3*(sin(x)-cos(x)*x)/x**3
  else
    tophat=1
  endif
endfunction tophat

function Dgrow(a)
  implicit none
  real, parameter :: om=omega_m
  real, parameter :: ol=omega_l
  real a
  real Dgrow
  real g,ga,hsq,oma,ola
  hsq=om/a**3+(1-om-ol)/a**2+ol
  oma=om/(a**3*hsq)
  ola=ol/hsq
  g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
  ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
  Dgrow=a*ga/g
end function Dgrow

function vfactor(a)
  implicit none
  real :: a
  real :: H,km,lm
  real :: vfactor
  lm=omega_l/omega_m
  km=(1-omega_m-omega_l)/omega_m
  H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
  vfactor=a**2*H
endfunction vfactor

Function lcg(s) !// Linear congruential generator
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
End Function lcg

end
