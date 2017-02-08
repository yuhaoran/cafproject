#define ZA
#define mkdir
#define PID
!#define READ_SEED
!#define READ_NOISE

program initial_conditions
  use pencil_fft
  !use powerspectrum
  use iso_fortran_env , only : int64
  implicit none
  save

  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  real,parameter :: a=1/(1+z_i)
  real,parameter :: Vphys2sim=1.0/(300.*sqrt(omega_m)*box/a/2/nf)
  integer,parameter :: nk=nk_tf   !1000
  integer i,j,k,seedsize
  real kmax,temp_r,temp_theta,pow,phi8,temp8[*]
  real(8) v8, norm
  integer(int64) :: time64

  integer nplocal
  ! power spectrum arrays
  real, dimension(7,nk) :: tf    !CAMB
  real, dimension(2,nc) :: pkm,pkn

  integer,allocatable :: iseed(:)
  real,allocatable :: rseed_all(:,:)

  integer ind,dx,dxy,kg,mg,jg,ig,ii,jj,kk,itx,ity,itz,idx,imove,nf_shake
  integer ileft,iright,nlen,nlast,g(3)
  real kr,kx,ky,kz
  !real xi(10,nbin)

  complex cx_temp(nf*nn/2+1,nf,npen)
  real phi(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]
#ifndef ZA
    real phi1(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]
    real,dimension(nf,nf,nf) :: phixx,phiyy,phizz,phixy,phiyz,phizx
#endif
  logical,parameter :: correct_kernel=.true.


  ! zip arrays
  integer,parameter :: npt=nt*np_nc ! np / tile / dim !64
  integer,parameter :: npb=ncb*np_nc !24
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

  !equivalence(phixx,phixy)
  !equivalence(phiyy,phiyz)
  !equivalence(phizz,phizx)

  call geometry

  if (head) then
    print*, 'Initial conditions on resolution', ng
    print*, 'Number of particles per side', np_nc*nc
    print*, 'np_2n3 =',np_2n3
    print*, 'output: ', opath
    print*, 'head image number',icx,icy,icz
  endif


#ifdef mkdir
    call system('mkdir -p '//opath//'node'//image2str(this_image()-1))
    !print*, 'mkdir -p '//opath//'node'//image2str(this_image()-1)
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
  call create_penfft_plan

  ! transferfnc
  open(11,file='../configs/mmh_transfer/simtransfer_bao.dat',form='formatted') ! for Xin
  read(11,*) tf
  close(11)

  ! normalization
  norm=2.*pi**2.*(h0/100.)**4*(h0/100./0.05)**(n_s-1)
  write(*,*) 'norm:', norm


  !Delta^2
  tf(2,:)=tf(2,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
  tf(3,:)=tf(3,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
  tf(6,:)=tf(6,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
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
  write(*,*) 's8**2/v8:', v8, s8**2/v8

  !!tf(2:3,:)=tf(2:3,:)*(s8**2/v8)*Dgrow(a)**2
  tf(2:3,:)= scalar_amp*tf(2:3,:)*Dgrow(a)**2

  ! noisemap
  print*,'Generating random noise'

  call random_seed(size=seedsize)
  print*,'min seedsize =', seedsize
  seedsize=max(seedsize,12)

  allocate(iseed(seedsize))
  allocate(rseed_all(seedsize,nn**3))

#ifdef READ_SEED
    ! Read seeds from ../configs

    call system('cp ../configs/seed_'//image2str(this_image()-1)//'.bin '//opath//'node'//image2str(this_image()-1)) ! for Xin

    ! need to use seed[image_number].dat for parallel

    open(11,file=output_dir()//'seed'//output_suffix(),status='old',access='stream')
    read(11) iseed
    close(11)
    ! Input iseed
    call random_seed(put=iseed)
    print*, 'iseed', iseed
#else
    ! Generate at least 12 seeds according to system clock
    call system_clock(time64)
    do i = 1, seedsize
      iseed(i) = lcg(time64)
      print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
    enddo
    ! Input iseed to system
    call random_seed(put=iseed)
    print*, 'iseed', iseed
    ! Write iseed into file
    open(11,file=output_dir()//'seed'//output_suffix(),status='replace',access='stream')
    write(11) iseed
    close(11)
#endif

  call random_number(r3)

#ifdef READ_NOISE
      !!!!! test only
      open(11,file='initnoise.bin',access='stream')
      read(11) r3
      close(11)
      print*, 'READ IN NOISE MAP:', r3(1,1,1), r3(ng,ng,ng)
      !!!!! test only
#endif

  deallocate(iseed)
  deallocate(rseed_all)

  ! to-do: rseed_all for image-dependent seed

  !! Box-Muller transform
  print*,'Box-Muller transform'
  do k=1,nf
  do j=1,nf
  do i=1,nf,2
    temp_theta=2*pi*r3(i,j,k)
    temp_r=sqrt(-2*log(1-r3(i+1,j,k)))
    r3(i,j,k)=temp_r*cos(temp_theta)
    r3(i+1,j,k)=temp_r*sin(temp_theta)
  enddo
  enddo
  enddo

  !call cross_power(xi,r3,r3)
  !open(11,file='initpower.dat',access='stream')
  !write(11) xi
  !close(11)

  print*, 'Start ftran'
  call pencil_fft_forward

  ! delta field
  print*, 'Wiener filter noise'
  do k=1,npen
  do j=1,nf
  do i=1,nf*nn/2+1
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nf/2-1,nf)-nf/2
    ky=mod(jg+nf/2-1,nf)-nf/2
    kx=ig-1
    kr=sqrt(kx**2+ky**2+kz**2)
    kr=max(kr,1.0)
    pow=interp_tf(2*pi*kr/box,1,2)/(4*pi*kr**3)
    cxyz(i,j,k)=cxyz(i,j,k)*sqrt(pow*nf**3)
  enddo
  enddo
  enddo
  if (head) cxyz(1,1,1)=0 ! DC frequency

  sync all

  cx_temp=cxyz ! backup delta

  print*,'Start btran'
  call pencil_fft_backward

  ! write initial overdensity
  print*,'Write delta_L into file'
  print*,'Growth factor Dgrow(a) =',Dgrow(a)
  open(11,file=output_dir()//'delta_L'//output_suffix(),status='replace',access='stream')
  write(11) r3/Dgrow(a)
  close(11)
  ! potential field

  do k=1,npen
  do j=1,nf
  do i=1,nf*nn+1,2
    kg=(nn*(icz-1)+icy-1)*npen+k
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
    cxyz((i+1)/2,j,k)=-4*pi/kr
  enddo
  enddo
  enddo
  if (head) cxyz(1,1,1)=0 ! DC frequency

  sync all

  if (correct_kernel) then

    call pencil_fft_backward

    !open(11,file='laplace.dat',status='replace',access='stream')
    !write(11) r3
    !close(11)

    sync all
    temp8=0
    if (rank==0) temp8=temp8+r3(9,1,1)+r3(1,9,1)+r3(1,1,9)
    sync all
    if (icx==nn .and. icy==1 .and. icz==1) temp8=temp8+r3(nf-7,1,1)
    sync all
    if (icx==1 .and. icy==nn .and. icz==1) temp8=temp8+r3(1,nf-7,1)
    sync all
    if (icx==1 .and. icy==1 .and. icz==nn) temp8=temp8+r3(1,1,nf-7)
    sync all

    do i=1,nn**3
      phi8=phi8+temp8[i]
    enddo
    sync all
    phi8=phi8/6
    !phi8=0.0

    !phi8=r3(9,1,1)[image1d(1,1,1)]+r3(1,9,1)[image1d(1,1,1)]+r3(1,1,9)[image1d(1,1,1)]
    !phi8=phi8+r3(nf-7,1,1)[image1d(nn,1,1)]
    !phi8=phi8+r3(1,nf-7,1)[image1d(1,nn,1)]
    !phi8=phi8+r3(1,1,nf-7)[image1d(1,1,nn)]
    !phi8=phi8/6

  !print*, 'phi8='
  !print*, r3(9,1,1)[image1d(1,1,1)]
  !print*, r3(1,9,1)[image1d(1,1,1)]
  !print*, r3(1,1,9)[image1d(1,1,1)]
  !print*, r3(nf-7,1,1)[image1d(nn,1,1)]
  !print*, r3(1,nf-7,1)[image1d(1,nn,1)]
  !print*, r3(1,1,nf-7)[image1d(1,1,nn)]
print*,'phi8 =',phi8
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
        r3(i,j,k)=r3(i,j,k)-(phi8+1/8.)
      elseif (kr>0) then
        r3(i,j,k)=-1/kr
      elseif (kr==0) then
        r3(i,j,k)=-2.5
      endif
    enddo
    enddo
    enddo
    sync all

    call pencil_fft_forward
  endif

  ! Complex multiply density field with potential kernel
  cxyz=real(cxyz)*cx_temp
  cx_temp=cxyz  ! backup phi(k)

  !!!!!! output phi in real space
  call pencil_fft_backward
  phi=0
  phi(1:nf,1:nf,1:nf)=r3 ! phi1
  !open(11,file='phi1.dat',status='replace',access='stream')
  !write(11) r3
  !close(11)
  call pencil_fft_forward
  !!!!!!

#ifndef ZA
    ! 2LPT

    ! phi_xx
    do k=1,npen
    do j=1,nf
    do i=1,nf*nn+1,2
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*nf+j
      ig=i
      kz=mod(kg+nf/2-1,nf)-nf/2
      ky=mod(jg+nf/2-1,nf)-nf/2
      kx=(ig-1)/2
      kz=2*sin(pi*kz/nf) ! ?
      ky=2*sin(pi*ky/nf) ! ?
      kx=2*sin(pi*kx/nf) ! ?
      kr=kx**2+ky**2+kz**2
      cxyz((i+1)/2,j,k)=(-kx**2)/(4*pi)*cx_temp((i+1)/2,j,k)
    enddo
    enddo
    enddo
    call pencil_fft_backward
    phixx=r3

    ! phi_yy
    do k=1,npen
    do j=1,nf
    do i=1,nf*nn+1,2
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*nf+j
      ig=i
      kz=mod(kg+nf/2-1,nf)-nf/2
      ky=mod(jg+nf/2-1,nf)-nf/2
      kx=(ig-1)/2
      kz=2*sin(pi*kz/nf) ! ?
      ky=2*sin(pi*ky/nf) ! ?
      kx=2*sin(pi*kx/nf) ! ?
      kr=kx**2+ky**2+kz**2
      cxyz((i+1)/2,j,k)=(-ky**2)/(4*pi)*cx_temp((i+1)/2,j,k)
    enddo
    enddo
    enddo
    call pencil_fft_backward
    phiyy=r3

    ! phi_zz
    do k=1,npen
    do j=1,nf
    do i=1,nf*nn+1,2
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*nf+j
      ig=i
      kz=mod(kg+nf/2-1,nf)-nf/2
      ky=mod(jg+nf/2-1,nf)-nf/2
      kx=(ig-1)/2
      kz=2*sin(pi*kz/nf) ! ?
      ky=2*sin(pi*ky/nf) ! ?
      kx=2*sin(pi*kx/nf) ! ?
      kr=kx**2+ky**2+kz**2
      cxyz((i+1)/2,j,k)=(-kz**2)/(4*pi)*cx_temp((i+1)/2,j,k)
    enddo
    enddo
    enddo
    call pencil_fft_backward
    phizz=r3

    ! phi_xy
    do k=1,npen
    do j=1,nf
    do i=1,nf*nn+1,2
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*nf+j
      ig=i
      kz=mod(kg+nf/2-1,nf)-nf/2
      ky=mod(jg+nf/2-1,nf)-nf/2
      kx=(ig-1)/2
      kz=2*sin(pi*kz/nf) ! ?
      ky=2*sin(pi*ky/nf) ! ?
      kx=2*sin(pi*kx/nf) ! ?
      kr=kx**2+ky**2+kz**2
      cxyz((i+1)/2,j,k)=(-kx*ky)/(4*pi)*cx_temp((i+1)/2,j,k)
    enddo
    enddo
    enddo
    call pencil_fft_backward
    phixy=r3

    ! phi_yz
    do k=1,npen
    do j=1,nf
    do i=1,nf*nn+1,2
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*nf+j
      ig=i
      kz=mod(kg+nf/2-1,nf)-nf/2
      ky=mod(jg+nf/2-1,nf)-nf/2
      kx=(ig-1)/2
      kz=2*sin(pi*kz/nf) ! ?
      ky=2*sin(pi*ky/nf) ! ?
      kx=2*sin(pi*kx/nf) ! ?
      kr=kx**2+ky**2+kz**2
      cxyz((i+1)/2,j,k)=(-ky*kz)/(4*pi)*cx_temp((i+1)/2,j,k)
    enddo
    enddo
    enddo
    call pencil_fft_backward
    phiyz=r3

    ! phi_zx
    do k=1,npen
    do j=1,nf
    do i=1,nf*nn+1,2
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*nf+j
      ig=i
      kz=mod(kg+nf/2-1,nf)-nf/2
      ky=mod(jg+nf/2-1,nf)-nf/2
      kx=(ig-1)/2
      kz=2*sin(pi*kz/nf) ! ?
      ky=2*sin(pi*ky/nf) ! ?
      kx=2*sin(pi*kx/nf) ! ?
      kr=kx**2+ky**2+kz**2
      cxyz((i+1)/2,j,k)=(-kz*kx)/(4*pi)*cx_temp((i+1)/2,j,k)
    enddo
    enddo
    enddo
    call pencil_fft_backward
    phizx=r3

    ! phi2_source in real space
    r3=phixx*phiyy+phiyy*phizz+phizz*phixx-phixy**2-phiyz**2-phizx**2

    !! diff in real space
    !phixx=phi(0:nf-1,1:nf,1:nf)-2*phi(1:nf,1:nf,1:nf)+phi(2:nf+1,1:nf,1:nf)
    !phiyy=phi(1:nf,0:nf-1,1:nf)-2*phi(1:nf,1:nf,1:nf)+phi(1:nf,2:nf+1,1:nf)
    !phizz=phi(1:nf,1:nf,0:nf-1)-2*phi(1:nf,1:nf,1:nf)+phi(1:nf,1:nf,2:nf+1)

    !r3=phixx*phiyy+phiyy*phizz+phizz*phixx

    !phixy=(phi(2:nf+1,2:nf+1,1:nf)-phi(0:nf-1,2:nf+1,1:nf)+phi(0:nf-1,0:nf-1,1:nf)-phi(2:nf+1,0:nf-1,1:nf))/4
    !phiyz=(phi(1:nf,2:nf+1,2:nf+1)-phi(1:nf,0:nf-1,2:nf+1)+phi(1:nf,0:nf-1,0:nf-1)-phi(1:nf,2:nf+1,0:nf-1))/4
    !phizx=(phi(2:nf+1,1:nf,2:nf+1)-phi(0:nf-1,1:nf,2:nf+1)+phi(0:nf-1,1:nf,0:nf-1)-phi(2:nf+1,1:nf,0:nf-1))/4

    !r3=r3-phixy**2-phiyz**2-phizx**2 ! source term of phi2

    open(11,file='phi2_source.dat',status='replace',access='stream')
    write(11) r3
    close(11)

    sync all

    call pencil_fft_forward

    ! solve phi2
    do k=1,npen
    do j=1,nf
    do i=1,nf*nn+1,2
      ! global grid in Fourier space for i,j,k
      kg=(nn*(icz-1)+icy-1)*npen+k
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
      cxyz((i+1)/2,j,k)=(-4*pi/kr)*cxyz((i+1)/2,j,k)
    enddo
    enddo
    enddo
    if (head) cxyz(1,1,1)=0 ! DC frequency

    call pencil_fft_backward

    open(11,file='phi2.dat',status='replace',access='stream')
    write(11) r3
    close(11)

    phixx=phi(1:nf,1:nf,1:nf) ! backup phi1
    phiyy=r3 ! backup phi2

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
  call destroy_penfft_plan

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

  open(10,file=ic_name('zip0'),status='replace',access='stream')
  open(11,file=ic_name('zip1'),status='replace',access='stream')
  open(12,file=ic_name('zip2'),status='replace',access='stream')
#ifdef PID
  open(14,file=ic_name('zipid'),status='replace',access='stream')
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
      kk=nft*(itz-1)+(ncell/np_nc)*(k-1)+1+imove
      jj=nft*(ity-1)+(ncell/np_nc)*(j-1)+1+imove
      ii=nft*(itx-1)+(ncell/np_nc)*(i-1)+1+imove
      g(1)=ceiling( real(i-1)/np_nc+(0.5+imove)/ncell + (phi(ii-1,jj,kk)-phi(ii+1,jj,kk))/2/(4*pi)/ncell )
      g(2)=ceiling( real(j-1)/np_nc+(0.5+imove)/ncell + (phi(ii,jj-1,kk)-phi(ii,jj+1,kk))/2/(4*pi)/ncell )
      g(3)=ceiling( real(k-1)/np_nc+(0.5+imove)/ncell + (phi(ii,jj,kk-1)-phi(ii,jj,kk+1))/2/(4*pi)/ncell )
      rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1
      !if (i==1 .and. j==1 .and. k==1) then
      !  print*, 4*( real(i-1)/np_nc+0.5/ncell + (phi(ii-1,jj,kk)-phi(ii+1,jj,kk))/2/(4*pi)/ncell )
      !  print*, 4*( real(j-1)/np_nc+0.5/ncell + (phi(ii,jj-1,kk)-phi(ii,jj+1,kk))/2/(4*pi)/ncell )
      !  print*, 4*( real(k-1)/np_nc+0.5/ncell + (phi(ii,jj,kk-1)-phi(ii,jj,kk+1))/2/(4*pi)/ncell )
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
      kk=nft*(itz-1)+(ncell/np_nc)*(k-1)+1+imove
      jj=nft*(ity-1)+(ncell/np_nc)*(j-1)+1+imove
      ii=nft*(itx-1)+(ncell/np_nc)*(i-1)+1+imove
      g(1)=ceiling( real(i-1)/np_nc+(0.5+imove)/ncell + (phi(ii-1,jj,kk)-phi(ii+1,jj,kk))/2/(4*pi)/ncell )
      g(2)=ceiling( real(j-1)/np_nc+(0.5+imove)/ncell + (phi(ii,jj-1,kk)-phi(ii,jj+1,kk))/2/(4*pi)/ncell )
      g(3)=ceiling( real(k-1)/np_nc+(0.5+imove)/ncell + (phi(ii,jj,kk-1)-phi(ii,jj,kk+1))/2/(4*pi)/ncell )
      rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
      idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3))
      x(1,idx)=floor(( real(i-1)/np_nc+(0.5+imove)/ncell + (phi(ii-1,jj,kk)-phi(ii+1,jj,kk))/2/(4*pi)/ncell )/x_resolution)
      x(2,idx)=floor(( real(j-1)/np_nc+(0.5+imove)/ncell + (phi(ii,jj-1,kk)-phi(ii,jj+1,kk))/2/(4*pi)/ncell )/x_resolution)
      x(3,idx)=floor(( real(k-1)/np_nc+(0.5+imove)/ncell + (phi(ii,jj,kk-1)-phi(ii,jj,kk+1))/2/(4*pi)/ncell )/x_resolution)
      v(1,idx)=nint( (phi(ii-1,jj,kk)-phi(ii+1,jj,kk))/2/(4*pi)*vf/v_i2r(1), kind=izipv)
      v(2,idx)=nint( (phi(ii,jj-1,kk)-phi(ii,jj+1,kk))/2/(4*pi)*vf/v_i2r(2), kind=izipv)
      v(3,idx)=nint( (phi(ii,jj,kk-1)-phi(ii,jj,kk+1))/2/(4*pi)*vf/v_i2r(3), kind=izipv)
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
  endfunction

end
