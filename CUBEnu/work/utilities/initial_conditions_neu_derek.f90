!#define READ_NOISE

program initial_conditions
  use pencil_fft
  !use powerspectrum
  use iso_fortran_env, only : int64
  implicit none
  save

  ! nc: coarse grid per image per dim
  ! nf: fine grid per image per dim, ng=nf
  ! nyquest: Nyquest frequency
  logical,parameter :: correct_kernel=.true.
  logical,parameter :: write_potential=.true.

  real,parameter :: vsim2phys_zi_nu=(150./a_i_nu)*box*h0*sqrt(omega_m)/nf_global
  real, parameter :: fdf = 25.8341 !kT/m for T=1K, m=1eV
  real, parameter :: fd = (1./vsim2phys_zi_nu)*fdf*maxval(Tnu/Mnu)/a_i_nu !kBcTnu/mass with temp in K and mass in eV
  real, parameter :: sigma_vi_nu = 3.59714*fd/sqrt(3.) !fd velocity dispersion (45/3 * Zeta(5)/Zeta(3))**0.5

  !Fermi-Dirac CDF
  integer, parameter :: ncdf = 10000
  real, dimension(2,ncdf) :: cdf

  integer(8),parameter :: nk=132 ! ???
  integer(8) i,j,k,n
  integer(4) seedsize
  real kmax,temp_r,temp_theta,pow,phi8,temp8[*]
  real(8) v8, norm, xq(3),gradphi(3),gradphiv(3),vreal(3), dvar[*], dvarg
  integer(int64) :: time64

  integer(8) nplocal[*],npglobal,ip,l,idx_v
  ! power spectrum arrays
  real, dimension(20,nk) :: tf, tfv, tfd !CLASS
  real, dimension(2,nc) :: pkm,pkn

  integer(4),allocatable :: iseed(:)
  real,allocatable :: rseed_all(:,:)

  integer(8) ind,dx,dxy,kg,mg,jg,ig,ii,jj,kk,itx,ity,itz,idx,imove,nf_shake
  integer(8) ileft,iright,nlen,nlast,g(3)
  real kr,kx,ky,kz
  !real xi(10,nbin)

  complex delta_k(nyquest+1,nf,npen)
  real phi(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]
  real phiv(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]

  ! zip arrays
  integer(8),parameter :: npt=nt*np_nc_nu ! np / tile / dim !64
  integer(8),parameter :: npb=ncb*np_nc_nu !24
  integer(8),parameter :: npmax=tile_buffer*(npt+2*npb)**3
  integer(4) rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(4) rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  real(4) vfield(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(8) cume(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(izipx) xp(3,npmax)
  integer(izipv) vp(3,npmax)
  real :: rng(3,npmax)
  real :: rng1(3)
#ifdef EID
    integer(8) pid(npmax)
    integer(8) iq(3)
#endif
  real grad_max(3)[*],gradv_max(3)[*],vmax(3),vf
  !real vdisp(506,2),sigma_vi
  real(4) svz(500,2),svr(100,2)
  real(8) sigma_vc,sigma_vf,sigma_vres
  real(8) std_vsim_c,std_vsim_res,std_vsim

  character (10) :: img_s, z_s

  !equivalence(phixx,phixy)
  !equivalence(phiyy,phiyz)
  !equivalence(phizz,phizx)

  call geometry
  call system('mkdir -p '//opath//'image'//image2str(image))

  if (head) then
    print*, ''
    print*, 'CUBE Initial Conditions -- NEUTRINOS'
    print*, 'on',nn**3,' images'
    print*, 'Resolution', ng*nn
    print*, 'Number of particles per side', np_nc_nu*nc*nn
    print*, 'Box size', box
    print*, 'body_centered_cubic =',body_centered_cubic
    print*, 'output: ', opath
    print*, 'head image number',icx,icy,icz
    print*, '-----------------------------------------'
  endif

  open(10,file=ic_name('info'),status='old',access='stream')
  read(10) sim
  close(10)

  sync all

  ! initialize variables ------------------------------
  if (body_centered_cubic) then
    nf_shake=1
  else
    nf_shake=0
  endif
  ! initvar
  phi=0
  tf=0;tfv=0;tfd=0
  pkm=0
  pkn=0
  sync all

  if (head) print*,'Creating FFT plans'
  call create_penfft_plan
  sync all

  ! transferfnc --------------------------------------
  ! remark: requires "CLASS" format for tf ("CAMB"="CLASS"/(-k^2) with k in 1/Mpc)

  open(11,file='../tf/cafs_z3_dtk.dat',form='formatted')
  read(11,*) tfd(:14,:)
  close(11)

  open(11,file='../tf/cafs_z3_vtk.dat',form='formatted')
  read(11,*) tfv(:14,:)
  close(11)

  tfd(2:3,:)=0
  do i=1,Nnu
     kr= Mnu(i)*(Tnu(i)/Tcnb)**3./Meff
     tfd(2,:)=tfd(2,:)+kr*tfd(5+i-1,:)
     tfd(3,:)=tfd(3,:)+kr*tfv(5+i-1,:) !velocity, units are km/s
  end do
  tfd(3,:) = -tfd(3,:)/vsim2phys_zi_nu

!!$  !!@cdm
!!$  tfd(2,:)=tfd(4,:)
!!$  tfd(3,:)=-tfv(4,:)*vphys2sim

  tfd(3,:)=tfd(3,:)/tfd(2,:) !ratio of v to do

  tf=tfd

  ! compute power spectrum @ z_tf
  tf(2,:) = A_s*(tf(1,:)/k_o)**(n_s-1.)*tf(2,:)**2

  ! propagate to starting redshift
  if (head .and. z_i_nu .ne. z_tf) write(*,*) 'Warning: z_i_nu != z_tf',z_i_nu,z_tf,DgrowRatio(z_i_nu,z_tf)
  tf(2,:) = tf(2,:)*DgrowRatio(z_i_nu,z_tf)**2
  sync all

  ! noisemap -------------------------------------
  if (head) print*,'Generating random noise'
  call random_seed(size=seedsize)
  if (head) print*,'min seedsize =', seedsize
  seedsize=max(seedsize,36)
  allocate(iseed(seedsize))
  allocate(rseed_all(seedsize,nn**3))

!    if (head) print*, 'Copy and read seeds from ../confings/'
!    call system('cp ../configs/seed_'//image2str(image)//'.bin '//opath//'image'//image2str(image))
  open(11,file=output_dir()//'seed'//output_suffix(),status='old',access='stream')
  read(11) iseed
  close(11)
  ! Input iseed
  call random_seed(put=iseed)
  if (head) print*, 'iseed', iseed

  call random_number(r3)
  deallocate(iseed)
  deallocate(rseed_all)
  sync all

# ifdef READ_NOISE
    open(11,file=output_dir()//'noise'//output_suffix(),access='stream')
    read(11) r3
    close(11)
    print*, 'READ IN NOISE MAP:', r3(1,1,1), r3(ng,ng,ng)
!# else
!    open(11,file=output_dir()//'noise'//output_suffix(),status='replace',access='stream')
!    write(11) r3
!    close(11)
!    print*, 'noise',int(image,1),r3(1:2,1,1)
# endif
  sync all

  ! Box-Muller transform ----------------------------------------------
  if (head) print*,'Box-Muller transform'
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
  sync all

  ! delta_field ----------------------------------------------------
  if (head) print*, ''
  if (head) print*, 'Delta field'
  if (head) print*, 'Start ftran'
  call pencil_fft_forward
  if (head) print*, 'Wiener filter on white noise'
  do k=1,npen
  do j=1,nf
  do i=1,nyquest+1
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nyquest-1,nf_global)-nyquest
    ky=mod(jg+nyquest-1,nf_global)-nyquest
    kx=ig-1
    kr=sqrt(kx**2+ky**2+kz**2)
    kr=max(kr,1.0)
    pow=interp_tf(2*pi*kr/box,1,2)/(4*pi*kr**3)
    cxyz(i,j,k)=cxyz(i,j,k)*sqrt(pow*nf_global*nf_global*nf_global)
  enddo
  enddo
  enddo
  if (head) cxyz(1,1,1)=0 ! DC frequency
  sync all
  delta_k=cxyz ! backup k-space delta_L

  if (head) print*,'Start btran'
  call pencil_fft_backward
  !print*,'r3',r3(1,1,1)
  print*,'rms of delta',sqrt(sum(r3**2*1.d0)/nf_global/nf_global/nf_global)

  if (head) print*,'Write delta_L_nu into file'
  if (head) print*,'Growth factor DgrowRatio(',a_i_nu,') =',DgrowRatio(z_i_nu,z_tf)
  open(11,file=output_dir()//'delta_L_nu'//output_suffix(),status='replace',access='stream')
  write(11) r3/DgrowRatio(z_i_nu,z_tf)
  close(11)

  ! Potential field ----------------------------------------------------
  if (head) print*, ''
  if (head) print*, 'Potential field'
  do k=1,npen
  do j=1,nf
  do i=1,nyquest+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nyquest-1,nf_global)-nyquest
    ky=mod(jg+nyquest-1,nf_global)-nyquest
    kx=ig-1
    kz=2*sin(pi*kz/nf_global)
    ky=2*sin(pi*ky/nf_global)
    kx=2*sin(pi*kx/nf_global)
    kr=kx**2+ky**2+kz**2
    kr=max(kr,1.0/nf_global**2) ! avoid kr being 0
    cxyz(i,j,k)=-4*pi/kr
  enddo
  enddo
  enddo
  if (head) cxyz(1,1,1)=0 ! DC frequency
  sync all

  if (correct_kernel) then
    if (head) print*, 'correct kernel'
    call pencil_fft_backward
    temp8=0
    if (image==1) temp8=temp8+r3(9,1,1)+r3(1,9,1)+r3(1,1,9)
    sync all
    if (icx==nn .and. icy==1 .and. icz==1) temp8=temp8+r3(nf-7,1,1)
    sync all
    if (icx==1 .and. icy==nn .and. icz==1) temp8=temp8+r3(1,nf-7,1)
    sync all
    if (icx==1 .and. icy==1 .and. icz==nn) temp8=temp8+r3(1,1,nf-7)
    sync all
    phi8=0
    do i=1,nn**3
      phi8=phi8+temp8[i]
    enddo
    sync all
    phi8=phi8/6
    if (head) print*,'phi8 =',phi8
    sync all
    if (head) print*, 'Construct Ewald potential kernel in real space'
    do k=1,nf
    do j=1,nf
    do i=1,nf
      kg=k+nf*(icz-1)
      jg=j+nf*(icy-1)
      ig=i+nf*(icx-1)
      kx=mod(kg+nyquest-1,nf_global)-nyquest
      ky=mod(jg+nyquest-1,nf_global)-nyquest
      kz=mod(ig+nyquest-1,nf_global)-nyquest
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

  ! Complex multiply delta_L with potential kernel
  cxyz=real(cxyz)*delta_k

  call pencil_fft_backward

  ! Primordial Non-Gaussianity
  dvar=sum((r3*1d0)**2)
  dvarg=0
  do i=1,nn**3 ! co_sum
    dvarg=dvarg+dvar[i]
  enddo
  dvarg=dvarg/nf_global/nf_global/nf_global

  r3=r3-f_nl*(r3**2-dvarg)

  phi=0
  phi(1:nf,1:nf,1:nf)=r3 ! phi1
  if (head) print*, 'Write phi1_nu into file'
  open(11,file=ic_name_nu('phi1_nu'),status='replace',access='stream')
  write(11) r3
  close(11)

  sync all

  ! buffer phi ---------------------------------------------------
  if (head) print*, 'Buffer phi'
  phi(:0,:,:)=phi(nf-nfb:nf,:,:)[image1d(inx,icy,icz)]
  phi(nf+1:,:,:)=phi(1:nfb+1,:,:)[image1d(ipx,icy,icz)]
  sync all
  phi(:,:0,:)=phi(:,nf-nfb:nf,:)[image1d(icx,iny,icz)]
  phi(:,nf+1:,:)=phi(:,1:nfb+1,:)[image1d(icx,ipy,icz)]
  sync all
  phi(:,:,:0)=phi(:,:,nf-nfb:nf)[image1d(icx,icy,inz)]
  phi(:,:,nf+1:)=phi(:,:,1:nfb+1)[image1d(icx,icy,ipz)]
  sync all

  ! potential field v ------------------------------------------------
  if (head) write(*,*) 'Computing velocity potential'

  if (head) write(*,*) 'ftran'
  call pencil_fft_forward
  !! ? cxyz=delta_k ! use backed-up delta_L in k-space
  do k=1,npen
  do j=1,nf
  do i=1,nyquest+1
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nyquest-1,nf_global)-nyquest
    ky=mod(jg+nyquest-1,nf_global)-nyquest
    kx=ig-1
    kr=sqrt(kx**2+ky**2+kz**2)
    kr=max(kr,1.0)
    pow=interp_tf(2*pi*kr/box,1,3)!/(4*pi*kr**3)
    cxyz(i,j,k)=cxyz(i,j,k)*pow*kr
  enddo
  enddo
  enddo
  sync all

  if (head) write(*,*) 'btran'
  call pencil_fft_backward
  sync all

  phiv=0
  phiv(1:nf,1:nf,1:nf)=r3
  sync all
  if (head) print*, 'Write phiv_nu into file'
  open(11,file=ic_name_nu('phiv_nu'),status='replace',access='stream')
  write(11) r3
  close(11)

  if (head) write(*,*) 'buffer'
  phiv(:0,:,:)=phiv(nf-nfb:nf,:,:)[image1d(inx,icy,icz)]
  phiv(nf+1:,:,:)=phiv(1:nfb+1,:,:)[image1d(ipx,icy,icz)]
  sync all
  phiv(:,:0,:)=phiv(:,nf-nfb:nf,:)[image1d(icx,iny,icz)]
  phiv(:,nf+1:,:)=phiv(:,1:nfb+1,:)[image1d(icx,ipy,icz)]
  sync all
  phiv(:,:,:0)=phiv(:,:,nf-nfb:nf)[image1d(icx,icy,inz)]
  phiv(:,:,nf+1:)=phiv(:,:,1:nfb+1)[image1d(icx,icy,ipz)]
  sync all

  ! destroy fft ------------------------------------------------

  if (head) print*, 'Destroying FFT plans'
  call destroy_penfft_plan
  sync all

  ! zip checkpoints ------------------------------------------------
  if (head) print*, 'zip checkpoints'
  vf=1.!vfactor(a)
  if (head) print*, 'vf',vf
  grad_max(1)=maxval(abs(phi(0:nf-1,1:nf,1:nf)-phi(2:nf+1,1:nf,1:nf)))
  grad_max(2)=maxval(abs(phi(1:nf,0:nf-1,1:nf)-phi(1:nf,2:nf+1,1:nf)))
  grad_max(3)=maxval(abs(phi(1:nf,1:nf,0:nf-1)-phi(1:nf,1:nf,2:nf+1)))
  gradv_max(1)=maxval(abs(phiv(0:nf-1,1:nf,1:nf)-phiv(2:nf+1,1:nf,1:nf)))
  gradv_max(2)=maxval(abs(phiv(1:nf,0:nf-1,1:nf)-phiv(1:nf,2:nf+1,1:nf)))
  gradv_max(3)=maxval(abs(phiv(1:nf,1:nf,0:nf-1)-phiv(1:nf,1:nf,2:nf+1)))
  sync all
  do i=1,nn**3 ! co_max
    grad_max=max(grad_max,grad_max(:)[i])
    gradv_max=max(gradv_max,gradv_max(:)[i])
  enddo
  !print*, 'grad_max',grad_max
  !print*, 'gradv_max',gradv_max
  !stop
  vmax=gradv_max/2/(4*pi)*vf
  sim%dt_vmax_nu=vbuf*20./maxval(abs(vmax))
  if (head) then
    print*, 'grad_max',grad_max
    print*, 'max dsp',grad_max/2/(4*pi)
    print*, 'vmax',vmax
    if (maxval(grad_max)/2/(4*pi)>=nfb) then
      print*, 'particle dsp > buffer'
      print*, maxval(grad_max)/2/(4*pi),nfb
      stop
    endif
  endif
  sync all

!!$  open(11,file='../velocity_conversion/sigmav_z.bin',access='stream')
!!$  read(11) svz
!!$  close(11)
!!$  open(11,file='../velocity_conversion/sigmav_r.bin',access='stream')
!!$  read(11) svr
!!$  close(11)

!!$  sigma_vf=interp_sigmav(a,box/nf_global) ! sigma(v) on scale of fine grid, in km/s
!!$  sigma_vc=interp_sigmav(a,box/nc_global) ! sigma(v) on scale of coarse grid, in km/s
!!$  sigma_vres=sqrt(sigma_vf**2-sigma_vc**2) ! sigma(v) residual, in km/s
!!$  sigma_vi=sigma_vres/sim%vsim2phys/sqrt(3.) ! sigma(v_i) residual, in sim unit

  !sigma_vf=sigma_vi_nu/vphys2sim
  !sigma_vc=sigma_vf/vphys2sim
  !sigma_vres=0.
  !sigma_vi=0.

  !sim%sigma_vres=sigma_vres
  !sim%sigma_vi=sigma_vi
  !sync all

  if (head) then
    print*,''
    print*,'Theory velocity dispersion prediction'
    print*,' =',real(sigma_vi_nu,4)*vsim2phys_zi_nu*sqrt(3.),'km/s'
    print*,'sigma_vi_nu =',real(sigma_vi_nu,4),'(simulation unit)'
    print*,''
  endif
  sim%sigma_vi_nu=sigma_vi_nu
  sync all

  if (head) write(*,*) 'Generating cdf'
  call compute_cdf
  if (head) then
     write(*,*) 'Writing CDF to file'
     write(*,*) 'FD factor used: ',fd
     open(11,file=opath//'cdf.txt')
     do i=1,ncdf
        write(11,*) cdf(1,i)/fd,cdf(2,i)
     end do
     close(11)
  end if

  ! create particles (no communication) ----------------------------
  if (head) print*,''
  if (head) print*, 'Create particles'
  open(11,file=ic_name('xp_nu'),status='replace',access='stream')
  open(12,file=ic_name('vp_nu'),status='replace',access='stream')
  open(13,file=ic_name('np_nu'),status='replace',access='stream')
  open(14,file=ic_name('vc_nu'),status='replace',access='stream')
#ifdef EID
  if (head) print*, '  also create EID'
  open(15,file=ic_name('id_nu'),status='replace',access='stream')
#endif

  vfield=0
  nplocal=0
  std_vsim_c=0; std_vsim_res=0; std_vsim=0;
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    print*, 'working on tile',itx,ity,itz

    !Neutrino random velocities
    call random_number(rng)
    do k=1,npmax
       n=nearest_loc(rng(1,k),cdf(2,:))
       kx=interp(rng(1,k),cdf(1,n),cdf(1,n+1),cdf(2,n),cdf(2,n+1))
       ky=2.*rng(2,k)-1.
       kz=2.*pi*rng(3,k)
       rng(1,k)=kx*sqrt(1.-ky**2.)*cos(kz)
       rng(2,k)=kx*sqrt(1.-ky**2.)*sin(kz)
       rng(3,k)=kx*ky
    enddo

    idx_v=0
    iright=0
    rhoce=0
    rholocal=0
    do k=1-npb,npt+npb ! calculate coarse mesh density
    do j=1-npb,npt+npb
    do i=1-npb,npt+npb
    do imove=0,merge(1,0,body_centered_cubic)
      kk=nft*(itz-1)+(ncell/np_nc_nu)*(k-1)+1+imove*(ncell/np_nc_nu/2)
      jj=nft*(ity-1)+(ncell/np_nc_nu)*(j-1)+1+imove*(ncell/np_nc_nu/2)
      ii=nft*(itx-1)+(ncell/np_nc_nu)*(i-1)+1+imove*(ncell/np_nc_nu/2)
      xq=((/i,j,k/)-1d0)/np_nc_nu + (0.5d0+imove*(ncell/np_nc_nu/2))/ncell ! Lagrangian position q

      gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
      gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
      gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)

      gradphiv(1)=phiv(ii+1,jj,kk)-phiv(ii-1,jj,kk)
      gradphiv(2)=phiv(ii,jj+1,kk)-phiv(ii,jj-1,kk)
      gradphiv(3)=phiv(ii,jj,kk+1)-phiv(ii,jj,kk-1)

      g=ceiling(xq-gradphi/(8*pi*ncell))
      rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1
      vreal=-gradphiv/(8*pi)*vf
      vfield(:,g(1),g(2),g(3))=vfield(:,g(1),g(2),g(3))+vreal ! record vfield according to real particles
    enddo
    enddo
    enddo
    enddo

    vfield(1,:,:,:)=vfield(1,:,:,:)/merge(1,rhoce,rhoce==0)
    vfield(2,:,:,:)=vfield(2,:,:,:)/merge(1,rhoce,rhoce==0)
    vfield(3,:,:,:)=vfield(3,:,:,:)/merge(1,rhoce,rhoce==0)
    cume=cumsum3(rhoce)

    if (body_centered_cubic .and. ncell/np_nc_nu/2==0) stop 'ncell/np_nc_nu/2 = 0, unsuitable for body centered cubic'
    do k=1-npb,npt+npb ! create particles in extended mesh
    do j=1-npb,npt+npb
    do i=1-npb,npt+npb
    do imove=0,merge(1,0,body_centered_cubic)
      idx_v=idx_v+1
      kk=nft*(itz-1)+(ncell/np_nc_nu)*(k-1)+1+imove*(ncell/np_nc_nu/2)
      jj=nft*(ity-1)+(ncell/np_nc_nu)*(j-1)+1+imove*(ncell/np_nc_nu/2)
      ii=nft*(itx-1)+(ncell/np_nc_nu)*(i-1)+1+imove*(ncell/np_nc_nu/2)
      xq=((/i,j,k/)-1d0)/np_nc_nu + (0.5d0+imove*(ncell/np_nc_nu/2))/ncell ! Lagrangian position q

      gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
      gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
      gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)
      gradphiv(1)=phiv(ii+1,jj,kk)-phiv(ii-1,jj,kk)
      gradphiv(2)=phiv(ii,jj+1,kk)-phiv(ii,jj-1,kk)
      gradphiv(3)=phiv(ii,jj,kk+1)-phiv(ii,jj,kk-1)

      g=ceiling(xq-gradphi/(8*pi*ncell))
      rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
      idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3))
      xp(:,idx)=floor((xq-gradphi/(8*pi*ncell))/x_resolution_nu,kind=8)
      rng1=rng(:,idx_v)
      vreal=rng1-gradphiv/(8*pi)*vf-vfield(:,g(1),g(2),g(3)) ! save relative velocity
      !if (imove==0) then
      !   !pos rng dir
      !   rng1=rng(:,idx)
      !   vreal=-gradphiv/(8*pi)*vf+rng1
      !else
      !   !true random velocities
      !   !rng1=rng(:,idx)
      !   !neg rng dir, conserves momentum
      !   rng1=-1.*rng1

      !   vreal=-gradphiv/(8*pi)*vf+rng1
      !end if
      vp(:,idx)=nint(real(nvbin_nu-1)*atan(sqrt(pi/2)/(sigma_vi_nu*vrel_boost)*vreal)/pi,kind=izipv)
      !print*, vreal
      !print*, sigma_vi_nu
      !print*, vp(:,idx)
#     ifdef EID
        iq = ((/icx,icy,icz/)-1)*nf + ((/itx,ity,itz/)-1)*nft + (ncell/np_nc_nu)*((/i,j,k/)-1)+imove
        iq = modulo(iq,nf_global)
        pid(idx)=iq(1)+nf_global*iq(2)+nf_global**2*iq(3)+1
#     endif
    enddo
    enddo
    enddo
    enddo

    do k=1,nt ! delete buffer particles
    do j=1,nt
      ileft=iright+1
      nlast=cume(nt,j,k)
      nlen=nlast-cume(0,j,k)
      iright=ileft+nlen-1
      xp(:,ileft:iright)=xp(:,nlast-nlen+1:nlast)
      vp(:,ileft:iright)=vp(:,nlast-nlen+1:nlast)
#     ifdef EID
        pid(ileft:iright)=pid(nlast-nlen+1:nlast)
#     endif
    enddo
    enddo

    ! velocity analysis
    ip=0
    do k=1,nt
    do j=1,nt
    do i=1,nt
      std_vsim_c=std_vsim_c+sum(vfield(:,i,j,k)**2)
      do l=1,rhoce(i,j,k)
        ip=ip+1
        vreal=tan(pi*real(vp(:,ip))/real(nvbin_nu-1))/(sqrt(pi/2)/(sigma_vi_nu*vrel_boost))
        std_vsim_res=std_vsim_res+sum(vreal**2)
        vreal=vreal+vfield(:,i,j,k)
        std_vsim=std_vsim+sum(vreal**2)
      enddo
    enddo
    enddo
    enddo

    write(11) xp(:,1:iright)
    write(12) vp(:,1:iright)
    write(13) rhoce(1:nt,1:nt,1:nt)
    write(14) vfield(:,1:nt,1:nt,1:nt)
#   ifdef EID
      write(15) pid(1:iright)
#   endif
    nplocal=nplocal+iright
    !stop 'here'
  enddo
  enddo
  enddo ! end of tile loop

  close(11)
  close(12)
  close(13)
  close(14)
# ifdef EID
    close(15)
# endif

  sim%nplocal_nu=nplocal

  if (head) then
    print*,''
    print*,'Velocity analysis on head node'
    std_vsim_res=sqrt(std_vsim_res/nplocal)
    std_vsim_c=sqrt(std_vsim_c/nc/nc/nc)
    std_vsim=sqrt(std_vsim/nplocal)
    print*,'  vsim2phys_zi_nu ',vsim2phys_zi_nu
    print*,'  std_vsim         ',real(std_vsim*sim%vsim2phys,4),'km/s'
    print*,'  std_vsim_c       ',real(std_vsim_c*sim%vsim2phys,4),'km/s'
    print*,'  std_vsim_res     ',real(std_vsim_res*sim%vsim2phys,4),'km/s'
    print*,'  std_vi (sim unit)',real(std_vsim_res/sqrt(3.),4),'(simulation unit)'
    print*,''
  endif
  sync all

  print*,'image',image,', nplocal',nplocal
  sync all
  npglobal=0
  do i=1,nn**3
    npglobal=npglobal+nplocal[i]
  enddo
  if (head) print*, 'npglobal =',npglobal
  sync all
  sim%npglobal_nu=npglobal
  !sim%mass_p=real(nf_global**3,kind=8)/npglobal
  call print_header(sim)

  sync all
  open(10,file=ic_name('info'),status='replace',access='stream')
  write(10) sim
  close(10)
  sync all
  if (head) print*, 'initial condition done'

  contains

  function cumsum3(input)
    implicit none
    integer(4) input(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8) cumsum3(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8) nsum,igx,igy,igz
    nsum=0
    do igz=1-2*ncb,nt+2*ncb
    do igy=1-2*ncb,nt+2*ncb
    do igx=1-2*ncb,nt+2*ncb
      nsum=nsum+input(igx,igy,igz)
      cumsum3(igx,igy,igz)=nsum
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
  endfunction

  real function interp_tf(kr,ix,iy)
    implicit none
    integer(4) ix,iy
    integer(8) ii,i1,i2
    real kr,xx,yy,x1,x2,y1,y2
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
  endfunction


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

  function DgrowRatio(z1,z2) result(Dgrow)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real, parameter :: np = -(1./4.)+(5./4.)*sqrt(1-24.*omega_nu/omega_m/25.) !~1-3f/5
    real z1,z2
    real Dgrow
    real hsq,oma,ola,a1,a2,ga1,ga2

    a1=1./(1.+z1)
    hsq=om/a1**3+(1-om-ol)/a1**2+ol
    oma=om/(a1**3*hsq)
    ola=ol/hsq
    ga1=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    a2=1./(1.+z2)
    hsq=om/a2**3+(1-om-ol)/a2**2+ol
    oma=om/(a2**3*hsq)
    ola=ol/hsq
    ga2=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    Dgrow=(a1*ga1)/(a2*ga2)
    Dgrow=Dgrow**np!(1.-3.*(omega_nu/omega_m)/5.)
  end function DgrowRatio

  function vfactor(a)
    implicit none
    real, parameter :: np = -(1./4.)+(5./4.)*sqrt(1-24.*omega_nu/omega_m/25.) !~1-3f/5
    real :: a
    real :: H,km,lm
    real :: vfactor
    lm=omega_l/omega_m
    km=(1-omega_m-omega_l)/omega_m
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H
    vfactor=vfactor*np!(1.-3.*(omega_nu/omega_m)/5.)
  endfunction vfactor

  function lcg(s) !// Linear congruential generator
    implicit none
    integer(4) :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  endfunction

  function approxCDF(v) result(c)
    implicit none
    real, intent(in) :: v
    real :: c
    real, parameter :: s=3.5
    c=1.-exp(-(v/s)**2.)
  end function approxCDF

  function invertCDF(c) result(v)
    implicit none
    real, intent(in) :: c
    real :: v
    real, parameter :: s=3.5
    v=s*sqrt(log(1./(1.-c)))
  end function invertCDF

  subroutine compute_cdf
    implicit none
    real(8), dimension(2,ncdf) :: cdf0
    real, dimension(2,ncdf) :: cdfn
    integer, parameter :: ni = 1000 !How many points to integrate per cdf
    real, dimension(ni) :: x,y
    real, parameter :: maxu = 15.0
    real, parameter :: cdfinf = 1.80309
    integer :: i,j,n
    real :: l,u,fnu,fdnu

    cdf0 = 0.
    do i=2,ncdf
       !Limits of integration
       l=maxu*(1.0*i-1.)/ncdf
       u=maxu*(1.0*i)/ncdf
       cdf0(1,i)=u
       !Integral
       do j=1,ni
          x(j)=l+(j-1)*(u-l)/(ni-1)
          y(j)=f0(x(j))
       enddo
       cdf0(2,i)=cdf0(2,i-1)+integrate(x,y)

    enddo

    write(*,*) 'cdf: u->inf = ',cdf0(2,ncdf),cdfinf
    cdf0(2,:) = cdf0(2,:)/cdf0(2,ncdf)

    !Now sum up for each neutrino
    cdf=0
    cdf(1,:)=fd*cdf0(1,:)
    cdfn(2,:)=cdf0(2,:)
    do n=1,Nnu
       fnu=Mnu(n)*(Tnu(n)/Tcnb)**3/Meff ! fraction of energy in this neutrino
       fdnu=(1./vsim2phys_zi_nu)*fdf*Tnu(n)/Mnu(n)/a_i_nu
       cdfn(1,:)=cdf0(1,:)*fdnu
       do i=1,ncdf
          j=nearest_loc(cdf(1,i),cdfn(1,:))
          cdf(2,i) = cdf(2,i)+fnu*interp(cdf(1,i),cdfn(1,j),cdfn(1,j+1),cdfn(2,j),cdfn(2,j+1))
       enddo
    enddo
  endsubroutine compute_cdf

  function f0(u) result(f)
    implicit none
    real, intent(in) :: u
    real :: f
    f= u**2./(exp(u)+1.)
  end function f0

  function integrate(x,y) result(s)
    implicit none
    real, dimension(:), intent(in) :: x,y
    real :: s
    integer :: i
    s=0
    do i=2,size(x)
       s=s+0.5*(x(i)-x(i-1))*(y(i)+y(i-1))
    enddo
  end function integrate

  function nearest_loc(u,c) result(nl)
    implicit none
    real, intent(in) :: u
    real, dimension(:) :: c
    integer :: b1,b2,b,nl
    b1=1
    b2=size(c)
    do while(b2-b1>1)
       b=(b1+b2)/2
       if ( u.gt.c(b) ) then
          b1=b
       else
          b2=b
       endif
    enddo
    nl=merge(b1,b2,b1<b2)
  end function nearest_loc

  function interp(x,x1,x2,y1,y2) result(y)
    implicit none
    real, intent(in) :: x,x1,x2,y1,y2
    real :: y
    y = (y1*(x2-x)+y2*(x-x1))/(x2-x1)
  end function interp

end program initial_conditions
