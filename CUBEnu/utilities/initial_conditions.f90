#define READ_SEED
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

  real,parameter :: a=1/(1+z_i)
  real,parameter :: Vphys2sim=1.0/(300.*sqrt(omega_m)*box/a/2/nf_global)
  integer(8),parameter :: nk=132 ! ???
  integer(8) i,j,k
  integer(4) seedsize
  real kmax,temp_r,temp_theta,pow,phi8,temp8[*]
  real(8) v8, norm, xq(3),gradphi(3),vreal(3), dvar[*], dvarg
  integer(int64) :: time64

  integer(8) nplocal[*],npglobal,ip,l
  ! power spectrum arrays
  real, dimension(14,nk) :: tf    !CAMB ! ???
  real, dimension(2,nc) :: pkm,pkn

  integer(4),allocatable :: iseed(:)
  real,allocatable :: rseed_all(:,:)

  integer(8) ind,dx,dxy,kg,mg,jg,ig,ii,jj,kk,itx,ity,itz,idx,imove,nf_shake
  integer(8) ileft,iright,nlen,nlast,g(3)
  real kr,kx,ky,kz
  !real xi(10,nbin)

  complex delta_k(nyquest+1,nf,npen)
  real phi(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]

  ! zip arrays
  integer(8),parameter :: npt=nt*np_nc ! np / tile / dim !64
  integer(8),parameter :: npb=ncb*np_nc !24
  integer(8),parameter :: npmax=2*(npt+2*npb)**3
  integer(4) rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(4) rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  real(4) vfield(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(8) cume(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(izipx) xp(3,npmax)
  integer(izipv) vp(3,npmax)
#ifdef PID
    integer(8) pid(npmax)
    integer(8) iq(3)
#endif
  real grad_max(3)[*],vmax(3),vf
  !real vdisp(506,2),sigma_vi
  real(4) svz(500,2),svr(100,2)
  real(8) sigma_vc,sigma_vf,sigma_vres,sigma_vi
  real(8) std_vsim_c,std_vsim_res,std_vsim

  character (10) :: img_s, z_s

  !equivalence(phixx,phixy)
  !equivalence(phiyy,phiyz)
  !equivalence(phizz,phizx)

  call geometry

  if (head) then
    print*, ''
    print*, 'CUBE Initial Conditions'
    print*, 'on',nn**3,' images'
    print*, 'Resolution', ng*nn
    print*, 'Number of particles per side', np_nc*nc*nn
    print*, 'Box size', box
    print*, 'np_2n3 =',np_2n3
    print*, 'output: ', opath
    print*, 'head image number',icx,icy,icz
    print*, '-----------------------------------------'
  endif
  sync all

  call system('mkdir -p '//opath//'image'//image2str(image))

  sim%nplocal=0
  sim%nplocal_nu=0
  sim%a=1./(1+z_i)
  sim%t=0
  sim%tau=0

  sim%istep=0

  sim%dt_f_acc=1000
  sim%dt_pp_acc=1000
  sim%dt_c_acc=1000

  sim%cur_checkpoint=0
  sim%cur_proj=0
  sim%cur_halo=0

  sim%mass_p=real(nf**3)/sim%nplocal ! will be overwritten

  sim%box=box
  sim%image=image
  sim%nn=nn
  sim%nnt=nnt
  sim%nt=nt
  sim%ncell=ncell
  sim%ncb=ncb
  sim%izipx=izipx
  sim%izipv=izipv

  sim%h0=h0
  sim%omega_m=omega_m
  sim%omega_l=omega_l
  sim%s8=s8
  sim%vsim2phys=(1.5/a)*box*h0*sqrt(omega_m)/nf_global
  sim%z_i=z_i
  sim%z_i_nu=z_i_nu
  sync all

  ! initialize variables ------------------------------
  if (np_2n3) then
    nf_shake=1
  else
    nf_shake=0
  endif
  ! initvar
  phi=0
  tf=0
  pkm=0
  pkn=0
  sync all

  if (head) print*,'Creating FFT plans'
  call create_penfft_plan
  sync all

  ! transferfnc --------------------------------------
  ! remark: requires "CLASS" format for tf ("CAMB"="CLASS"/(-k^2) with k in 1/Mpc)
  open(11,file='../tf/caf_z10_tk.dat',form='formatted')
  read(11,*) !header
  read(11,*) tf
  close(11)

  ! replace T_g with T_cb = f_c T_c + f_b T_b
  tf(2,:) = (omega_bar*tf(3,:)+omega_cdm*tf(4,:))/(omega_bar+omega_cdm)

  ! compute power spectrum @ z_tf
  tf(2,:) = A_s*(tf(1,:)/k_o)**(n_s-1.)*tf(2,:)**2

  ! propagate to starting redshift
  tf(2,:) = tf(2,:)*DgrowRatio(z_i,z_tf)**2

  sync all

!!$  ! transferfnc --------------------------------------
!!$  !open(11,file='../tf/ith2_nu0p05_z5_tk.dat',form='formatted')
!!$  open(11,file='../tf/nu100_onu3/nu100_onu3_transfer_out_z10.dat',form='formatted')
!!$  !open(11,file='../configs/mmh_transfer/simtransfer_bao.dat',form='formatted') ! for Xin
!!$  read(11,*) tf
!!$  close(11)
!!$  ! normalization
!!$  ! norm=2.*pi**2.*(h0/100.)**4*(h0/100./0.05)**(n_s-1) ! for Xin
!!$  norm=1
!!$  if (head) print*, 'Normalization factor: norm =', norm
!!$  !Delta^2
!!$  do i=2,size(tf,dim=1)
!!$     tf(i,:)=tf(i,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
!!$  end do
!!$!  tf(2,:)=tf(2,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
!!$!  tf(3,:)=tf(3,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
!!$!  tf(6,:)=tf(6,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
!!$  !dk
!!$  tf(4,1)=tf(1,2)/2
!!$  do k=2,nk-1
!!$    tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
!!$  enddo
!!$  tf(4,nk)=tf(1,nk)-tf(1,nk-1)
!!$  v8=0
!!$  kmax=2*pi*sqrt(3.)*nyquest/box
!!$  do k=1,nk
!!$    if (tf(1,k)>kmax) exit
!!$    v8=v8+tf(7,k)*tophat(tf(1,k)*8)**2*tf(4,k)/tf(1,k)
!!$  enddo
!!$!  if (head) print*, 's8**2/v8:', v8, s8**2/v8,nyquest ;stop
!!$!  tf(2:3,:)=tf(2:3,:)*(s8**2/v8)*Dgrow(a)**2
!!$  tf(2,:)=tf(6,:)*(s8**2/v8)*DgrowRatio(z_i,z_tf)**2 ! T_cb rather than T_c
!!$!  tf(2:3,:)= scalar_amp*tf(2:3,:)*Dgrow(a)**2 ! for Xin
!!$  sync all

!print*, tf(1,:)
!print*, ''
!print*, tf(2,:)
!stop

  ! noisemap -------------------------------------
  if (head) print*,'Generating random noise'
  call random_seed(size=seedsize)
  if (head) print*,'min seedsize =', seedsize
  seedsize=max(seedsize,36)
  allocate(iseed(seedsize))
  allocate(rseed_all(seedsize,nn**3))
#ifdef READ_SEED
    if (head) print*, 'Copy and read seeds from ../confings/'
    call system('cp ../configs/seed_'//image2str(image)//'.bin '//opath//'image'//image2str(image))
    open(11,file=output_dir()//'seed'//output_suffix(),status='old',access='stream')
    read(11) iseed
    close(11)
    ! Input iseed
    call random_seed(put=iseed)
    if (head) print*, 'iseed', iseed
#else
    ! Generate at least 12 seeds according to system clock
    call system_clock(time64)
    do i = 1, seedsize
      iseed(i) = lcg(time64) + image*137
      !print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
    enddo
    ! Input iseed to system
    call random_seed(put=iseed)
    !print*, 'iseed', iseed
    ! Write iseed into file
    open(11,file=output_dir()//'seed'//output_suffix(),status='replace',access='stream')
    write(11) iseed
    close(11)
    ! execute the following line if you want to save seeds
    !call system('cp ../output/universe1/image*/seed* ../configs/')
#endif

  call random_number(r3)
  deallocate(iseed)
  deallocate(rseed_all)
  sync all

# ifdef READ_NOISE
    open(11,file=output_dir()//'noise'//output_suffix(),access='stream')
    read(11) r3
    close(11)
    print*, 'READ IN NOISE MAP:', r3(1,1,1), r3(ng,ng,ng)
# else
    open(11,file=output_dir()//'noise'//output_suffix(),status='replace',access='stream')
    write(11) r3
    close(11)
    print*, 'noise',int(image,1),r3(1:2,1,1)
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
  !print*,'rms of delta',sqrt(sum(r3**2*1.d0)/nf_global/nf_global/nf_global)

  if (head) print*,'Write delta_L into file'
  if (head) print*,'Growth factor Dgrow(',a,') =',Dgrow(a)
  open(11,file=output_dir()//'delta_L'//output_suffix(),status='replace',access='stream')
  write(11) r3/Dgrow(a)
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
  delta_k=cxyz  ! backup phi(k)
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
  if (head) print*, 'Write phi1 into file'
  open(11,file=ic_name('phi1'),status='replace',access='stream')
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

  if (head) print*, 'Destroying FFT plans'
  call destroy_penfft_plan
  sync all

  ! zip checkpoints ------------------------------------------------
  if (head) print*, 'zip checkpoints'
  vf=vfactor(a)
  if (head) print*, 'vf',vf
  grad_max(1)=maxval(abs(phi(0:nf-1,1:nf,1:nf)-phi(2:nf+1,1:nf,1:nf)))
  grad_max(2)=maxval(abs(phi(1:nf,0:nf-1,1:nf)-phi(1:nf,2:nf+1,1:nf)))
  grad_max(3)=maxval(abs(phi(1:nf,1:nf,0:nf-1)-phi(1:nf,1:nf,2:nf+1)))
  sync all
  do i=1,nn**3 ! co_max
    grad_max=max(grad_max,grad_max(:)[i])
  enddo
  !print*, grad_max
  vmax=grad_max/2/(4*pi)*vf
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

  open(23,file='../velocity_conversion/sigmav_z.bin',access='stream')
  read(23) svz
  close(23)
  open(24,file='../velocity_conversion/sigmav_r.bin',access='stream')
  read(24) svr
  close(43)

  sigma_vf=interp_sigmav(a,box/nf_global) ! sigma(v) on scale of fine grid, in km/s
  sigma_vc=interp_sigmav(a,box/nc_global) ! sigma(v) on scale of coarse grid, in km/s
  sigma_vres=sqrt(sigma_vf**2-sigma_vc**2) ! sigma(v) residual, in km/s
  sigma_vi=sigma_vres/sim%vsim2phys/sqrt(3.) ! sigma(v_i) residual, in sim unit

  sim%sigma_vres=sigma_vres
  sim%sigma_vi=sigma_vi
  sync all
  if (head) then
    print*, ''
    print*,'sigma_vf(a=',a,', r=',box/nf_global,'Mpc/h)=',sigma_vf,'km/s'
    print*,'sigma_vc(a=',a,', r=',box/nc_global,'Mpc/h)=',sigma_vc,'km/s'
    print*,'sigma_vres=',sigma_vres,'km/s'
    print*,'sigma_vi =',sigma_vi,'(simulation unit)'
  endif


  sync all

  ! create particles (no communication) ----------------------------
  if (head) print*,''
  if (head) print*, 'Create particles'
  open(11,file=ic_name('xp'),status='replace',access='stream')
  open(12,file=ic_name('vp'),status='replace',access='stream')
  open(13,file=ic_name('np'),status='replace',access='stream')
  open(14,file=ic_name('vc'),status='replace',access='stream')
#ifdef PID
  if (head) print*, '  also create PID'
  open(15,file=ic_name('id'),status='replace',access='stream')
#endif


  vfield=0
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
      xq=((/i,j,k/)-1d0)/np_nc + (0.5d0+imove)/ncell ! Lagrangian position q
      gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
      gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
      gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)
      g=ceiling(xq-gradphi/(8*pi*ncell))
      rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1
      vreal=-gradphi/(8*pi)*vf
      vfield(:,g(1),g(2),g(3))=vfield(:,g(1),g(2),g(3))+vreal ! record vfield according to real particles
    enddo
    enddo
    enddo
    enddo
    vfield(1,:,:,:)=vfield(1,:,:,:)/rhoce
    vfield(2,:,:,:)=vfield(2,:,:,:)/rhoce
    vfield(3,:,:,:)=vfield(3,:,:,:)/rhoce

    cume=cumsum3(rhoce)

    do k=1-npb,npt+npb ! create particles in extended mesh
    do j=1-npb,npt+npb
    do i=1-npb,npt+npb
    do imove=0,nf_shake
      kk=nft*(itz-1)+(ncell/np_nc)*(k-1)+1+imove
      jj=nft*(ity-1)+(ncell/np_nc)*(j-1)+1+imove
      ii=nft*(itx-1)+(ncell/np_nc)*(i-1)+1+imove
      xq=((/i,j,k/)-1d0)/np_nc + (0.5d0+imove)/ncell
      gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
      gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
      gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)
      g=ceiling(xq-gradphi/(8*pi*ncell))
      rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
      idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3))
      xp(:,idx)=floor((xq-gradphi/(8*pi*ncell))/x_resolution,kind=8)
      vreal=-gradphi/(8*pi)*vf
      vreal=vreal-vfield(:,g(1),g(2),g(3)) ! save relative velocity
      vp(:,idx)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
#     ifdef PID
        iq = ((/icx,icy,icz/)-1)*nf + ((/itx,ity,itz/)-1)*nft + (ncell/np_nc)*((/i,j,k/)-1)+imove
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
#     ifdef PID
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
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
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
#   ifdef PID
      write(15) pid(1:iright)
#   endif
    nplocal=nplocal+iright

  enddo
  enddo
  enddo ! end of tile loop

  close(11)
  close(12)
  close(13)
  close(14)
# ifdef PID
    close(15)
# endif

  sim%nplocal=nplocal

  if (head) then
    print*,''
    print*,'Velocity analysis on head node'
    std_vsim_res=sqrt(std_vsim_res/nplocal)
    std_vsim_c=sqrt(std_vsim_c/nc/nc/nc)
    std_vsim=sqrt(std_vsim/nplocal)
    print*,'  std_vsim',std_vsim*sim%vsim2phys,'km/s'
    print*,'  std_vsim_c',std_vsim_c*sim%vsim2phys,'km/s'
    print*,'  std_vsim_res',std_vsim_res*sim%vsim2phys,'km/s'
    print*,'  std_vi (sim unit)',std_vsim_res/sqrt(3.),'(simulation unit)'
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
  sim%mass_p=real(nf_global**3,kind=8)/npglobal
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
end
