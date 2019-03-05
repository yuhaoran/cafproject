!#define specify_potential
program lpt
!  use omp_lib
  use variables, only: spine_tile
  use pencil_fft
  use iso_fortran_env, only : int64
  use halo_output, only: type_halo_catalog, type_halo_info
  implicit none
  save

  integer(4) t1,t2,tt1,tt2,ttt1,ttt2,t_rate,nhalo,ihalo
  integer(8) kg,jg,ig,ii,jj,imassbin
  real kr,kx,ky,kz,pow,r_filter,rm,masstemp
  real spin_u(3),spin_v(3),spin_r(3),spin_e(3,3),spin_vor(3)


  real phi(0:nf+1,0:nf+1,0:nf+1)[*]
  real phi_large(0:nf+1,0:nf+1,0:nf+1)[*]
  complex phi_k(nyquest+1,nf,npen)
  real delta2(ng,ng,ng)
  integer i,j,k,n_rsmall,n_ratio,nmassbin,itemp
  real t11,t22,t33,t12,t23,t31
  real tsmall(3,3),tlarge(3,3),torque(3,3)
  real spin(3,ng,ng,ng),spinmag(ng,ng,ng)
  type(type_halo_info) halo_info
  type(type_halo_catalog) halo_catalog
  real,allocatable :: spin_x(:,:),spin_q(:,:),spin_t(:,:),theta(:,:),imass(:),imass_info(:,:)
  real,allocatable :: corr_x(:,:,:),corr_q(:,:,:),corr_t(:,:,:),r_small(:),ratio_scale(:)
  integer,allocatable :: ind(:,:),isort_mass(:),i1(:),i2(:),ii1,ii2

#ifdef specify_potential
  character(*),parameter :: fn='/mnt/raid-cita/haoran/spin/cafproject/CUBEnu/output/universe38/image1/0.000_phiE_1.bin'
#endif

  ! nc: coarse grid per image per dim
  ! nf: fine grid per image per dim, ng=nf
  ! nyquest: Nyquest frequency
  call omp_set_num_threads(ncore)
  call geometry
  open(16,file='../main/z_checkpoint.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  cur_checkpoint=n_checkpoint ! read halos first
  if (head) then
    print*, ''
    print*, 'program LPT'
    print*, 'on',nn**3,' images and',ncore,' cores'
    print*, 'Resolution', ng*nn
    print*, 'at redshift=',z_checkpoint(cur_checkpoint)
    print*, 'Box size', box
    print*, 'output: ', opath
    print*, '-----------------------------------------'
  endif
  sync all
  if (head) print*,''
  if (head) print*,'Creating FFT plans'
  call system_clock(t1,t_rate)
  call create_penfft_plan
  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all
  ! open halo catalog
  print*,''
  print*,'read halo catalog'
  open(11,file=output_name('halo'),status='old',access='stream')
  open(12,file=output_name('halo_init_spin'),status='old',access='stream')
  open(13,file=output_name('halo_sort_index'),status='old',access='stream')
  read(11) halo_catalog
  nhalo=halo_catalog%nhalo
  allocate(isort_mass(nhalo),imass(nhalo))
  allocate(spin_x(3,nhalo),spin_q(3,nhalo),spin_t(3,nhalo),ind(3,nhalo),theta(3,nhalo))
  ! read halo q-pos & spin
  do ihalo=1,nhalo
    read(11) halo_info
    read(12) spin_q(:,ihalo),spin_t(:,ihalo),spin_x(:,ihalo),spin_u,spin_v,spin_r,spin_e(:,1:3)
    imass(ihalo)=halo_info%mass_odc
    ind(:,ihalo)=floor(halo_info%q_mean)+1
  enddo
  read(13) isort_mass ! read index by halo mass sort
  close(11)
  close(12)
  close(13)
  ! sort according to halo mass
  spin_x=spin_x(:,isort_mass)
  spin_q=spin_q(:,isort_mass)
  spin_t=spin_t(:,isort_mass)
  imass=imass(isort_mass)
  ind=ind(:,isort_mass)
  ! construct halo mass bins
  rm=2.0 ! mass bin ratio
  n_rsmall=50
  n_ratio=1
  nmassbin=ceiling(log(imass(1)/imass(nhalo))/log(rm))
  allocate(imass_info(4,nmassbin),i1(nmassbin),i2(nmassbin))
  allocate(corr_x(n_rsmall,n_ratio,nmassbin),corr_q(n_rsmall,n_ratio,nmassbin),corr_t(n_rsmall,n_ratio,nmassbin))
  allocate(r_small(n_rsmall),ratio_scale(n_ratio))
  print*,' N_halos_global =',halo_catalog%nhalo_tot
  print*,' N_halos_local  =',halo_catalog%nhalo
  print*,' Overdensity    =',halo_catalog%den_odc
  print*,' nmassbin =',nmassbin
  imassbin=nmassbin
  i2(nmassbin)=nhalo
  i1(1)=1
  itemp=1
  do ihalo=nhalo,1,-1
    if (imass(ihalo)<imass(nhalo)*rm**itemp .or. imassbin==1) then
      i1(imassbin)=ihalo ! still in the bin
    else
      imassbin=imassbin-1 ! assign to previous bin
      i2(imassbin)=ihalo
      itemp=itemp+1
    endif
  enddo
  ! construct scale bins
  do i=1,n_rsmall
    r_small(i)=0.1+0.2*(i-1)
  enddo
  do i=1,n_ratio
    ratio_scale(i)=1.1+0.1*(i-1)
  enddo

  ! open potential file
#ifdef specify_potential
  open(11,file=fn,access='stream')
#else
  cur_checkpoint=1;            open(11,file=output_name('phi1'),access='stream')
  !cur_checkpoint=n_checkpoint; open(11,file=output_name('phiE'),access='stream')
#endif
  read(11) r3
  close(11)
  call pencil_fft_forward ! complex phi stored in cxyz
  phi_k=cxyz ! store raw complex phi into phi_k


  cur_checkpoint=n_checkpoint
  open(11,file=output_name('lptcorr'),status='replace',access='stream')
  write(11),nmassbin,n_rsmall,n_ratio,imass_info(:,:),r_small(:),ratio_scale(:)
  do jj=1,n_ratio
    print*, jj,'/',n_ratio,' rs, ratio, t, q, x'
    do ii=1,n_rsmall
      call correlate_spin(r_small(ii),ratio_scale(jj)) ! loop over mass bins
    enddo
    print*,''
  enddo

  do imassbin=1,nmassbin
    ii1=i1(imassbin)
    ii2=i2(imassbin)
    imass_info(1,imassbin)=minval(imass(ii1:ii2))
    imass_info(2,imassbin)=maxval(imass(ii1:ii2))
    imass_info(3,imassbin)=sum(imass(ii1:ii2))/(ii2-ii1+1)
    imass_info(4,imassbin)=ii2-ii1+1
    write(11) corr_t(:,:,imassbin)
    write(11) corr_q(:,:,imassbin)
    write(11) corr_x(:,:,imassbin)
  enddo
  rewind(11)
  write(11) nmassbin,n_rsmall,n_ratio,imass_info(:,:)
  close(11)
call destroy_penfft_plan


contains

  subroutine correlate_spin(r_small,ratio_scale)
    implicit none
    save
    real r_small,ratio_scale
    call gaussian_fourier_filter(phi_k,r_small)
    call pencil_fft_backward
    phi(1:nf,1:nf,1:nf)=r3
    call buffer_1layer(phi)
    call gaussian_fourier_filter(phi_k,r_small*ratio_scale)
    call pencil_fft_backward
    phi_large(1:nf,1:nf,1:nf)=r3
    call buffer_1layer(phi_large)
    call spinfield
    do ihalo=1,nhalo
      theta(1,ihalo)=ccc(spin_t(:,ihalo),spin(:,ind(1,ihalo),ind(2,ihalo),ind(3,ihalo)))
      theta(2,ihalo)=ccc(spin_q(:,ihalo),spin(:,ind(1,ihalo),ind(2,ihalo),ind(3,ihalo)))
      theta(3,ihalo)=ccc(spin_x(:,ihalo),spin(:,ind(1,ihalo),ind(2,ihalo),ind(3,ihalo)))
    enddo
    do imassbin=1,nmassbin
      ii1=i1(imassbin)
      ii2=i2(imassbin)
      corr_t(ii,jj,imassbin)=sum(theta(1,ii1:ii2))/(ii2-ii1+1)
      corr_q(ii,jj,imassbin)=sum(theta(2,ii1:ii2))/(ii2-ii1+1)
      corr_x(ii,jj,imassbin)=sum(theta(3,ii1:ii2))/(ii2-ii1+1)
    enddo
    print*, r_small, ratio_scale, sum(theta,2)/nhalo
  endsubroutine

  real function ccc(vec1,vec2)
    implicit none
    real vec1(3),vec2(3)
    vec1=vec1/norm2(vec1)
    vec2=vec2/norm2(vec2)
    ccc=sum(vec1*vec2)
  endfunction

  subroutine buffer_1layer(phi)
    implicit none
    real phi(0:nf+1,0:nf+1,0:nf+1)[*]
    sync all
    phi(0,:,:)=phi(nf,:,:)[image1d(inx,icy,icz)]; sync all
    phi(nf+1,:,:)=phi(1,:,:)[image1d(ipx,icy,icz)]; sync all
    phi(:,0,:)=phi(:,nf,:)[image1d(icx,iny,icz)]; sync all
    phi(:,nf+1,:)=phi(:,1,:)[image1d(icx,ipy,icz)]; sync all
    phi(:,:,0)=phi(:,:,nf)[image1d(icx,icy,inz)]; sync all
    phi(:,:,nf+1)=phi(:,:,1)[image1d(icx,icy,ipz)]; sync all
  endsubroutine

  subroutine spinfield
    implicit none
    do k=1,ng
    do j=1,ng
    do i=1,ng
      tsmall(1,1)=phi(i+1,j,k)-2*phi(i,j,k)+phi(i-1,j,k)
      tsmall(2,2)=phi(i,j+1,k)-2*phi(i,j,k)+phi(i,j-1,k)
      tsmall(3,3)=phi(i,j,k+1)-2*phi(i,j,k)+phi(i,j,k-1)
      tsmall(1,2)=(phi(i+1,j+1,k)+phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))/4
      tsmall(2,3)=(phi(i,j+1,k+1)+phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))/4
      tsmall(3,1)=(phi(i+1,j,k+1)+phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))/4
      tsmall(2,1)=tsmall(1,2)
      tsmall(3,2)=tsmall(2,3)
      tsmall(1,3)=tsmall(3,1)

      tlarge(1,1)=phi_large(i+1,j,k)-2*phi_large(i,j,k)+phi_large(i-1,j,k)
      tlarge(2,2)=phi_large(i,j+1,k)-2*phi_large(i,j,k)+phi_large(i,j-1,k)
      tlarge(3,3)=phi_large(i,j,k+1)-2*phi_large(i,j,k)+phi_large(i,j,k-1)
      tlarge(1,2)=(phi_large(i+1,j+1,k)+phi_large(i-1,j-1,k)-phi_large(i+1,j-1,k)-phi_large(i-1,j+1,k))/4
      tlarge(2,3)=(phi_large(i,j+1,k+1)+phi_large(i,j-1,k-1)-phi_large(i,j+1,k-1)-phi_large(i,j-1,k+1))/4
      tlarge(3,1)=(phi_large(i+1,j,k+1)+phi_large(i-1,j,k-1)-phi_large(i+1,j,k-1)-phi_large(i-1,j,k+1))/4
      tlarge(2,1)=tlarge(1,2)
      tlarge(3,2)=tlarge(2,3)
      tlarge(1,3)=tlarge(3,1)

      torque=-matmul(tsmall,tlarge)
      spin(1,i,j,k)=-torque(2,3)+torque(3,2);
      spin(2,i,j,k)=-torque(3,1)+torque(1,3);
      spin(3,i,j,k)=-torque(1,2)+torque(2,1)
      spinmag(i,j,k)=norm2(spin(:,i,j,k))
    enddo
    enddo
    enddo

  endsubroutine

  subroutine gaussian_fourier_filter(phi_k,r_filter)
    ! apply Gaussian window function with radius r_filter to Fourier space phi_k
    ! returns Fourier space cxyz
    implicit none
    complex phi_k(nyquest+1,nf,npen)
    real r_filter
    cxyz=0
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
      kr=2*pi*kr/box
      pow=exp(-kr**2*r_filter**2/2)**0.25 ! apply E-mode window function
      cxyz(i,j,k)=phi_k(i,j,k)*pow
    enddo
    enddo
    enddo
    if (head) cxyz(1,1,1)=0 ! DC frequency
  endsubroutine

endprogram
