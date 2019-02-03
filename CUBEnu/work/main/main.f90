!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define record_Fc
program main
  use omp_lib
  use variables
  use neutrinos
  use buffer_grid_subroutines
  use buffer_particle_subroutines
  use update_particle
  use pp_force
# ifdef halo_spin_correlation
    use halo_output, only: type_halo_info
# endif
  implicit none

# ifdef halo_spin_correlation
    integer(2) halo_id,hid(nf**3),ih
    integer(4) nwrite,i0,j0,ipos(3)
    integer(4),allocatable :: isort_mass(:)
    real,allocatable :: hpos(:,:),spin_q(:,:,:),x_mean(:,:),v_mean(:,:),hp(:),scale_fac(:),spin_x(:,:),halonp(:)
    real,allocatable :: dspinx(:,:),corrdj(:),ratiodj(:),jtot(:,:,:)
    real,allocatable :: inertia(:,:,:),tide(:,:,:),tidu(:,:,:),torque(:,:,:),eIT(:,:,:)
    real,allocatable :: haloloc(:,:),haloforce(:,:)
    real,allocatable :: force_cdm(:,:),force_neu(:,:)
    real,allocatable :: phi(:,:,:)[:],phu(:,:,:)[:]
    real,allocatable :: hinfo(:,:,:)
    real ztfnu(10),atfnu(10)
    real den_odc,dx(3),dv(3)
    type(type_halo_info) halo_info
# endif
  save

  call initialize
  call particle_initialization
# ifdef halo_spin_correlation
    open(21,file='tf_nu.txt',status='old')
    read(21,fmt='(f8.4)') ztfnu
    close(21)
    atfnu=1.0/(1+ztfnu)

    open(19,file=output_dir()//'0.000_hid'//output_suffix(),access='stream',status='old')
    read(19) hid
    close(19)

    open(19,file=output_dir()//'0.000_halo'//output_suffix(),access='stream',status='old')
    read(19) nhalo_tot,nhalo,den_odc
    ! allocate arrays for halos
    allocate(isort_mass(nhalo),halonp(nhalo),hpos(3,nhalo),spin_q(3,9,nhalo),x_mean(3,nhalo),v_mean(3,nhalo),hp(nhalo))
    allocate(hinfo(15,nhalo,istep_max),scale_fac(istep_max),spin_x(3,nhalo),dspinx(3,nhalo))
    allocate(corrdj(nhalo),ratiodj(nhalo))
    allocate(jtot(3,2,nhalo))
    allocate(inertia(3,3,nhalo),torque(3,3,nhalo),eIT(3,2,nhalo))
    allocate(force_cdm(3,nhalo),force_neu(3,nhalo))
    allocate(haloloc(3,nhalo),haloforce(3,nhalo))
    ! read halo_info and get their positions and masses
    do ih=1,nhalo
      read(19) halo_info
      hpos(:,ih)=halo_info%x_mean
      halonp(ih)=halo_info%mass_odc
    enddo
    close(19)
    open(19,file=output_dir()//'0.000_halo_init_spin'//output_suffix(),access='stream',status='old')
    read(19) spin_q
    close(19)
    do ih=1,nhalo
      jtot(:,1,ih)=0
      jtot(:,2,ih)=0
    enddo
    nwrite=0
# endif

  call buffer_grid
  call buffer_x
  call buffer_v
  cur_checkpoint=cur_checkpoint+1
  cur_halofind=cur_checkpoint+1
  if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')
#ifdef record_Fc
  !if (head) open(65,file=output_dir()//'zFc'//output_suffix(),form='FORMATTED',status='replace')
#endif


  if (head) print*, '---------- starting main loop ----------'
  DO istep=1,istep_max
    call system_clock(ttt1,t_rate)

#   ifdef halo_spin_correlation
      dt_old=0
      call update_x
      call buffer_grid
      call buffer_x
      call buffer_v
#    ifdef halo_spin_correlation
      call spin_analysis
#    endif
      dt=0
#   endif
    call timestep
    call update_x
#   ifdef FORCETEST
      sim%nplocal=2
      sim%mass_p_cdm=1
      rhoc=0;vfield=0;xp=0;vp=0
      rhoc(1,1,1,1,1,1)=2
      xp(:,1:2)=16383
      xp(1,2)=-16385 ! offset the second particle to +x by 1 fine cell.
#   endif
    call buffer_grid
    call buffer_x
    if (Extended_pp_force) then
      call ext_pp_force
    endif
    call particle_mesh !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call buffer_v
#   ifdef FORCETEST
      stop
#   endif
    if (checkpoint_step .or. halofind_step) then
      dt_old=0
      call update_x
      if (checkpoint_step) then
        call checkpoint
        cur_checkpoint=cur_checkpoint+1
      endif
      call buffer_grid
      call buffer_x
      call buffer_v
#     ifdef halo_spin_correlation
        call spin_analysis
#     endif
      if (halofind_step) then
        call halofind
        cur_halofind=cur_halofind+1
      endif
      call print_header(sim)
      if (final_step) exit
      dt=0
    endif
    call system_clock(ttt2,t_rate)
    print*, 'total elapsed time =',real(ttt2-ttt1)/t_rate,'secs';
  ENDDO

#ifdef halo_spin_correlation
  isort_mass(:nhalo)=[(i0,i0=1,nhalo)]
  call indexedsort(nhalo,halonp,isort_mass)
  do i0=1,istep_max
  do j0=1,15
    hinfo(j0,:,i0)=hinfo(j0,isort_mass(nhalo:1:-1),i0)
  enddo
  enddo
  open(19,file=output_dir()//'halo_spin_corr'//output_suffix(),access='stream',status='replace')
  write(19) nhalo,nwrite
  write(19) scale_fac(:nwrite),halonp(nhalo:1:-1)
  write(19) hinfo(:,:,:nwrite)

  write(19) spin_x(:,isort_mass(nhalo:1:-1))
  write(19) jtot(:,2,isort_mass(nhalo:1:-1))
  write(19) spin_q(:,6,isort_mass(nhalo:1:-1))
  write(19) spin_q(:,7,isort_mass(nhalo:1:-1))
  write(19) spin_q(:,8,isort_mass(nhalo:1:-1))
  write(19) spin_q(:,9,isort_mass(nhalo:1:-1))

  close(19)
#endif

  if (head) close(77)
#ifdef record_Fc
  !if (head) close(65)
#endif
  call finalize


contains













#ifdef halo_spin_correlation
  subroutine spin_analysis
    integer iztf(1)
    character(:),allocatable :: prefix_nu,prefix_dm
    character(20) :: str_zdm,str_znu,str_i

    print*,''
    print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    print*,'spin_analysis'
    a_mid=max(a_mid,1e-6)
    dt=max(dt,1e-3)
    print*,'  use scale factor a_mid =',a_mid

    allocate(tide(3,3,nhalo),tidu(3,3,nhalo))
    allocate(phi(0:nf+1,0:nf+1,0:nf+1)[*])
    allocate(phu(0:nf+1,0:nf+1,0:nf+1)[*])
    iztf=minloc(abs(atfnu-a_mid))
    write(str_i,'(i6)') image
    write(str_zdm,'(f7.3)') z_checkpoint(1)
    write(str_znu,'(f7.3)') ztfnu(iztf)
    !write(str_znu,'(f7.3)') 5.0

!  open(21,file='../../output/universe2/image1/0.000_phiE_1.bin',status='old',access='stream')

!    prefix_dm=opath//'image'//trim(adjustl(str_i))//'/0.000_phiE' ! E-mode
    prefix_dm=opath//'image'//trim(adjustl(str_i))//'/'//trim(adjustl(str_zdm))//'_phi1' ! filtered
    prefix_nu=opath//'image'//trim(adjustl(str_i))//'/5.000_tf'//trim(adjustl(str_znu))//'_'
    print*, 'phi_dm ',prefix_dm//output_suffix()
    print*, 'phi_nu ',prefix_nu//'phi1_nu'//output_suffix()

    open(21,file=prefix_dm//output_suffix(),status='old',access='stream')
    read(21) phi(1:nf,1:nf,1:nf)
    close(21)
    open(21,file=prefix_nu//'phi1_nu'//output_suffix(),status='old',access='stream')
    read(21) phu(1:nf,1:nf,1:nf)
    close(21)
    phi=phi*Dgrow(a_mid)/Dgrow(1/(1+z_checkpoint(1))) ! for filtered phi1
!    phi=phi*Dgrow(a_mid) ! for E-mode
    phu=phu*Dgrow(a_mid)/Dgrow(1./6.)

    if (head) print*, '  buffer phi'
    phi(0,:,:)=phi(nf,:,:)[image1d(inx,icy,icz)]
    phi(nf+1,:,:)=phi(1,:,:)[image1d(ipx,icy,icz)]
    sync all
    phi(:,0,:)=phi(:,nf,:)[image1d(icx,iny,icz)]
    phi(:,nf+1,:)=phi(:,1,:)[image1d(icx,ipy,icz)]
    sync all
    phi(:,:,0)=phi(:,:,nf)[image1d(icx,icy,inz)]
    phi(:,:,nf+1)=phi(:,:,1)[image1d(icx,icy,ipz)]
    sync all

    if (head) print*, '  buffer phu'
    phu(0,:,:)=phu(nf,:,:)[image1d(inx,icy,icz)]
    phu(nf+1,:,:)=phu(1,:,:)[image1d(ipx,icy,icz)]
    sync all
    phu(:,0,:)=phu(:,nf,:)[image1d(icx,iny,icz)]
    phu(:,nf+1,:)=phu(:,1,:)[image1d(icx,ipy,icz)]
    sync all
    phu(:,:,0)=phu(:,:,nf)[image1d(icx,icy,inz)]
    phu(:,:,nf+1)=phu(:,:,1)[image1d(icx,icy,ipz)]
    sync all

    dspinx=spin_x;
    spin_x=0; hp=0; x_mean=0; v_mean=0
    inertia=0; tide=0; tidu=0; force_cdm=0; force_neu=0
    nwrite=nwrite+1
    ! loop over all particles, link to halo_id.
    ! first loop, get x_mean & v_mean for each halo_id.
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np
          ip=nzero+l
          halo_id=hid(pid(ip))
          if (halo_id==0) cycle
          hp(halo_id)=hp(halo_id)+1

          xq=nt*([itx,ity,itz]-1d0)+([i,j,k]-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          dx=xq*ncell-hpos(:,halo_id)
          dx=modulo(dx+nf/2,real(nf))-nf/2
          if (maxval(abs(dx))>64) stop "warning 1"
          x_mean(:,halo_id)=x_mean(:,halo_id)+dx ! x_mean is relative to hpos at z=0

          vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          v_mean(:,halo_id)=v_mean(:,halo_id)+vreal
        enddo
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo

    do ih=1,nhalo
      x_mean(:,ih)=x_mean(:,ih)/hp(ih)+hpos(:,ih)
      v_mean(:,ih)=v_mean(:,ih)/hp(ih)
    enddo
    x_mean(:,:nhalo)=modulo(x_mean(:,:nhalo),real(nf)) ! x_mean is the global coordinate

    ! second loop, get spin_x, inertia, tide, tidu, torque.
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np
          ip=nzero+l
          if (pid(ip) == 0) then
            print*, 'warning pid=0',ip
            cycle
          endif 
          halo_id=hid(pid(ip))
          if (halo_id==0) cycle ! particle does not belong to a halo

          xq=nt*([itx,ity,itz]-1d0)+([i,j,k]-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          dx=xq*ncell-x_mean(:,halo_id)
          dx=modulo(dx+nf/2,real(nf))-nf/2

          vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          dv=vreal!-v_mean(:,halo_id)

          spin_x(1,halo_id)=spin_x(1,halo_id)+dx(2)*dv(3)-dx(3)*dv(2)
          spin_x(2,halo_id)=spin_x(2,halo_id)+dx(3)*dv(1)-dx(1)*dv(3)
          spin_x(3,halo_id)=spin_x(3,halo_id)+dx(1)*dv(2)-dx(2)*dv(1)

          inertia(1,1,halo_id)=inertia(1,1,halo_id)+dx(1)**2
          inertia(2,2,halo_id)=inertia(2,2,halo_id)+dx(2)**2
          inertia(3,3,halo_id)=inertia(3,3,halo_id)+dx(3)**2
          inertia(1,2,halo_id)=inertia(1,2,halo_id)+dx(1)*dx(2)
          inertia(2,3,halo_id)=inertia(2,3,halo_id)+dx(2)*dx(3)
          inertia(3,1,halo_id)=inertia(3,1,halo_id)+dx(3)*dx(1)

          !ipos=ceiling(hpos(:,halo_id)) ! center part of halo
          ipos=ceiling(xq*ncell) ! actual shape
          tide(1,1,halo_id)=tide(1,1,halo_id)+&
            phi(ipos(1)+1,ipos(2),ipos(3))-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1)-1,ipos(2),ipos(3))
          tide(2,2,halo_id)=tide(2,2,halo_id)+&
            phi(ipos(1),ipos(2)+1,ipos(3))-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1),ipos(2)-1,ipos(3))
          tide(3,3,halo_id)=tide(3,3,halo_id)+&
            phi(ipos(1),ipos(2),ipos(3)+1)-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1),ipos(2),ipos(3)-1)
          tide(1,2,halo_id)=tide(1,2,halo_id)+(phi(ipos(1)+1,ipos(2)+1,ipos(3))+phi(ipos(1)-1,ipos(2)-1,ipos(3))&
            -phi(ipos(1)+1,ipos(2)-1,ipos(3))-phi(ipos(1)-1,ipos(2)+1,ipos(3)))/4
          tide(2,3,halo_id)=tide(2,3,halo_id)+(phi(ipos(1),ipos(2)+1,ipos(3)+1)+phi(ipos(1),ipos(2)-1,ipos(3)-1)&
            -phi(ipos(1),ipos(2)+1,ipos(3)-1)-phi(ipos(1),ipos(2)-1,ipos(3)+1))/4
          tide(3,1,halo_id)=tide(3,1,halo_id)+(phi(ipos(1)+1,ipos(2),ipos(3)+1)+phi(ipos(1)-1,ipos(2),ipos(3)-1)&
            -phi(ipos(1)+1,ipos(2),ipos(3)-1)-phi(ipos(1)-1,ipos(2),ipos(3)+1))/4
          haloforce(1,halo_id)=haloforce(1,halo_id)-&
            (phi(ipos(1)+1,ipos(2),ipos(3))-phi(ipos(1)-1,ipos(2),ipos(3)))/2
          haloforce(2,halo_id)=haloforce(2,halo_id)-&
            (phi(ipos(1),ipos(2)+1,ipos(3))-phi(ipos(1),ipos(2)-1,ipos(3)))/2
          haloforce(3,halo_id)=haloforce(3,halo_id)-&
            (phi(ipos(1),ipos(2),ipos(3)+1)-phi(ipos(1),ipos(2),ipos(3)-1))/2

          !ipos=ceiling(hpos(:,halo_id)) ! center part of halo
          ipos=ceiling(xq*ncell) ! actual shape
          tidu(1,1,halo_id)=tidu(1,1,halo_id)+&
            phu(ipos(1)+1,ipos(2),ipos(3))-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1)-1,ipos(2),ipos(3))
          tidu(2,2,halo_id)=tidu(2,2,halo_id)+&
            phu(ipos(1),ipos(2)+1,ipos(3))-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1),ipos(2)-1,ipos(3))
          tidu(3,3,halo_id)=tidu(3,3,halo_id)+&
            phu(ipos(1),ipos(2),ipos(3)+1)-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1),ipos(2),ipos(3)-1)
          tidu(1,2,halo_id)=tidu(1,2,halo_id)+(phu(ipos(1)+1,ipos(2)+1,ipos(3))+phu(ipos(1)-1,ipos(2)-1,ipos(3))&
            -phu(ipos(1)+1,ipos(2)-1,ipos(3))-phu(ipos(1)-1,ipos(2)+1,ipos(3)))/4
          tidu(2,3,halo_id)=tidu(2,3,halo_id)+(phu(ipos(1),ipos(2)+1,ipos(3)+1)+phu(ipos(1),ipos(2)-1,ipos(3)-1)&
            -phu(ipos(1),ipos(2)+1,ipos(3)-1)-phu(ipos(1),ipos(2)-1,ipos(3)+1))/4
          tidu(3,1,halo_id)=tidu(3,1,halo_id)+(phu(ipos(1)+1,ipos(2),ipos(3)+1)+phu(ipos(1)-1,ipos(2),ipos(3)-1)&
            -phu(ipos(1)+1,ipos(2),ipos(3)-1)-phu(ipos(1)-1,ipos(2),ipos(3)+1))/4
        enddo
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    inertia(2,1,:)=inertia(1,2,:)
    inertia(3,2,:)=inertia(2,3,:)
    inertia(1,3,:)=inertia(3,1,:)
    tide(2,1,:)=tide(1,2,:)
    tide(3,2,:)=tide(2,3,:)
    tide(1,3,:)=tide(3,1,:)
    tidu(2,1,:)=tidu(1,2,:)
    tidu(3,2,:)=tidu(2,3,:)
    tidu(1,3,:)=tidu(3,1,:)
    dspinx(:,:ih)=spin_x(:,:ih)-dspinx(:,:ih)
    do ih=1,nhalo
      tide(:,:,ih)=tide(:,:,ih)/hp(ih)
      tidu(:,:,ih)=tidu(:,:,ih)/hp(ih)
      haloforce(:,ih)=haloforce(:,ih)/hp(ih)
      torque(:,:,ih)=matmul(inertia(:,:,ih),tide(:,:,ih))
      eIT(1,1,ih)=-torque(2,3,ih)+torque(3,2,ih)
      eIT(2,1,ih)=-torque(3,1,ih)+torque(1,3,ih)
      eIT(3,1,ih)=-torque(1,2,ih)+torque(2,1,ih)
      torque(:,:,ih)=matmul(inertia(:,:,ih),tidu(:,:,ih))
      eIT(1,2,ih)=-torque(2,3,ih)+torque(3,2,ih)
      eIT(2,2,ih)=-torque(3,1,ih)+torque(1,3,ih)
      eIT(3,2,ih)=-torque(1,2,ih)+torque(2,1,ih)

      jtot(:,:,ih)=jtot(:,:,ih)+eIT(:,:,ih)*a_mid*dt/6/pi
      corrdj(ih)=sum(eIT(:,1,ih)*dspinx(:,ih))/sqrt(sum(eIT(:,1,ih)**2))/sqrt(sum(dspinx(:,ih)**2))
      ratiodj(ih)=sqrt( sum((eIT(:,1,ih)*a_mid*dt/6/pi)**2) / sum(dspinx(:,ih)**2))
    enddo

    if (.false.) then
      ih=1
      print*,'halo1, np',halonp(ih),hp(ih)
      print*,'x_mean',x_mean(:,ih)
      print*,'force',haloforce(:,ih)
      print*,'inertia'
      print*,inertia(:,:,ih)
      print*,'tide'
      print*,tide(:,:,ih)
      print*,'tide from neutrinos'
      print*,tidu(:,:,ih)

      print*,'j1'
      print*,spin_q(:,1,ih)
      print*,'eIT * a_mid * dt / ( 6 pi) ='
      print*,eIT(:,1,ih)*a_mid*dt/6/pi
      print*,'corr j1-eIT'
      print*,sum(spin_q(:,1,ih)*eIT(:,1,ih))/sqrt(sum(spin_q(:,1,ih)**2))/sqrt(sum(eIT(:,1,ih)**2))
      print*,'SPIN_x ='
      print*,spin_x(:,ih)
      print*,'eIT - delta_j correlation'
      print*,sum(corrdj(:nhalo))/nhalo
      print*,'ratio of amplitude'
      print*,sum(ratiodj(:nhalo))/nhalo
      print*,'corr CDM and neutrino eIT'
      print*,sum(eIT(:,1,ih)*eIT(:,2,ih))/sqrt(sum(eIT(:,1,ih)**2))/sqrt(sum(eIT(:,2,ih)**2))
    endif

    scale_fac(nwrite)=a
    do ih=1,nhalo
      ! correlation := sum(A*B)/sqrt(sum(A**2)*sum(B**2))
      hinfo(1,ih,nwrite)=sum(spin_q(:,1,ih)*spin_x(:,ih))/sqrt(sum(spin_q(:,1,ih)**2)*sum(spin_x(:,ih)**2))
      hinfo(2,ih,nwrite)=sum(spin_q(:,2,ih)*spin_x(:,ih))/sqrt(sum(spin_q(:,2,ih)**2)*sum(spin_x(:,ih)**2))
      hinfo(3,ih,nwrite)=sum(spin_q(:,3,ih)*spin_x(:,ih))/sqrt(sum(spin_q(:,3,ih)**2)*sum(spin_x(:,ih)**2))
      hinfo(4,ih,nwrite)=sum(jtot(:,1,ih)*jtot(:,2,ih))/sqrt(sum(jtot(:,1,ih)**2)*sum(jtot(:,2,ih)**2))

      hinfo(5,ih,nwrite)=sum(jtot(:,1,ih)*spin_x(:,ih))/sqrt(sum(jtot(:,1,ih)**2)*sum(spin_x(:,ih)**2))
      hinfo(6,ih,nwrite)=sum(jtot(:,2,ih)*spin_x(:,ih))/sqrt(sum(jtot(:,2,ih)**2)*sum(spin_x(:,ih)**2))
      hinfo(7,ih,nwrite)=sum(jtot(:,2,ih)*spin_q(:,5,ih))/sqrt(sum(jtot(:,2,ih)**2)*sum(spin_q(:,5,ih)**2))
      hinfo(8,ih,nwrite)=sum(jtot(:,2,ih)*spin_q(:,6,ih))/sqrt(sum(jtot(:,2,ih)**2)*sum(spin_q(:,6,ih)**2))

      !hinfo(9,ih,nwrite)=sum(spin_x(:,ih)*spin_q(:,6,ih))/sqrt(sum(spin_x(:,ih)**2)*sum(spin_q(:,6,ih)**2))
      !hinfo(10,ih,nwrite)=sqrt(sum(spin_x(:,ih)**2))
      hinfo(9,ih,nwrite)=sum(eIT(:,1,ih)*spin_q(:,2,ih))/sqrt(sum(eIT(:,1,ih)**2)*sum(spin_q(:,2,ih)**2))
      hinfo(10,ih,nwrite)=sum(eIT(:,2,ih)*spin_q(:,2,ih))/sqrt(sum(eIT(:,2,ih)**2)*sum(spin_q(:,2,ih)**2))
      hinfo(11,ih,nwrite)=sum(eIT(:,2,ih)*spin_q(:,5,ih))/sqrt(sum(eIT(:,2,ih)**2)*sum(spin_q(:,5,ih)**2))
      !hinfo(11,ih,nwrite)=sqrt(sum(jtot(:,1,ih)**2))
      hinfo(12,ih,nwrite)=sum(eIT(:,2,ih)*spin_q(:,6,ih))/sqrt(sum(eIT(:,2,ih)**2)*sum(spin_q(:,6,ih)**2))
      !hinfo(12,ih,nwrite)=sqrt(sum(jtot(:,2,ih)**2))

      hinfo(13,ih,nwrite)=sum(jtot(:,2,ih)*spin_q(:,7,ih))/sqrt(sum(jtot(:,2,ih)**2)*sum(spin_q(:,7,ih)**2))
      hinfo(14,ih,nwrite)=sum(jtot(:,2,ih)*spin_q(:,8,ih))/sqrt(sum(jtot(:,2,ih)**2)*sum(spin_q(:,8,ih)**2))
      hinfo(15,ih,nwrite)=sum(jtot(:,2,ih)*spin_q(:,9,ih))/sqrt(sum(jtot(:,2,ih)**2)*sum(spin_q(:,9,ih)**2))

    enddo
    ih=1
    print*,'=========== halo average ==========='
    print*,'qx',sum(hinfo(1,:nhalo,nwrite))/nhalo
    print*,'tx',sum(hinfo(2,:nhalo,nwrite))/nhalo
    print*,'xx',sum(hinfo(3,:nhalo,nwrite))/nhalo
    print*,'cu',sum(hinfo(4,:nhalo,nwrite))/nhalo
    print*,sum(hinfo(5,:nhalo,nwrite))/nhalo
    print*,sum(hinfo(6,:nhalo,nwrite))/nhalo
    print*,'<Jv . eIqTv> =',sum(hinfo(7,:nhalo,nwrite))/nhalo
    print*,'<Jv . eTcTv> =',sum(hinfo(8,:nhalo,nwrite))/nhalo
    print*,sum(hinfo(9,:nhalo,nwrite))/nhalo
    print*,sum(hinfo(10,:nhalo,nwrite))/nhalo
    print*,'<eI(z)Tv . eIqTv> =',sum(hinfo(11,:nhalo,nwrite))/nhalo
    print*,'<eI(z)Tv . eTcTv> =',sum(hinfo(12,:nhalo,nwrite))/nhalo
    print*,'<J_nu . eTc1Tv> =',sum(hinfo(13,:nhalo,nwrite))/nhalo
    print*,'<J_nu . eTc2Tv> =',sum(hinfo(14,:nhalo,nwrite))/nhalo
    print*,'<J_nu . eTc3Tv> =',sum(hinfo(15,:nhalo,nwrite))/nhalo

    print*,'at scale factor',a
    print*,'=============================='
    deallocate(phi,phu,tide,tidu)
  endsubroutine
#endif

  function Dgrow(scale_factor)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real scale_factor
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/scale_factor**3+(1-om-ol)/scale_factor**2+ol
    oma=om/(scale_factor**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=scale_factor*ga/g
  end function Dgrow

endprogram
