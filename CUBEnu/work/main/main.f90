!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    integer,parameter :: nhmax=3000
    integer(2) halo_id,hid(nf**3),ih
    integer(4) nhalo,nhalo_tot,nwrite
    real den_odc,hpos(3,nhmax),spin_q(3,nhmax),x_mean(3,nhmax),v_mean(3,nhmax),hp(nhmax)
    real hcorr(nhmax,istep_max),scale_fac(istep_max)
    real spin_x(3,nhmax),dx(3),dv(3)
    type(type_halo_info) halo_info
# endif
  save

  call initialize
  call particle_initialization
# ifdef halo_spin_correlation
    open(19,file=output_dir()//'0.000_hid'//output_suffix(),access='stream',status='old')
    read(19) hid
    close(19)
    open(19,file=output_dir()//'0.000_halo'//output_suffix(),access='stream',status='old')
    read(19) nhalo_tot,nhalo,den_odc
    do ih=1,nhalo
      read(19) halo_info
      hpos(:,ih)=halo_info%x_mean
    enddo
    close(19)
    open(19,file=output_dir()//'0.000_halo_init_spin'//output_suffix(),access='stream',status='old')
    read(19) spin_q(:,:nhalo)
    close(19)
    nwrite=0
# endif
  call buffer_grid
  call buffer_x
  call buffer_v
  cur_checkpoint=cur_checkpoint+1
  cur_halofind=cur_checkpoint+1
  if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')




  if (head) print*, '---------- starting main loop ----------'
  DO istep=1,istep_max
    call system_clock(ttt1,t_rate)

!#   ifdef halo_spin_correlation
      dt_old=0
      call update_x
      call buffer_grid
      call buffer_x
      call buffer_v
#    ifdef halo_spin_correlation
      call spin_analysis
#    endif
      dt=0
!#   endif
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
    call particle_mesh
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
  open(19,file=output_dir()//'halo_spin_corr'//output_suffix(),access='stream',status='replace')
  write(19) nhalo,nwrite
  write(19) scale_fac(:nwrite),hcorr(:nhalo,:nwrite)
  close(19)
#endif

  if (head) close(77)

  call finalize


contains

#ifdef halo_spin_correlation
  subroutine spin_analysis
#   ifdef HALOFIND
      !stop "disable HALOFIND when using halo_spin_correlation"
#   endif
#   ifndef PID
      stop "enable PID when using halo_spin_correlation"
#   endif
    !! halo particle cross correlation
    spin_x=0; hp=0; x_mean=0; v_mean=0
    nwrite=nwrite+1
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
          dx=xq*ncell-x_mean(:,halo_id)
          dx=modulo(dx+nf/2,real(nf))-nf/2

          vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          dv=vreal!-v_mean(:,halo_id)

          spin_x(1,halo_id)=spin_x(1,halo_id)+dx(2)*dv(3)-dx(3)*dv(2)
          spin_x(2,halo_id)=spin_x(2,halo_id)+dx(3)*dv(1)-dx(1)*dv(3)
          spin_x(3,halo_id)=spin_x(3,halo_id)+dx(1)*dv(2)-dx(2)*dv(1)
        enddo
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo

    scale_fac(nwrite)=a
    do ih=1,nhalo
      hcorr(ih,nwrite)=sum(spin_q(:,ih)*spin_x(:,ih))/sqrt(sum(spin_q(:,ih)**2))/sqrt(sum(spin_x(:,ih)**2))
    enddo
    print*,'=========== halo 4 ==========='
    print*,x_mean(:,4)
    print*,v_mean(:,4)
    print*,spin_q(:,4)
    print*,spin_x(:,4)
    print*,hcorr(4,nwrite)
    print*,'at scale factor',a
    print*,'=============================='
  endsubroutine
#endif


endprogram
