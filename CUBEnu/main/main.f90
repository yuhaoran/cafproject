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
  implicit none
  save

  call initialize

  call particle_initialization
  call buffer_grid
  call buffer_x
  call buffer_v
  if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')

  if (head) print*, '---------- starting main loop ----------'
  DO istep=1,istep_max
    call system_clock(ttt1,t_rate)
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
      if (checkpoint_step) call checkpoint
      call buffer_grid
      call buffer_x
      call buffer_v
      if (halofind_step) call halofind
      !call projection
      call print_header(sim)
      if (final_step) exit

      dt=0
    endif
    call system_clock(ttt2,t_rate)
    print*, 'total elapsed time =',real(ttt2-ttt1)/t_rate,'secs';
  ENDDO

  if (head) close(77)

  call finalize

endprogram
