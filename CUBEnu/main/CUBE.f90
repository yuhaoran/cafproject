!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program CUBE
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
    call timestep
    call update_x
    call buffer_grid
    call buffer_x
    if (Extended_pp_force) then
      call ext_pp_force
    endif
    call particle_mesh
    call buffer_v
    if (checkpoint_step) then
      dt_old=0
      call update_x
      call checkpoint
      call projection
      call print_header(sim)
      if (final_step) exit
      call buffer_grid
      call buffer_x
      call buffer_v
      cur_checkpoint=cur_checkpoint+1
      checkpoint_step=.false.
      dt=0
    endif
  ENDDO

  if (head) close(77)

  call finalize

endprogram
