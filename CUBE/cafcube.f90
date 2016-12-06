!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Coarray CubeP3M       !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cafcube
use variables
use pm
implicit none
save

integer istep, hostnm,ierr
character(100) :: myhost

call initialize
call particle_initialization

sync all
print*,'buffer density'
call buffer_density
call buffer_x
call buffer_v

sync all

if (head) print*, '---------- starting main loop ----------'

DO istep=1,1000
  call timestep
  call update_particle
  call buffer_density
  call buffer_x
  call particle_mesh
  call buffer_v
  if (checkpoint_step) then
    dt_old=0
    call update_particle
    call checkpoint
    call projection
    call print_header(sim)
    if (final_step) exit
    call buffer_density
    call buffer_x
    call buffer_v
    cur_checkpoint=cur_checkpoint+1
    checkpoint_step=.false.
    dt=0
  endif
ENDDO

call finalize

endprogram
