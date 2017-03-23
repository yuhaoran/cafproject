!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Coarray CubeP3M       !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cafcube
use variables
use pm
implicit none
save

real abs_vsim(nf**3), std_vsim, max_vsim, kurt_vsim

if (this_image()==1) print*, 'Coarray CUBE on',nn**3,'  images'
sync all
call system('hostname')
sync all

call initialize
call particle_initialization

call buffer_density
call buffer_x
call buffer_v

#ifdef LINEAR_V
  if (head) open(77,file='vel_info.bin',access='stream',status='replace')
#endif

if (head) print*, '---------- starting main loop ----------'
DO istep=1,1000
  call timestep
  call update_particle

#ifdef LINEAR_V
    ! velocity analysis
    print*,'velocity analysis'
    print*,'  scale factor',a,a_mid

    abs_vsim= sqrt((v(1,:nplocal)*v_i2r(1))**2 &
                  +(v(2,:nplocal)*v_i2r(2))**2 &
                  +(v(3,:nplocal)*v_i2r(3))**2)
    max_vsim=maxval(abs_vsim)
    std_vsim=sqrt(sum(abs_vsim**2/nplocal*1d0))
    kurt_vsim=sum(abs_vsim**4*1d0/nplocal)/std_vsim**4
    print*,'  vmax',max_vsim
    print*,'   std',std_vsim
    print*,'  kurt',kurt_vsim

    write(77) a-da,std_vsim,max_vsim,kurt_vsim
#endif

  sync all
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

#ifdef LINEAR_V
  if (head) close(77)
#endif

call finalize

endprogram
