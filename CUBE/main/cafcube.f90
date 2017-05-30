!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cafcube
  use variables
  implicit none
  save

  real(8) abs_vsim, std_vsim, max_vsim, kurt_vsim

  if (this_image()==1) print*, 'Coarray CUBE on',nn**3,'  images'
  sync all
  call system('hostname')
  sync all

  call initialize
  call particle_initialization
  if (head) print*,v(:,1),v(:,nplocal)
  call buffer_density
  if (head) print*,v(:,1),v(:,nplocal)
  call buffer_x
  if (head) print*,v(:,1),v(:,nplocal)
  call buffer_v
  if (head) print*,v(:,1),v(:,nplocal)
  if (head) print*,' buffer before main loop done'

  if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')

  if (head) print*, '---------- starting main loop ----------'
  DO istep=1,istep_max
    call timestep
    call update_particle

      ! velocity analysis
      if (head) print*,'velocity analysis'
      if (head) print*,'  scale factor',a,a_mid
      max_vsim=0; std_vsim=0; kurt_vsim=0
      if (head) print*, 'nplocal',nplocal
      do ip=1,nplocal
        vreal=tan(pi*real(v(:,ip))/real(nvbin-1))/(sqrt(pi/2)/sigma_vi_old)
        abs_vsim=sqrt(sum(vreal**2))
        max_vsim=max(max_vsim,abs_vsim)
        std_vsim=std_vsim+abs_vsim**2!         sqrt(sum(abs_vsim**2/nplocal*1d0))
        kurt_vsim=kurt_vsim+abs_vsim**4!       sum(abs_vsim**4*1d0/nplocal)/std_vsim**4
        !print*,abs_vsim,max_vsim,std_vsim,kurt_vsim
      enddo
      std_vsim=sqrt(std_vsim/nplocal)
      kurt_vsim=kurt_vsim/nplocal/std_vsim**4
      print*,'  vmax',max_vsim
      print*,'   std',std_vsim
      print*,' std_L',sigma_vi_old*sqrt(3.)
      !print*,'  kurt',kurt_vsim

      write(77) a-da,real(std_vsim,4),real(max_vsim,4),real(kurt_vsim,4)

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

  if (head) close(77)

  call finalize

endprogram
