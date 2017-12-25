!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cafcube
  use variables
  use neutrinos
  implicit none
  save

  if (this_image()==1) print*, 'Coarray CUBE on',nn**3,'  images'
  sync all
  call system('hostname')
  sync all

  call initialize
  call particle_initialization
  call buffer_np(rhoc)
  call buffer_np(rhoc_nu)
  call buffer_vc(vfield)
  call buffer_vc(vfield_nu)
  call redistribute_cdm()
  call redistribute_nu()
  call buffer_xp
  call buffer_xp_nu
  call buffer_vp
  call buffer_vp_nu
  if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')

  if (head) print*, '---------- starting main loop ----------'
  DO istep=1,istep_max
    call timestep
    call update_vp()
    call update_vp_nu()
    !call buffer_density
    ! the following to be collected
    call buffer_np(rhoc)
    call buffer_np(rhoc_nu)
    call buffer_vc(vfield)
    call buffer_vc(vfield_nu)
    call redistribute_cdm()
    call redistribute_nu()
    call buffer_xp
    call buffer_xp_nu

    !call buffer_x
    call particle_mesh ! to include nu
    !call buffer_v
    call buffer_vp
    call buffer_vp_nu
    if (checkpoint_step) then
      dt_old=0
      call update_vp()
      call update_vp_nu()
      call checkpoint ! to include nu
      call projection
      call print_header(sim)
      if (final_step) exit
      !call buffer_density
      !call buffer_x
      !call buffer_v
      call buffer_np(rhoc)
      call buffer_np(rhoc_nu)
      call buffer_vc(vfield)
      call buffer_vc(vfield_nu)
      call redistribute_cdm()
      call redistribute_nu()
      call buffer_xp
      call buffer_xp_nu
      call buffer_vp
      call buffer_vp_nu

      cur_checkpoint=cur_checkpoint+1
      checkpoint_step=.false.
      dt=0
    endif
  ENDDO

  if (head) close(77)

  call finalize

endprogram
