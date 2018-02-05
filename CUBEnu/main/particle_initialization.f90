subroutine particle_initialization
  use variables
  use neutrinos

  implicit none
  save

  character(100) fn10,fn11,fn12,fn13,fn14,fn15,fn21,fn22,fn23,fn24,fn25

  if (head) then
    print*, 'particle_initialization'
    call system_clock(t1,t_rate)
  endif
  fn10=ic_name('info')
  fn11=ic_name('xp')
  fn12=ic_name('vp')
  fn13=ic_name('np')
  fn14=ic_name('vc')
  fn15=ic_name('id')
  fn21=ic_name_nu('xp_nu')
  fn22=ic_name_nu('vp_nu')
  fn23=ic_name_nu('np_nu')
  fn24=ic_name_nu('vc_nu')
  fn25=ic_name_nu('id_nu')

  open(10,file=fn10,status='old',access='stream')
  read(10) sim
  close(10)
  nplocal=sim%nplocal
  nplocal_nu=sim%nplocal_nu
  sigma_vi=sim%sigma_vi
  sigma_vi_nu=sim%sigma_vi_nu
  dt_vmax=sim%dt_vmax
  dt_vmax_nu=sim%dt_vmax_nu

  if (sim%izipx/=izipx .or. sim%izipv/=izipv .or. sim%izipx_nu/=izipx_nu .or. sim%izipv_nu/=izipv_nu) then
    print*, '  zip format incompatable'
    close(12)
    stop
  endif

! may create separate file names
  !$omp parallelsections default(shared)
  !$omp section
    open(11,file=fn11,status='old',access='stream'); read(11) xp(:,:nplocal); close(11)
  !$omp section
    open(12,file=fn12,status='old',access='stream'); read(12) vp(:,:nplocal); close(12)
  !$omp section
    open(13,file=fn13,status='old',access='stream'); read(13) rhoc(1:nt,1:nt,1:nt,:,:,:); close(13)
  !$omp section
    open(14,file=fn14,status='old',access='stream'); read(14) vfield(:,1:nt,1:nt,1:nt,:,:,:); close(14)
#ifdef PID
  !$omp section
    open(15,file=fn15,status='old',access='stream'); read(15) pid(:nplocal); close(15)
#endif

  !$omp section
    open(21,file=fn21,status='old',access='stream'); read(21) xp_nu(:,:nplocal_nu); close(21)
  !$omp section
    open(22,file=fn22,status='old',access='stream'); read(22) vp_nu(:,:nplocal_nu); close(22)
  !$omp section
    open(23,file=fn23,status='old',access='stream'); read(23) rhoc_nu(1:nt,1:nt,1:nt,:,:,:); close(23)
  !$omp section
    open(24,file=fn24,status='old',access='stream'); read(24) vfield_nu(:,1:nt,1:nt,1:nt,:,:,:); close(24)
#ifdef EID
  !$omp section
    open(25,file=fn25,status='old',access='stream'); read(25) pid_nu(:nplocal_nu); close(25)
#endif
  !$omp endparallelsections

  a=a_i
  npglobal=0
  npglobal_nu=0
  sync all

  do i=1,nn**3
    npglobal=npglobal+nplocal[i]
    npglobal_nu=npglobal_nu+nplocal_nu[i]
  enddo
  mass_p=real((nf*nn)**3*f_cdm)/npglobal ! if a CDM-only simulation
  mass_p_cdm=real((nf*nn)**3*f_cdm)/npglobal
  mass_p_nu=real((nf*nn)**3*f_nu)/npglobal_nu

  print*,'  from image',this_image(),'read',nplocal,nplocal_nu,' CDM/nu particles'
  sync all

  if (head) then
    print*,'  npglobal    =', npglobal
    print*,'  npglobal_nu =', npglobal_nu
    print*,'  omega_cdm   =', omega_m-omega_nu
    print*,'  omega_nu    =', omega_nu
    print*,'  f_cdm   =', f_cdm
    print*,'  f_nu    =', f_nu
    print*,'  mass_p     =', mass_p
    print*,'  mass_p_cdm =', mass_p_cdm
    print*,'  mass_p_nu  =', mass_p_nu

    print*,'  vsim2phy =',sim%vsim2phys, ' (km/s)/(1.0)'
    !print*,'  std_vf(a=',a_i,', r=',box/nf_global,'Mpc/h)',sqrt(3.)*sigma_vfi*sim%vsim2phys,'km/s'
    !print*,'  std_vc(a=',a_i,', r=',box/nc_global,'Mpc/h)',sqrt(3.)*sigma_vci*sim%vsim2phys,'km/s'
    !print*,'  std_vres',sqrt(3.)*sigma_vres*sim%vsim2phys,'km/s'
    print*,'  sigma_vi    =',sigma_vi,'(simulation unit)'
    print*,'  sigma_vi_nu =',sigma_vi_nu,'(simulation unit)'
  endif
  if (head) then
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, ''
  endif
  sync all
endsubroutine
