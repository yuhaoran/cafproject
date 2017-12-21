subroutine particle_initialization
  use variables
  use neutrinos

  implicit none
  save

  if (head) print*, 'particle_initialization'

  open(10,file=ic_name('info'),status='old',access='stream')
  read(10) sim
  close(10)
  nplocal=sim%nplocal
  nplocal_nu=sim%nplocal_nu
  sigma_vi=sim%sigma_vi
  sigma_vi_nu=sim%sigma_vi_nu

  if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
    print*, '  zip format incompatable'
    close(12)
    stop
  endif

  open(11,file=ic_name('xp'),status='old',access='stream')
  read(11) xp(:,:nplocal)
  close(11)

  open(11,file=ic_name('vp'),status='old',access='stream')
  read(11) vp(:,:nplocal)
  close(11)

  open(11,file=ic_name('np'),status='old',access='stream')
  read(11) rhoc(1:nt,1:nt,1:nt,:,:,:)
  close(11)

  open(11,file=ic_name('vc'),status='old',access='stream')
  read(11) vfield(:,1:nt,1:nt,1:nt,:,:,:)
  close(11)

#ifdef PID
    open(11,file=ic_name('id'),status='old',access='stream')
    read(11) pid(:nplocal)
    close(11)
#endif

  open(11,file=ic_name('xp_nu'),status='old',access='stream')
  read(11) xp_nu(:,:nplocal_nu)
  close(11)

  open(11,file=ic_name('vp_nu'),status='old',access='stream')
  read(11) vp_nu(:,:nplocal_nu)
  close(11)

  open(11,file=ic_name('np_nu'),status='old',access='stream')
  read(11) rhoc_nu(1:nt,1:nt,1:nt,:,:,:)
  close(11)

  open(11,file=ic_name('vc_nu'),status='old',access='stream')
  read(11) vfield_nu(:,1:nt,1:nt,1:nt,:,:,:)
  close(11)

#ifdef PID
      open(11,file=ic_name('id_nu'),status='old',access='stream')
      read(11) pid_nu(:nplocal_nu)
      close(11)
#endif

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
  sync all
endsubroutine
