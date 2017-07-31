subroutine particle_initialization
  use variables
  implicit none
  save

  integer(8),parameter :: blocksize=1024**2
  integer(8) nplow,nphigh,num_io

  if (head) print*, 'particle_initialization'

  open(12,file=ic_name('zip2'),status='old',access='stream')
    read(12) sim ! 128 bytes header
    ! check zip format
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, '  zip format incompatable'
      close(12)
      stop
    endif
    read(12) rhoc(1:nt,1:nt,1:nt,:,:,:)
  close(12)

  open(15,file=ic_name('vfield'),status='old',access='stream')
    read(15) vfield(:,1:nt,1:nt,1:nt,:,:,:)
  close(15)

  nplocal=sim%nplocal

  num_io=(nplocal-1)/blocksize+1
  if (head) print*, 'number of io',num_io

  open(10,file=ic_name('zip0'),status='old',access='stream')
  !read(10) x(:,:nplocal)
  do i=1,num_io
    nplow=(i-1)*blocksize+1
    nphigh=min(i*blocksize,nplocal)
    read(10) xp(:,nplow:nphigh)
  enddo
  close(10)

  open(11,file=ic_name('zip1'),status='old',access='stream')
  !read(11) v(:,:nplocal)
  do i=1,num_io
    nplow=(i-1)*blocksize+1
    nphigh=min(i*blocksize,nplocal)
    read(11) vp(:,nplow:nphigh)
  enddo
  close(11)

  !print*,nplocal,ip
  !print*,x(:,1),x(:,nplocal)
  !print*,v(:,1),v(:,nplocal)
  !print*,'done'
  !stop

#ifdef PID
  open(14,file=ic_name('zipid'),status='old',access='stream')
  !read(14) pid(:,:nplocal)
  do i=1,num_io
    nplow=(i-1)*blocksize+1
    nphigh=min(i*blocksize,nplocal)
    read(14) pid(nplow:nphigh)
  enddo
  close(14)
#endif
  a=a_i
  npglobal=0
  sync all

  do i=1,nn**3
    npglobal=npglobal+nplocal[i]
  enddo
  mass_p=real((nf*nn)**3)/npglobal

  print*,'  from image',this_image(),'read',nplocal,' particles'
  sync all

  !sigma_vfi=interp_sigmav(a_i,box/nf_global) ! sigma(v) on scale of fine grid
  !sigma_vci=interp_sigmav(a_i,box/nc_global) ! sigma(v) on scale of coarse grid

  !print*, sigma_vfi, sim%vsim2phys; stop

  sigma_vi=sim%sigma_vi
  !sigma_vres=sqrt(sigma_vfi**2-sigma_vci**2) ! sigma(v) residual
  !sigma_vci_old=sigma_vci
  !sigma_vfi_old=sigma_vfi
  !sigma_vres_old=sigma_vres

  if (head) then
    print*,'  npglobal =', npglobal
    print*,'  mass_p=', mass_p
    print*,'  vsim2phy =',sim%vsim2phys, ' (km/s)/(1.0)'
    !print*,'  std_vf(a=',a_i,', r=',box/nf_global,'Mpc/h)',sqrt(3.)*sigma_vfi*sim%vsim2phys,'km/s'
    !print*,'  std_vc(a=',a_i,', r=',box/nc_global,'Mpc/h)',sqrt(3.)*sigma_vci*sim%vsim2phys,'km/s'
    !print*,'  std_vres',sqrt(3.)*sigma_vres*sim%vsim2phys,'km/s'
    print*,'  sigma_vi =',sigma_vi,'(simulation unit)'
  endif
  sync all
endsubroutine
