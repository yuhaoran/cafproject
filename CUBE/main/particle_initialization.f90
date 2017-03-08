subroutine particle_initialization
  use variables
  implicit none
  save

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
  nplocal=sim%nplocal
  v_i2r=sim%v_i2r

  open(10,file=ic_name('zip0'),status='old',access='stream')
  read(10) x(:,:nplocal)
  close(10)
  open(11,file=ic_name('zip1'),status='old',access='stream')
  read(11) v(:,:nplocal)
  close(11)

#ifdef PID
    open(14,file=ic_name('zipid'),status='old',access='stream')
    read(14) pid(:,:nplocal)
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
  if (head) then
    print*, '  npglobal =', npglobal
    print*, '  mass_p=', mass_p
    !print*, '  max np per cell in this image =', maxval(rhoc)
    !print*, '  npmax per image =', npmax
    !print*, 'min v_i =', minval(v)*v_i2r
    !print*, 'max v_i =', maxval(v)*v_i2r
    print*, '  v_ir2 =',v_i2r,' (1.0)/1'
    print*, '  vsim2phy =',sim%vsim2phys, ' (km/s)/(1.0)'
    print*, '  velocity resolution =',v_i2r*sim%vsim2phys,' (km/s)'
    print*, '  velocity max =',2**(8*izipv)*v_i2r*sim%vsim2phys,' (km/s)'
    !print*, v_i2r
  endif
  sync all

endsubroutine
