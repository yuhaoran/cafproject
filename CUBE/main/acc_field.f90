! compute the acceleration field for all the checkpoints in redshifts.txt
program acc_field
  use variables
  implicit none
  save

  call initialize

  do cur_checkpoint=1,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    open(12,file=output_name('zip2'),status='old',action='read',access='stream')
    read(12) sim
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, 'zip format incompatable'
      close(12)
      stop
    endif
    read(12) rhoc(1:nt,1:nt,1:nt,:,:,:) ! coarse grid density
    close(12)
    mass_p=sim%mass_p
    if (head) print*, 'mass_p =',mass_p
    nplocal=sim%nplocal
    if (head) print*, 'nplocal =',nplocal
    open(10,file=output_name('zip0'),status='old',action='read',access='stream')
    read(10) xp(:,:nplocal) ! particle Eulerian positions
    close(10)

    call buffer_density
    call buffer_x
    a_mid=1./(1+z_checkpoint(cur_checkpoint))
    dt=1
    call pm_acceleration
  enddo

end
