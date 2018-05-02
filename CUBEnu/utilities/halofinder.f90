program halofinder
  use omp_lib
  use parameters
  use buffer_grid_subroutines
  use buffer_particle_subroutines
  use variables, only : xp,vp,rhoc,vfield,i
  implicit none

  call geometry
  if (head) then
    print*, 'halofinder at:'
    open(16,file='../main/z_checkpoint.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_halofind(i)
      print*, z_halofind(i)
    enddo
    71 n_halofind=i-1
    close(16)
    print*,''
  endif

  sync all
  n_halofind=n_halofind[1]
  z_halofind(:)=z_halofind(:)[1]
  n_checkpoint=n_halofind
  z_checkpoint=z_halofind
  sync all

  do cur_halofind= 1,n_halofind
    cur_checkpoint=cur_halofind
    xp=0; vp=0; rhoc=0; vfield=0
    call particle_initialization
    call buffer_grid
    call buffer_x
    call buffer_v
    call halofind
  enddo
  sync all
  if (head) print*, 'halofinder done'

end
