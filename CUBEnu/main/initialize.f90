subroutine initialize
  use variables
  use cubefft
  use pencil_fft
  implicit none
  save
  include 'fftw3.f'

  if (this_image()==1) print*, 'initialize'
  if (this_image()==1) print*, 'call geometry'
  call geometry

  if (head) print*, 'call create_cubefft_plan ng = ',ng
  call create_cubefft_plan

  if (head) print*, 'call create_penfft_plan nfe = ',nfe
  call create_penfft_plan
  sync all
  ! omp_init


  call kernel_f
  call kernel_c

  istep=0
  t=0
  dt=0
  dt_old=0
  dt_vmax=1000
  dt_pp=1000
  dt_fine=1000
  dt_coarse=1000
  da=0
  tau=-3/sqrt(a_i)
  cur_checkpoint=1
  checkpoint_step=.false.
  final_step=.false.

  if (head) then
    print*, 'output: ', opath
    print*, 'checkpoint at:'
    open(16,file='redshifts.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      print*, z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
  endif
  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  call system('mkdir -p '//opath//'/image'//image2str(image))
  sync all

endsubroutine
