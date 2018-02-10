subroutine initialize
  use variables
  use cubefft
  use pencil_fft
  implicit none
  save
  include 'fftw3.f'

  !call system('hostname')
  if (this_image()==1) then
    print*, ''
    print*, 'Coarray CUBE on',int(nn**3,2),'images  x',int(ncore,1),'cores'
    print*, ''
    print*, 'initialize'
    print*, '  call geometry'
  endif
  sync all

  call geometry

  if (head) print*, '  call create_cubefft_plan ng = ',ng
  call system_clock(t1,t_rate)
  call create_cubefft_plan
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
  sync all

  if (head) print*, '  call create_penfft_plan nfe = ',nfe
  call system_clock(t1,t_rate)
  call create_penfft_plan
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
  sync all

#ifdef NEUTRINOS
  neutrino_flag=.true.
#else
  neutrino_flag=.false.
#endif

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
    print*, 'checkpoint information'
    print*, '  output: ', opath
    print*, '  checkpoint at:'
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

  print*,'OpenMP information'
  nth=omp_get_num_procs()
  print*,'  omp_get_num_procs =',nth
  nth=omp_get_thread_limit()
  print*,'  omp_get_thread_limit =',nth
  call omp_set_num_threads(ncore)
  nth=omp_get_max_threads()
  print*,'  omp_get_max_threads =',nth
  print*, ''
endsubroutine
