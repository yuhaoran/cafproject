!#define debug
subroutine kernel_f
  use variables
  implicit none
  save
  include 'fftw3.f'

  character(*),parameter :: dir_kern='../kernels/'
  integer itemp(3)
  real rtemp(2)
#ifdef debug
    real rho_f0(nfe+2,nfe,nfe,ncore)
#endif

  if (head) print*, 'fine kernel initialization'

  open(20,file=dir_kern//'wfxyzf.3.ascii',status='old',action='read')

  !x
  rho_f=0
  do k=1,nf_cutoff
  do j=1,nf_cutoff
  do i=1,nf_cutoff
    read(20,'(3i4,3e16.8)') itemp(1),itemp(2),itemp(3),rho_f(i,j,k,1),rtemp(1),rtemp(2)
    if (itemp(1).ne.i.or.itemp(2).ne.j.or.itemp(3).ne.k) stop 'error'
  enddo
  enddo
  enddo
  do j=2,nf_cutoff
    rho_f(1:nf_cutoff,nfe-j+2,1:nf_cutoff,1)=rho_f(1:nf_cutoff,j,1:nf_cutoff,1)
  enddo
  do i=2,nf_cutoff
    rho_f(nfe-i+2,1:nfe,1:nf_cutoff,1)=-rho_f(i,1:nfe,1:nf_cutoff,1)
  enddo
  do k=1,nf_cutoff
    rho_f(1:nfe,1:nfe,nfe-k+2,1)=rho_f(1:nfe,1:nfe,k,1)
  enddo
#ifdef debug
    rho_f0(:nfe,:,:,1)=rho_f(:nfe,:,:,1)
    sync all
#endif
  !print*, 'sum of rho_f', sum(rho_f)
  call sfftw_execute(plan_fft_fine)
  !print*, 'sum of rho_f', sum(rho_f)

#ifdef debug
    call sfftw_execute(plan_ifft_fine)
    rho_f(:nfe,:,:,1)=rho_f(:nfe,:,:,1)/real(nfe)**3
    print*, 'error on fft3d =',maxval(abs(rho_f(:nfe,:,:,1)-rho_f0(:nfe,:,:,1)))
    rho_f(nfe+1:,:,:,:)=0
    call sfftw_execute(plan_fft_fine)
    sync all
#endif
  do k=1,nfe
  do j=1,nfe
  do i=1,nfe/2+1
    kern_f(1,i,j,k)=rho_f(2*i,j,k,1)
  enddo
  enddo
  enddo

  !y
  rho_f=0
  rewind(20)
  do k=1,nf_cutoff
  do j=1,nf_cutoff
  do i=1,nf_cutoff
    read(20,*) itemp(1),itemp(2),itemp(3),rtemp(1),rho_f(i,j,k,1),rtemp(2)
    if (itemp(1).ne.i.or.itemp(2).ne.j.or.itemp(3).ne.k) stop 'error'
  enddo
  enddo
  enddo
  do j=2,nf_cutoff
    rho_f(1:nf_cutoff,nfe-j+2,1:nf_cutoff,1)=-rho_f(1:nf_cutoff,j,1:nf_cutoff,1)
  enddo
  do i=2,nf_cutoff
    rho_f(nfe-i+2,1:nfe,1:nf_cutoff,1)=rho_f(i,1:nfe,1:nf_cutoff,1)
  enddo
  do k=1,nf_cutoff
    rho_f(1:nfe,1:nfe,nfe-k+2,1)=rho_f(1:nfe,1:nfe,k,1)
  enddo
  call sfftw_execute(plan_fft_fine)
  do k=1,nfe
  do j=1,nfe
  do i=1,nfe/2+1
    kern_f(2,i,j,k)=rho_f(2*i,j,k,1)
  enddo
  enddo
  enddo

  !z
  rho_f=0
  rewind(20)
  do k=1,nf_cutoff
  do j=1,nf_cutoff
  do i=1,nf_cutoff
    read(20,*) itemp(1),itemp(2),itemp(3),rtemp(1),rtemp(2),rho_f(i,j,k,1)
    if (itemp(1).ne.i.or.itemp(2).ne.j.or.itemp(3).ne.k) stop 'error'
  enddo
  enddo
  enddo
  do j=2,nf_cutoff
    rho_f(1:nf_cutoff,nfe-j+2,1:nf_cutoff,1)=rho_f(1:nf_cutoff,j,1:nf_cutoff,1)
  enddo
  do i=2,nf_cutoff
    rho_f(nfe-i+2,1:nfe,1:nf_cutoff,1)=rho_f(i,1:nfe,1:nf_cutoff,1)
  enddo
  do k=1,nf_cutoff
    rho_f(1:nfe,1:nfe,nfe-k+2,1)=-rho_f(1:nfe,1:nfe,k,1)
  enddo
  call sfftw_execute(plan_fft_fine)
  do k=1,nfe
  do j=1,nfe
  do i=1,nfe/2+1
    kern_f(3,i,j,k)=rho_f(2*i,j,k,1)
  enddo
  enddo
  enddo

  close(20)

  !open(21,file=output_dir()//'kern_f'//output_suffix(),access='stream',status='replace')
  !write(21) kern_f
  !close(21)
  !print*, 'sum of kern_f', sum(kern_f)

  sync all

endsubroutine kernel_f
