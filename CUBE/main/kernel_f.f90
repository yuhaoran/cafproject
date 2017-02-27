!#define FPPKCORR
subroutine kernel_f
  use variables
  implicit none
  save
  include 'fftw3.f'

  character(*),parameter :: dir_kern='../kernels/'
  integer itemp(3),mfactor(3)
  real fk_table(nf_cutoff,nf_cutoff,nf_cutoff,3)

  if (head) print*, 'fine kernel initialization'

  open(20,file=dir_kern//'wfxyzf.3.ascii',status='old',action='read')
  do k=1,nf_cutoff
  do j=1,nf_cutoff
  do i=1,nf_cutoff
    read(20,*) itemp(1),itemp(2),itemp(3),fk_table(i,j,k,:)
    if (itemp(1).ne.i.or.itemp(2).ne.j.or.itemp(3).ne.k) stop 'error'
  enddo
  enddo
  enddo
  close(20)

  do i_dim=1,3
    rho_f=0
    mfactor=merge(-1,1,(/1,2,3/)==i_dim)
    print*, mfactor
    rho_f(:nf_cutoff,:nf_cutoff,:nf_cutoff,1)=fk_table(:,:,:,i_dim)
    rho_f(nfe-nf_cutoff+2:nfe,:,:,1)=mfactor(1)*rho_f(nf_cutoff:2:-1,:,:,1)
    rho_f(:,nfe-nf_cutoff+2:nfe,:,1)=mfactor(2)*rho_f(:,nf_cutoff:2:-1,:,1)
    rho_f(:,:,nfe-nf_cutoff+2:nfe,1)=mfactor(3)*rho_f(:,:,nf_cutoff:2:-1,1)
    call sfftw_execute(plan_fft_fine)
    kern_f(i_dim,:,:,:)=rho_f(2::2,:,:,1)
  enddo

  !open(21,file=output_dir()//'kern_f'//output_suffix(),access='stream',status='replace')
  !write(21) kern_f
  !close(21)
  !print*, 'sum of kern_f',&
  ! sum(kern_f(1,:,:,:)*1.d0),&
  ! sum(kern_f(2,:,:,:)*1.d0),&
  ! sum(kern_f(3,:,:,:)*1.d0)
  !print*, kern_f(1,:,1,1)
  sync all
endsubroutine kernel_f
