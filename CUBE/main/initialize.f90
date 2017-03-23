!#define readkernel
!#define penffttest
#define mkdir

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

#ifdef penffttest
  if (head) print*, 'pencil_fft test...'
  call random_number(r3)
  call pencil_fft_forward
  call pencil_fft_backward
  if (head) print*, 'penfft test done.'
  sync all
#endif


#ifndef readkernel
  call kernel_f
  call kernel_c
#else
! now works only for single node
  open(14,file='cubep3m_kern_f.dat',status='old',access='stream')
  read(14) kern_f
  close(14)
  open(15,file='cubep3m_kern_c.dat',status='old',access='stream')
  read(15) kern_c
  close(15)
#endif

its=0
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
sync all

#ifdef mkdir
  call system('mkdir -p '//opath//'/node'//image2str(this_image()-1))
#endif
sync all

open(23,file='../velocity_conversion/vdisp.bin',access='stream')
read(23) vdisp
close(23)
vdisp(:,2)=vdisp(:,2)*vdisp(:,1)/1.5/box/h0/sqrt(omega_m)*real(nf_global)/sqrt(3.)
sigma_vi_old=interp_vdisp(a_i)
sigma_vi=interp_vdisp(a_i)
!if (head) print*,'sigma_vi',sigma_vi
sync all

endsubroutine
