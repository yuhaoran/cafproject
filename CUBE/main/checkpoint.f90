subroutine checkpoint
  use variables
  implicit none
  save

  integer(8),parameter :: blocksize=1024**2
  integer(8) nplow,nphigh,num_io

  if (head) print*, 'checkpoint'

  sim%nplocal=nplocal
  sim%a=a
  sim%t=t
  sim%tau=tau

  sim%istep=istep

  sim%dt_f_acc=dt_fine
  sim%dt_pp_acc=dt_pp
  sim%dt_c_acc=dt_coarse

  sim%cur_checkpoint=cur_checkpoint
  sim%cur_proj=cur_checkpoint
  sim%cur_halo=cur_checkpoint

  sim%mass_p=mass_p
  sim%vsim2phys=(1.5/a)*box*h0*sqrt(omega_m)/nf_global

  num_io=(nplocal-1)/blocksize+1

  if (head) print*, '  write in file:',output_name('zip2')
  open(12,file=output_name('zip2'),status='replace',access='stream')
  write(12) sim
  write(12) rhoc(1:nt,1:nt,1:nt,:,:,:)
  close(12)

  if (head) print*, '  write in file:',output_name('vfield')
  open(12,file=output_name('vfield'),status='replace',access='stream')
  write(12) vfield(:,1:nt,1:nt,1:nt,:,:,:)
  close(12)

  if (head) print*, '  write in file:',output_name('zip0')
  open(10,file=output_name('zip0'),status='replace',access='stream')
  !write(10) x(:,:nplocal)
  do i=1,num_io
    nplow=(i-1)*blocksize+1
    nphigh=min(i*blocksize,nplocal)
    write(10) xp(:,nplow:nphigh)
  enddo
  close(10)

  if (head) print*, '  write in file:',output_name('zip1')
  open(11,file=output_name('zip1'),status='replace',access='stream')
  !write(11) v(:,:nplocal)
  do i=1,num_io
    nplow=(i-1)*blocksize+1
    nphigh=min(i*blocksize,nplocal)
    write(11) vp(:,nplow:nphigh)
  enddo
  close(11)

#ifdef PID
  if (head) print*, '  write in file:',output_name('zipid')
  open(14,file=output_name('zipid'),status='replace',access='stream')
  !write(14) pid(:nplocal)
  do i=1,num_io
    nplow=(i-1)*blocksize+1
    nphigh=min(i*blocksize,nplocal)
    write(14) pid(nplow:nphigh)
  enddo
  close(14)
#endif

sync all
print*,'  image',this_image(),'write',nplocal,'particles'
npglobal=0
do i=1,nn**3
  npglobal=npglobal+nplocal[i]
enddo
sync all
if (head) print*, '  npglobal =',npglobal
sync all

endsubroutine
