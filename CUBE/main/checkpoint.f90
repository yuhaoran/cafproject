subroutine checkpoint
use variables
implicit none
save

integer nfirst,nlast

if (head) print*, 'checkpoint'

sim%nplocal=nplocal
sim%a=a
sim%t=t
sim%tau=tau

sim%nts=its

sim%dt_f_acc=dt_fine
sim%dt_pp_acc=dt_pp
sim%dt_c_acc=dt_coarse

sim%cur_checkpoint=cur_checkpoint
sim%cur_proj=cur_checkpoint
sim%cur_halo=cur_checkpoint

sim%mass_p=mass_p
sim%v_i2r=v_i2r

if (head) print*, '  write in file:',output_name('zip2')
open(12,file=output_name('zip2'),status='replace',access='stream')
write(12) sim
write(12) rhoc(1:nt,1:nt,1:nt,:,:,:)
close(12)

if (head) print*, '  write in file:',output_name('zip0')
open(10,file=output_name('zip0'),status='replace',access='stream')
write(10) x(:,:nplocal)
close(10)

if (head) print*, '  write in file:',output_name('zip1')
open(11,file=output_name('zip1'),status='replace',access='stream')
write(11) v(:,:nplocal)
close(11)

#ifdef PID
if (head) print*, '  write in file:',output_name('zipid')
open(14,file=output_name('zipid'),status='replace',access='stream')
write(14) pid(:,:nplocal)
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
