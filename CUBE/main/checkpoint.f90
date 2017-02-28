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

sim%dt_f_acc=dt_fine(1)
sim%dt_pp_acc=dt_pp(1)
sim%dt_c_acc=dt_coarse(1)

sim%cur_checkpoint=cur_checkpoint
sim%cur_proj=cur_checkpoint
sim%cur_halo=cur_checkpoint

sim%mass_p=mass_p
sim%v_i2r=v_i2r
sim%shake_offset=0

sim%box=box
sim%rank=this_image()-1
sim%nn=int(nn,2)
sim%nnt=int(nnt,2)
sim%nt=int(nt,2)
sim%ncell=int(ncell,2)
sim%ncb=int(ncb,2)
sim%izipx=int(izipx,1)
sim%izipv=int(izipv,1)

sim%h0=h0
sim%omega_m=omega_m
sim%omega_l=omega_l
!sim%s8

sim%m_neu(1:3)=0
sim%vsim2phys=1.0/(300.*sqrt(omega_m)*box/a/2./ nf/nn)
sim%z_i=z_i

open(12,file=output_name('zip2'),status='replace',access='stream')
write(12) sim
write(12) rhoc(1:nt,1:nt,1:nt,:,:,:)
close(12)

open(10,file=output_name('zip0'),status='replace',access='stream')
open(11,file=output_name('zip1'),status='replace',access='stream')
!do itz=1,nnt ! loop over tile
!do ity=1,nnt
!do itx=1,nnt
!  do iz=1,nt
!  do iy=1,nt
!    nlast=cum(nt,iy,iz,itx,ity,itz)
!    nfirst=cum(0,iy,iz,itx,ity,itz)+1
!    write(10) x(:,nfirst:nlast)
!    write(11) v(:,nfirst:nlast)
!  enddo
!  enddo
!enddo
!enddo
!enddo
write(10) x(:,:nplocal)
write(11) v(:,:nplocal)

close(10)
close(11)

#ifdef PID
open(14,file=output_name('zipid'),status='replace',access='stream')
write(14) pid(:,:nplocal)
close(14)
#endif


endsubroutine
