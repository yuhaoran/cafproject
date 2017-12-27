subroutine checkpoint
  use variables
  use neutrinos
  implicit none
  save

  if (head) print*, 'checkpoint'

  sim%nplocal=nplocal
  sim%nplocal_nu=nplocal_nu
  sim%a=a
  sim%t=t
  sim%tau=tau

  sim%istep=istep

  sim%dt_f_acc=dt_fine
  sim%dt_pp_acc=dt_pp
  sim%dt_c_acc=dt_coarse
  sim%dt_vmax=dt_vmax
  sim%dt_vmax_nu=dt_vmax_nu

  sim%cur_checkpoint=cur_checkpoint

  sim%mass_p=mass_p
  sim%vsim2phys=(1.5/a)*box*h0*100.*sqrt(omega_m)/nf_global

  open(11,file=output_name('info'),status='replace',access='stream')
  write(11) sim
  close(11)

  open(11,file=output_name('xp'),status='replace',access='stream')
  write(11) xp(:,:nplocal)
  close(11)

  open(11,file=output_name('vp'),status='replace',access='stream')
  write(11) vp(:,:nplocal)
  close(11)

  open(11,file=output_name('np'),status='replace',access='stream')
  write(11) rhoc(1:nt,1:nt,1:nt,:,:,:)
  close(11)

  open(11,file=output_name('vc'),status='replace',access='stream')
  write(11) vfield(:,1:nt,1:nt,1:nt,:,:,:)
  close(11)

#ifdef PID
  open(11,file=output_name('id'),status='replace',access='stream')
  write(11) pid(:nplocal)
  close(11)
#endif

  open(11,file=output_name('xp_nu'),status='replace',access='stream')
  write(11) xp_nu(:,:nplocal)
  close(11)

  open(11,file=output_name('vp_nu'),status='replace',access='stream')
  write(11) vp_nu(:,:nplocal)
  close(11)

  open(11,file=output_name('np_nu'),status='replace',access='stream')
  write(11) rhoc_nu(1:nt,1:nt,1:nt,:,:,:)
  close(11)

  open(11,file=output_name('vc_nu'),status='replace',access='stream')
  write(11) vfield_nu(:,1:nt,1:nt,1:nt,:,:,:)
  close(11)

#ifdef EID
  open(11,file=output_name('id_nu'),status='replace',access='stream')
  write(11) pid_nu(:nplocal)
  close(11)
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
