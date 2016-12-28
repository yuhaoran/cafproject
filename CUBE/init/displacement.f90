#define linear_kbin

!#define remove_ny
program displacement
use penfft_fine
use powerspectrum
implicit none
save

! nc: coarse grid per node per dim
! nf: fine grid per node per dim

integer i,j,k,l,i_dim,dim_1,dim_2,dim_3
integer nplocal
integer,parameter :: npnode=nf**3 ! only true for this project
real,parameter :: density_buffer=1.2
integer,parameter :: npmax=npnode*density_buffer
integer ind,dx,dxy,kg,mg,jg,ig,i0,j0,k0,itx,ity,itz,idx,imove,nf_shake,ibin
integer nshift,ifrom,ileft,iright,nlen,nlast,g(3)
real kr,kx(3), sincx,sincy,sincz,sinc, rbin

real power(4,nbin), amp
real xi(10,nbin)[*],xi0(10,nbin)[*],xi1(10,nbin)[*],xi2(10,nbin)[*],xi3(10,nbin)[*]
real bk(nbin),n2k(nbin),wk(nbin)

integer rhoc(nt,nt,nt,nnt,nnt,nnt)
real rho_0(0:ng+1,0:ng+1,0:ng+1)[*]
real rho_grid(0:ng+1,0:ng+1,0:ng+1)[*]
real dsp(3,0:ng+1,0:ng+1,0:ng+1)[*]

real cube1(ng,ng,ng), cube2(ng,ng,ng), cube0(ng,ng,ng)

integer idx1(3),idx2(3),ip,np,icnode(3),irank
real mass_p,dx1(3),dx2(3),pos0(3),pos1(3),dpos(3)

integer(izipx) x(3,npmax)
integer(2)   pid(4,npmax)

! checkpoint variables
integer,parameter :: nmax_redshift=100
integer cur_checkpoint, n_checkpoint[*]
real z_checkpoint(nmax_redshift)[*]

character (10) :: img_s, z_s
character (200) :: fn0,fn1,fn2,fn3,fn4
integer,parameter :: nexp=4

complex cphi(ng*nn/2+1,ng,ngpen)
complex cdiv(ng*nn/2+1,ng,ngpen)
complex pdim, ekx(3)

if (head) then
  print*, 'Displacement field analysis on resolution:'
  print*, 'ng=',ng
endif

sync all

call create_penfft_fine_plan

if (head) then
  print*, 'checkpoint at:'
  open(16,file='../redshifts.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    print*, z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  print*,''
endif
sync all
n_checkpoint=n_checkpoint[1]
z_checkpoint(:)=z_checkpoint(:)[1]
sync all

do cur_checkpoint= n_checkpoint,n_checkpoint
  fn0='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip0_'//image2str(this_image()-1)//'.dat'
  fn1='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip1_'//image2str(this_image()-1)//'.dat'
  fn2='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip2_'//image2str(this_image()-1)//'.dat'
  fn3='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip3_'//image2str(this_image()-1)//'.dat'
  fn4='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zipid_'//image2str(this_image()-1)//'.dat'

  if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

  open(12,file=fn2,status='old',action='read',access='stream')
  read(12) sim
  ! check zip format and read rhoc
  if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
    print*, 'zip format incompatable'
    close(12)
    stop
  endif
  read(12) rhoc ! coarse grid density
  close(12)

  mass_p=sim%mass_p
  if (head) print*, 'mass_p =',mass_p

  !print*, sum(rhoc), sim%nplocal
  nplocal=sim%nplocal
  if (head) print*, 'nplocal =',nplocal

  open(10,file=fn0,status='old',action='read',access='stream')
  read(10) x(:,:nplocal) ! particle Eulerian positions
  close(10)

  open(14,file=fn4,status='old',action='read',access='stream')
  read(14) pid(:,:nplocal) ! particle Lagrangian positions
  close(14)

  ! mesh dsp and rho
  rho_0=0 ! CIC number count
  rho_grid=0 ! CIC final density
  dsp=0 ! CIC disp
  !dsp_e=0
  nlast=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
  if (head) print*, 'CIC interpolation on tile'
  do k=1,nt
  do j=1,nt
  do i=1,nt
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np
      ip=nlast+l

      irank=pid(1,ip)
      icnode(3)=irank/(nn**2)
      icnode(2)=(irank-nn**2*icnode(3))/nn
      icnode(1)=mod(irank,nn)

      pos0=real(ng)*icnode + real(ng)*(real(pid(2:4,ip))+0.5+2**15)/2**16
      pos1=nt*((/itx,ity,itz/)-1) + (/i,j,k/)-1 + ((x(:,ip)+ishift)+rshift)*x_resolution
      pos1=real(ng)*((/icx,icy,icz/)-1) + pos1*real(ng)/real(nc)
      dpos=pos1-pos0
      dpos=modulo(dpos+ng*nn/2,real(ng*nn))-ng*nn/2
      !dpos=dpos*real(nf)/real(ng)

!#ifdef remove_ny
!      pos0=pos0-0.5
!#else
      pos0=pos0-0.5 +0.5 ! dsp- "+0.5" assigns particle to the half-right grid
                         ! then the divergence will be done between current- and right-grid.
!#endif
      idx1=floor(pos0)+1
      idx2=idx1+1
      dx1=idx1-pos0
      dx2=1-dx1

      ! CIC dsp_{-}
      dsp(:,idx1(1),idx1(2),idx1(3))=dsp(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*dpos
      dsp(:,idx2(1),idx1(2),idx1(3))=dsp(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*dpos
      dsp(:,idx1(1),idx2(2),idx1(3))=dsp(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*dpos
      dsp(:,idx1(1),idx1(2),idx2(3))=dsp(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*dpos
      dsp(:,idx1(1),idx2(2),idx2(3))=dsp(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*dpos
      dsp(:,idx2(1),idx1(2),idx2(3))=dsp(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*dpos
      dsp(:,idx2(1),idx2(2),idx1(3))=dsp(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*dpos
      dsp(:,idx2(1),idx2(2),idx2(3))=dsp(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*dpos
      ! CIC number count
      rho_0(idx1(1),idx1(2),idx1(3))=rho_0(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
      rho_0(idx2(1),idx1(2),idx1(3))=rho_0(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
      rho_0(idx1(1),idx2(2),idx1(3))=rho_0(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
      rho_0(idx1(1),idx1(2),idx2(3))=rho_0(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
      rho_0(idx1(1),idx2(2),idx2(3))=rho_0(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
      rho_0(idx2(1),idx1(2),idx2(3))=rho_0(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
      rho_0(idx2(1),idx2(2),idx1(3))=rho_0(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
      rho_0(idx2(1),idx2(2),idx2(3))=rho_0(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
      ! CIC final density
      pos1=pos1-0.5
      idx1=floor(pos1)+1
      idx2=idx1+1
      dx1=idx1-pos1
      dx2=1-dx1
      rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
      rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
      rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
      rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
      rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
      rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
      rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
      rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p

    enddo
    nlast=nlast+np
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  
  print*, 'sum of rho_grid',sum(rho_grid(1:ng,1:ng,1:ng)*1d0)
  print*, 'Start sync from buffer regions'  
  print*, '  dsp'
  sync all
  ! buffer dsp
  dsp(:,1,:,:)=dsp(:,1,:,:)+dsp(:,ng+1,:,:)[image1d(inx,icy,icz)]
  dsp(:,ng,:,:)=dsp(:,ng,:,:)+dsp(:,0,:,:)[image1d(ipx,icy,icz)]
  sync all
  dsp(:,:,1,:)=dsp(:,:,1,:)+dsp(:,:,ng+1,:)[image1d(icx,iny,icz)]
  dsp(:,:,ng,:)=dsp(:,:,ng,:)+dsp(:,:,0,:)[image1d(icx,ipy,icz)]
  sync all
  dsp(:,:,:,1)=dsp(:,:,:,1)+dsp(:,:,:,ng+1)[image1d(icx,icy,inz)]
  dsp(:,:,:,ng)=dsp(:,:,:,ng)+dsp(:,:,:,0)[image1d(icx,icy,ipz)]
  sync all

  ! buffer fine density
  print*, '  rho_0'
  rho_0(1,:,:)=rho_0(1,:,:)+rho_0(ng+1,:,:)[image1d(inx,icy,icz)]
  rho_0(ng,:,:)=rho_0(ng,:,:)+rho_0(0,:,:)[image1d(ipx,icy,icz)]
  sync all
  rho_0(:,1,:)=rho_0(:,1,:)+rho_0(:,ng+1,:)[image1d(icx,iny,icz)]
  rho_0(:,ng,:)=rho_0(:,ng,:)+rho_0(:,0,:)[image1d(icx,ipy,icz)]
  sync all
  rho_0(:,:,1)=rho_0(:,:,1)+rho_0(:,:,ng+1)[image1d(icx,icy,inz)]
  rho_0(:,:,ng)=rho_0(:,:,ng)+rho_0(:,:,0)[image1d(icx,icy,ipz)]
  sync all
  print*, '  rho_grid'
  rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)[image1d(inx,icy,icz)]
  rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)[image1d(ipx,icy,icz)]
  sync all
  rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)[image1d(icx,iny,icz)]
  rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)[image1d(icx,ipy,icz)]
  sync all
  rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)[image1d(icx,icy,inz)]
  rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)[image1d(icx,icy,ipz)]
  sync all
  
  print*, 'sum of rho_grid = '
  print*, sum(rho_grid(1:ng,1:ng,1:ng)*1d0)
  rho_grid=rho_grid/(sum(rho_grid(1:ng,1:ng,1:ng)*1d0)/ng**3)-1
  cube2=rho_grid(1:ng,1:ng,1:ng)

  if (head) print*,'Write delta_N into file'
  open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/delta_nbody.dat',access='stream')
  write(15) cube2
  close(15)

  do i_dim=1,3
    dsp(i_dim,1:ng,1:ng,1:ng)=dsp(i_dim,1:ng,1:ng,1:ng)/rho_0(1:ng,1:ng,1:ng)
    print*, 'dsp: dim',int(i_dim,1),'min,max values ='
    print*, minval(dsp(i_dim,1:ng,1:ng,1:ng)), maxval(dsp(i_dim,1:ng,1:ng,1:ng))
  enddo
  
  if (head) print*,'Write dsp into file'
  open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/dsp.dat',access='stream')
  write(15) dsp(1,1:ng,1:ng,1:ng)
  write(15) dsp(2,1:ng,1:ng,1:ng)
  write(15) dsp(3,1:ng,1:ng,1:ng)
  close(15)

  if (head) print*,'Start reconstructing delta_R'
  cphi=0
  cdiv=0
  do i_dim=1,3
    if (head) print*,'Start working on dim',int(i_dim,1)
    cube=dsp(i_dim,1:ng,1:ng,1:ng)
    if (head) print*,'start forward tran'
    call fft_cube2pencil_fine
    if (head) print*,'start transpose'
    call trans_zxy2xyz_fine
    if (head) print*,'loop over k'
    ! cx is the fourier of dsp(i_dim,1:ng,1:ng,1:ng)
    do k=1,ngpen
    do j=1,ng
    do i=1,ng*nn/2+1
      kg=(nn*(icz-1)+icy-1)*ngpen+k
      jg=(icx-1)*ng+j
      ig=i
      kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
      kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
      ekx=exp(2*pi*(0,1)*kx/ng)
      dim_1=i_dim
      dim_2=mod(dim_1,3)+1
      dim_3=mod(dim_2,3)+1
      pdim=(ekx(dim_1)-1)*(ekx(dim_2)+1)*(ekx(dim_3)+1)/4
      cphi(i,j,k)=cphi(i,j,k)+cx(i,j,k)*pdim/(-4*sum(sin(pi*kx/ng)**2)+0.000001)
      cdiv(i,j,k)=cdiv(i,j,k)+cx(i,j,k)*pdim
    enddo
    enddo
    enddo
  enddo ! i_dim
  if (head) then
    cphi(1,1,1)=0
    cdiv(1,1,1)=0
  endif
  sync all


  !! reconstructed delta
  cx=cdiv
  if (head) print*,'start backward tran'
  call trans_xyz2zxy_fine
  if (head) print*,'start transpose'
  call ifft_pencil2cube_fine
  cube1=-cube
  if (head) print*,'Write delta_R into file'
  open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/delta_reco.dat',status='replace',access='stream')
  write(15) cube1
  close(15)
  sync all

  cx=cphi
  if (head) print*,'start backward tran'
  call trans_xyz2zxy_fine
  if (head) print*,'start transpose'
  call ifft_pencil2cube_fine
  cube1=cube
  if (head) print*,'Write delta_R into file'
  open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/phi_E.dat',status='replace',access='stream')
  write(15) cube1
  close(15)
  sync all

  if (head) print*,'Read delta_L from file'
  !! linear delta
  open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/delta_L.dat',status='old',access='stream')
  read(15) cube0
  close(15)
  sync all

  if (head) print*,'Main: call cross_power LR____________________'
  xi=0 ! force cross_power use generated Wiener filter
  call cross_power(xi,cube0,cube1)
  open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/xi_LR.dat',status='replace',access='stream')
  write(15) xi
  close(15)
  call system('mv power_fields.dat delta_wiener_LR.dat')

  if (head) print*,'Main: call cross_power LN____________________'
  ! xi from last step is input to cross_power for filtering delta_N
  call cross_power(xi,cube0,cube2)
  open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/xi_LN.dat',status='replace',access='stream')
  write(15) xi
  close(15)
  call system('mv power_fields.dat delta_wiener_LN.dat')

enddo !cur_checkpoint
if (head) print*, 'destroying fft plans'
call destroy_penfft_fine_plan
print*,'displacement done'

endprogram
