program rhof
use parameters
implicit none
save

integer i,j,k,l
integer nplocal
integer,parameter :: ng=nf/8
integer,parameter :: npnode=nf**3
real,parameter :: density_buffer=1.5
integer,parameter :: npmax=npnode*density_buffer
integer ind,dx,dxy,kg,mg,jg,ig,i0,j0,k0,itx,ity,itz,idx,imove,nf_shake,ibin
integer nshift,ifrom,ileft,iright,nlen,nlast,g(3)
real kr,kx(3), sincx,sincy,sincz,sinc, rbin

integer rhoc(nt,nt,nt,nnt,nnt,nnt)
real rho_f(0:ng+1,0:ng+1,0:ng+1)[*]
integer idx1(3),idx2(3),ip,np
real mass_p,dx1(3),dx2(3),tempx(3)

integer(izipx) x(3,npmax)[*]

! checkpoint variables
integer,parameter :: nmax_redshift=100
integer cur_checkpoint, n_checkpoint[*]
real z_checkpoint(nmax_redshift)[*]

character (10) :: img_s, z_s
character (200) :: fn0,fn1,fn2,fn3,fn4
integer,parameter :: nexp=4


if (head) print*, 'rho_f'
sync all


if (head) then
  print*, 'checkpoint at:'
  open(16,file='../redshifts.txt',status='old')
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

do cur_checkpoint=n_checkpoint,n_checkpoint
  fn0='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip0_'//image2str(this_image()-1)//'.dat'
  fn1='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip1_'//image2str(this_image()-1)//'.dat'
  fn2='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip2_'//image2str(this_image()-1)//'.dat'
  fn3='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip3_'//image2str(this_image()-1)//'.dat'
  
  open(12,file=fn2,status='old',action='read',access='stream')
  read(12) sim
  ! check zip format and read rhoc
  if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
    print*, 'zip format incompatable'
    close(12)
    stop
  endif
  read(12) rhoc
  close(12)

  mass_p=sim%mass_p

  !print*, sum(rhoc), sim%nplocal
  nplocal=sim%nplocal

  open(10,file=fn0,status='old',action='read',access='stream')
  read(10) x(:,:nplocal)
  close(10)

  ! particle mesh
  rho_f=0
  nlast=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
  do k=1,nt
  do j=1,nt
  do i=1,nt
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np
      ip=nlast+l
      tempx = nt*((/itx,ity,itz/)-1) + (/i,j,k/)-1 + (x(:,ip)+ishift+rshift)*x_resolution
      tempx = tempx * real(ng)/real(nc) - 0.5
      idx1=floor(tempx)+1
      idx2=idx1+1
      dx1=idx1-tempx
      dx2=1-dx1

      rho_f(idx1(1),idx1(2),idx1(3))=rho_f(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
      rho_f(idx2(1),idx1(2),idx1(3))=rho_f(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
      rho_f(idx1(1),idx2(2),idx1(3))=rho_f(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
      rho_f(idx1(1),idx1(2),idx2(3))=rho_f(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
      rho_f(idx1(1),idx2(2),idx2(3))=rho_f(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
      rho_f(idx2(1),idx1(2),idx2(3))=rho_f(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
      rho_f(idx2(1),idx2(2),idx1(3))=rho_f(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
      rho_f(idx2(1),idx2(2),idx2(3))=rho_f(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p

    enddo
    nlast=nlast+np
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  
  print*, sum(rho_f*1d0)

  sync all

  ! buffer fine density
  rho_f(1,:,:)=rho_f(1,:,:)+rho_f(ng+1,:,:)[image1d(inx,icy,icz)]
  rho_f(ng,:,:)=rho_f(ng,:,:)+rho_f(0,:,:)[image1d(ipx,icy,icz)]
  sync all
  rho_f(:,1,:)=rho_f(:,1,:)+rho_f(:,ng+1,:)[image1d(icx,iny,icz)]
  rho_f(:,ng,:)=rho_f(:,ng,:)+rho_f(:,0,:)[image1d(icx,ipy,icz)]
  sync all
  rho_f(:,:,1)=rho_f(:,:,1)+rho_f(:,:,ng+1)[image1d(icx,icy,inz)]
  rho_f(:,:,ng)=rho_f(:,:,ng)+rho_f(:,:,0)[image1d(icx,icy,ipz)]
  sync all
  print*, sum(rho_f(1:ng,1:ng,1:ng)*1d0)/ng**3
  rho_f(1:ng,1:ng,1:ng)=rho_f(1:ng,1:ng,1:ng)/(sum(rho_f(1:ng,1:ng,1:ng)*1d0)/ng**3)-1

  open(15,file='delta_cdm.dat',access='stream')
  write(15) rho_f(1:ng,1:ng,1:ng)
  close(15)

enddo !cur_checkpoint




contains


end
