! convert zip format into xv format (for DTFE & Voronoi)
! compute CIC density fields: delta_nbody.dat
! compute Voronoi density fields: delta_voronoi.dat
program convert_zip_to_xv
use parameters
implicit none

integer i,j,k,l
integer nplocal
integer,parameter :: ng=nf ! number of grids = number of fine cells
integer,parameter :: npnode=nf**3
real,parameter :: density_buffer=1.1
integer,parameter :: npmax=npnode*density_buffer
integer ind,dx,dxy,kg,mg,jg,ig,i0,j0,k0,itx,ity,itz,idx,imove,nf_shake,ibin
integer nshift,ifrom,ileft,iright,nlen,nlast,g(3)
real kr,kx(3), sincx,sincy,sincz,sinc, rbin

integer rhoc(nt,nt,nt,nnt,nnt,nnt)
real rho_f(0:ng+1,0:ng+1,0:ng+1)
integer idx1(3),idx2(3),ip,np
real mass_p,dx1(3),dx2(3),xv(6,npmax),xpos(3)
integer(izipx) x(3,npmax)
integer(izipv) v(3,npmax)

! checkpoint variables
integer,parameter :: nmax_redshift=100
integer cur_checkpoint, n_checkpoint
real z_checkpoint(nmax_redshift)

character (10) :: img_s, z_s
character (200) :: fn0,fn1,fn2,fn3,fn4
integer,parameter :: nexp=4

integer pp,i_s,j_s,k_s,ii_s,jj_s,kk_s,r
integer wp(npmax),pp_s
integer hoc_g(ng,ng,ng),ll_p(npmax) ! linked list
integer hoc_p(npmax),ll_g(ng**3) ! inverse linked list
real r2,r2min,gpos(3),hpos(3),dpos(3)
real(8) rho_tot
logical search

call geometry

if (head) print*, 'rho_f'


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


do cur_checkpoint=n_checkpoint,n_checkpoint
  fn0='.'//opath//'/node0/'//z2str(z_checkpoint(cur_checkpoint))//'zip0_0.dat'
  fn1='.'//opath//'/node0/'//z2str(z_checkpoint(cur_checkpoint))//'zip1_0.dat'
  fn2='.'//opath//'/node0/'//z2str(z_checkpoint(cur_checkpoint))//'zip2_0.dat'
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
  print*, 'nplocal =', sim%nplocal
  nplocal=sim%nplocal
  print*, 'v_r2i =', sim%v_r2i

  open(10,file=fn0,status='old',action='read',access='stream')
  read(10) x(:,:nplocal)
  close(10)
  open(11,file=fn1,status='old',action='read',access='stream')
  read(11) v(:,:nplocal)
  close(11)

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
      xv(1:3,ip) = nt*((/itx,ity,itz/)-1) + (/i,j,k/)-1 + (x(:,ip)+ishift+rshift)*x_resolution ! in unit of nc
      xv(4:6,ip) = v(:,ip)/sim%v_r2i
      xpos = xv(1:3,ip) * real(ng)/real(nc) - 0.5
      idx1=floor(xpos)+1
      idx2=idx1+1
      dx1=idx1-xpos
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

  ! buffer fine density
  rho_f(1,:,:)=rho_f(1,:,:)+rho_f(ng+1,:,:)
  rho_f(ng,:,:)=rho_f(ng,:,:)+rho_f(0,:,:)

  rho_f(:,1,:)=rho_f(:,1,:)+rho_f(:,ng+1,:)
  rho_f(:,ng,:)=rho_f(:,ng,:)+rho_f(:,0,:)

  rho_f(:,:,1)=rho_f(:,:,1)+rho_f(:,:,ng+1)
  rho_f(:,:,ng)=rho_f(:,:,ng)+rho_f(:,:,0)

  print*, sum(rho_f(1:ng,1:ng,1:ng)*1d0)/ng**3
  !rho_f(1:ng,1:ng,1:ng)=rho_f(1:ng,1:ng,1:ng)/(sum(rho_f(1:ng,1:ng,1:ng)*1d0)/ng**3)-1

  open(15,file='.'//opath//'node0/'//z2str(z_checkpoint(cur_checkpoint))//'delta_nbody.dat',status='replace',access='stream')
  write(15) rho_f(1:ng,1:ng,1:ng)
  close(15)
  print*, 'wrote file delta_nbody.dat'

  !write(16) xv(1:3,:nplocal) * sim%box*1000*(sim%h0/100)/real(sim%nt*sim%nnt*sim%nn) ! in unit of kpc

  print*, 'v_r2i =',sim%v_r2i

  open(16,file='.'//opath//'node0/'//z2str(z_checkpoint(cur_checkpoint))//'xv.dat',status='replace',access='stream')
  write(16) xv(:,:nplocal) / real(nc) ! in unit of box size
  close(16)
  print*, 'complete writing xv.dat'

  ! Voronoi
  print*, 'Voronoi'
  xv(1:3,:nplocal)=xv(1:3,:nplocal)*ng/nc ! in unit of required grid
  wp=1 ! weight -- number of grids assigned to the particle

  print*,  '  assign particles to grids'
  do ip=1,nplocal
    i=floor(xv(1,ip))+1
    j=floor(xv(2,ip))+1
    k=floor(xv(3,ip))+1
    ll_p(ip)=hoc_g(i,j,k)
    hoc_g(i,j,k)=ip
  enddo

  print*, '  assign non-empty grids to particles'
  do k=1,ng
  do j=1,ng
  do i=1,ng
    ig=ng**2*(k-1)+ng*(j-1)+i
    pp=hoc_g(i,j,k)
    do
      if (pp==0) exit
      hoc_p(pp)=ig
      pp=ll_p(pp)
    enddo
  enddo
  enddo
  enddo

  print*, '  assign empty grids to particles'
  do k=1,ng
  do j=1,ng
  do i=1,ng
    gpos=(/i,j,k/)-0.5
    ig=ng**2*(k-1)+ng*(j-1)+i
    if (hoc_g(i,j,k)==0) then
      r=0
      search=.true.
      r2min=1000*ng**2
      do ! search larger layer
        if (.not.search) exit
        r=r+1
        do k_s=k-r,k+r
        do j_s=j-r,j+r
        do i_s=i-r,i+r
          kk_s=modulo(k_s-1,ng)+1
          jj_s=modulo(j_s-1,ng)+1
          ii_s=modulo(i_s-1,ng)+1
          pp=hoc_g(ii_s,jj_s,kk_s)
          do
            if (pp==0) exit ! find next cell
            search=.false. ! found particle
            hpos=xv(1:3,pp)
            dpos=hpos-gpos
            dpos=modulo(dpos+ng/2,real(ng))-ng/2
            r2=sum(dpos**2)
            if (r2<r2min) then
              r2min=r2
              pp_s=pp
            endif
            pp=ll_p(pp)
          enddo
        enddo
        enddo
        enddo
      enddo ! finished searching current layer
      ! found nearest particle, exit searching loop
      ll_g(ig)=hoc_p(pp_s) ! let current grid point to pp_s's chain
      hoc_p(pp_s)=ig ! replace pp_s's hoc to be current grid
      wp(pp_s)=wp(pp_s)+1
    endif
  enddo
  enddo
  enddo

  print*, '  assign density'
  rho_f=0
  do ip=1,nplocal
    ig=hoc_p(ip)
    do
      if (ig==0) exit
      k=(ig-1)/ng**2+1
      j=(ig-1-ng**2*(k-1))/ng+1
      i=modulo(ig-1,ng)+1
      rho_f(i,j,k)=rho_f(i,j,k)+1./real(wp(ip));
      ig=ll_g(ig);
    enddo
  enddo
  rho_tot=sum(rho_f(1:ng,1:ng,1:ng)*1.d0)

  print*, 'Voronoi done'
  print*, 'nplocal =',nplocal
  print*, 'sum =',rho_tot
  print*, 'minval =',minval(rho_f(1:ng,1:ng,1:ng))
  rho_f(1:ng,1:ng,1:ng)=rho_f(1:ng,1:ng,1:ng)/(rho_tot/ng**3)-1
  open(11,file='.'//opath//'node0/'//z2str(z_checkpoint(cur_checkpoint))//'delta_voronoi.dat',status='replace',access='stream')
  write(11) rho_f(1:ng,1:ng,1:ng)
  close(11)
  print*, 'wrote into file delta_voronoi.dat'
enddo !cur_checkpoint

end
