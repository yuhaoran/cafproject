!#define twoptest
!#define zipconvert
!#define READCUBEP3M
subroutine particle_initialization
use variables
implicit none
save

integer np255,ileft,iright,ntemp
#ifdef READCUBEP3M
  real scalefactor0,time0,tau0,nts0,dt0(1:3),mass0,v_r2i,shake_offset(3)
  integer cur_check(1:3), cum_cubep3m(nc,nc,nc), nlen, nlast
  integer(1) cubep3m_zip2(nc,nc,nc)
#endif


if (head) print*, 'particle_initialization'

#ifdef READCUBEP3M
  print*, 'read in cubep3m checkpoint'
  open(10,file='./init/LOS1/node0/20.000zip0_0.dat',status='old',access='stream')
  open(11,file='./init/LOS1/node0/20.000zip1_0.dat',status='old',access='stream')
  open(12,file='./init/LOS1/node0/20.000zip2_0.dat',status='old',access='stream')
  open(13,file='./init/LOS1/node0/20.000zip3_0.dat',status='old',access='stream')
  read(11) nplocal,scalefactor0,time0,tau0,nts0,dt0(1:3),cur_check(1:3),mass0,v_r2i,shake_offset
  read(10) nplocal,scalefactor0,time0,tau0,nts0,dt0(1:3),cur_check(1:3),mass0,v_r2i,shake_offset
  print*, 'header:'
  print*, 'nplocal,scalefactor0,time0,tau0,nts0'
  print*, nplocal,scalefactor0,time0,tau0,nts0
  print*, 'dt'
  print*, dt0(1:3)
  print*, 'cur_checkpoint'
  print*, cur_check(1:3)
  print*,'mass0,v_r2i,shake_offset'
  print*, mass0,v_r2i,shake_offset

  ! assume zip3 is empty
  x_new=0; v_new=0; x=0; v=0
  read(12) cubep3m_zip2
# ifndef zipconvert
    v_i2r=1./v_r2i
    read(10) x_new(:,:nplocal)
    read(11) v_new(:,:nplocal)
# else
    v_i2r=1./v_r2i/2**((izipv-2)*8)
    read(10) xic_new(:,:nplocal)
    read(11) vic_new(:,:nplocal)
    x_new(:,:nplocal)=int(xic_new(:,:nplocal),kind=izipx)*2**((izipx-1)*8)+(2**((izipx-1)*7)-1)
    print*, ''
    print*, xic_new(:,1),'---->',x_new(:,1)
    v_new(:,:nplocal)=vic_new(:,:nplocal)*2**((izipv-2)*8)+(2**((izipv-2)*7)-1)
# endif
!vmax=2**(izipv*8-1)*v_i2r
vmax=maxval(v_new(:,:nplocal),2)*v_i2r
print*, 'vmax_ic',maxval(vic_new(:,:nplocal),2)/v_r2i
print*, 'vmax',vmax
  ! convert density and xv
  rhoc=0
  cum_cubep3m=cumsum_cubep3m(cubep3m_zip2)
  iright=0
  do itz=1,nnt!*0+1
  do ity=1,nnt!*0+1
  do itx=1,nnt!*0+1
    rhoc(1:nt,1:nt,1:nt,itx,ity,itz)=cubep3m_zip2(nt*(itx-1)+1:nt*itx,nt*(ity-1)+1:nt*ity,nt*(itz-1)+1:nt*itz)
    nptile(itx,ity,itz)=sum(rhoc(:,:,:,itx,ity,itz))
    do iz=1,nt!*0+1
    do iy=1,nt!*0+1
      ileft=iright+1
      nlen=sum(rhoc(1:nt,iy,iz,itx,ity,itz))
      iright=ileft+nlen-1
      nlast=cum_cubep3m(nt*itx,nt*(ity-1)+iy,nt*(itz-1)+iz)
      x(:,ileft:iright)=x_new(:,nlast-nlen+1:nlast)
      v(:,ileft:iright)=v_new(:,nlast-nlen+1:nlast)
    enddo
    enddo
  enddo
  enddo
  enddo
!v=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !print*, 'check', iright
  !print*, 'check', cum_cubep3m(nc,nc,nc)
  !print*, sum(int(x(:,:nplocal),8)), sum(int(x_new(:,:nplocal),8))
  !print*, sum(int(v(:,:nplocal),8)), sum(int(v_new(:,:nplocal),8))
  !print*, v(:,1), v_new(:,1)
  print*, 'nptile, sum(nptile) ='
  print*,  nptile
  print*,  sum(nptile)
#else
  fn0=ipath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip0_'//image2str(this_image()-1)//'.dat'
  fn1=ipath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip1_'//image2str(this_image()-1)//'.dat'
  fn2=ipath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip2_'//image2str(this_image()-1)//'.dat'
  fn3=ipath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zip3_'//image2str(this_image()-1)//'.dat'
  fn4=ipath//'/node'//image2str(this_image()-1)//'/'//z2str(z_i)//'zipid_'//image2str(this_image()-1)//'.dat'
  open(12,file=fn2,access='stream')
  !open(13,file='./init/2lpt/z_50_ZA/zip3.dat',status='old',access='stream')

  read(12) sim ! 128 bytes header
  ! check zip format
  if (sim%izipx/=izipx .or. sim%izipv/=izipv .or. sim%izip2/=4) then
    print*, 'zip format incompatable'
    close(12)
    stop
  endif

  read(12) rhoc(1:nt,1:nt,1:nt,:,:,:)
  close(12)
  nplocal=sim%nplocal
  v_i2r=1/sim%v_r2i
  !read(12) rhoc_i1(:,:,:,:,:,:)
  !if (count(rhoc_i1==-1)>0) then
  !  read(13) rhoc_i4(:count(rhoc_i1==-1))
  !endif
  !close(13)
  open(10,file=fn0,status='old',access='stream')
  open(11,file=fn1,status='old',access='stream')
  read(10) x(:,:nplocal)
  read(11) v(:,:nplocal)
  close(10)
  close(11)
#ifdef PID
  open(14,file=fn4,status='old',access='stream')
  read(14) pid(:,:nplocal)
  close(14)
#endif

  !iright=0
  !do itz=1,nnt
  !do ity=1,nnt
  !do itx=1,nnt
  !  ileft=iright+1
  !  np255=count(rhoc_i1(:,:,:,itx,ity,itz)==-1)
  !  iright=ileft+np255-1
  !  nptile(itx,ity,itz)=sum(rhoc_i1(:,:,:,itx,ity,itz)+int(-128,1)+128)-255*np255+sum(rhoc_i4(ileft:iright))
  !enddo
  !enddo
  !enddo

  !print*, 'iright =', iright
  !print*, 'nptile, sum(nptile) ='
  !print*,  nptile
  !print*,  sum(nptile)

  !rhoc(1:nt,1:nt,1:nt,:,:,:)=(rhoc_i1+int(-128,1))+128
  !ntemp=count(rhoc(1:nt,1:nt,1:nt,:,:,:)==255)
  !rhoc(1:nt,1:nt,1:nt,:,:,:)=unpack(rhoc_i4(:ntemp),rhoc(1:nt,1:nt,1:nt,:,:,:)==255,rhoc(1:nt,1:nt,1:nt,:,:,:))

#endif

#ifdef twoptest
  nplocal=2
  x(:,1)=-33
  x(:,2)=32
  v=0
  rhoc_i1=0
  rhoc_i1(12,12,12,1,1,1)=1
  rhoc_i1(1,1,1,2,2,2)=1
  rhoc_i4=0
  np255=0
  nptile(1,1,1)=1
  nptile(2,2,2)=1

  rhoc(1:nt,1:nt,1:nt,:,:,:)=(rhoc_i1+int(-128,1))+128
  ntemp=count(rhoc(1:nt,1:nt,1:nt,:,:,:)==255)
  rhoc(1:nt,1:nt,1:nt,:,:,:)=unpack(rhoc_i4(:ntemp),rhoc(1:nt,1:nt,1:nt,:,:,:)==255,rhoc(1:nt,1:nt,1:nt,:,:,:))
  mass_p=884736
  v_i2r=1
#endif


a=a_i

nptotal=0
do i=1,nn**3
  nptotal=nptotal+nplocal[i]
enddo
mass_p=real((nf*nn)**3)/nptotal

if (head) then
	print*, 'nptotal =', nptotal
  print*, 'mass_p=', mass_p
	print*, 'max np per cell =', maxval(rhoc)
	!print*, 'min v_i =', minval(v)*v_i2r
	!print*, 'max v_i =', maxval(v)*v_i2r
	!print*, 'velocity resolution ='
	!print*, v_i2r
	print*, 'particle_initialization done'
	print*, ''
endif


contains

! for converting cubep3m zip checkpoint
function cumsum_cubep3m(input)
implicit none
integer(1) input(nc,nc,nc)
integer cumsum_cubep3m(nc,nc,nc)
integer nsum,igx,igy,igz
nsum=0
do igz=1,nc
do igy=1,nc
do igx=1,nc
  nsum=nsum+input(igx,igy,igz)
  cumsum_cubep3m(igx,igy,igz)=nsum
enddo
enddo
enddo
endfunction cumsum_cubep3m

endsubroutine
