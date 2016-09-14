program init
use parameters

implicit none
integer, parameter :: npnode=(nf/2)**3
integer, parameter :: nptile=npnode/nnt**3
integer, parameter ::  nseedmax=200
integer  iseed(nseedmax), iseedsize,itile
logical head

real  x(3,npnode)[*], v(3,npnode)[*]
integer rank,mycoord3d(3)
integer i,j,k,l,ic(3),pp,hoc(nc,nc,nc),ll(npnode)
real, parameter:: v_resolution=1./32767.49
real v_i2r(3)[*],v_r2i(3)[*],vmax(3)[*]

integer(1) zip0_i1(3)
integer(2) zip1_i2(3)
integer(4) rhoc(nc,nc,nc)

rank=this_image()-1
mycoord3d(3)=rank/nn**2 ! 0,1,2,...
mycoord3d(2)=(rank-nn**2*mycoord3d(3))/nn
mycoord3d(1)=mod(rank,nn)

head=(this_image()==1)

call system('hostname')
print*, 'mycoord3d =', mycoord3d

call random_seed()
call random_seed(size=iseedsize)
if (iseedsize>nseedmax) stop 'need bigger nseedmax'
!print*, iseedsize
call random_seed(get=iseed(:iseedsize))
!print*, iseed(:iseedsize)
iseed=iseed+this_image()
call random_seed(put=iseed(:iseedsize))
x=0; v=0
call random_number(x) ! assign random locations and velocities
call random_number(v) 
x=(x+spread(mycoord3d,2,npnode))*nf
v=v-0.5

if (head) then
	print*, 'np on this node', size(x,2)
	print*, 'location min =', minval(x)
	print*, 'location max =', maxval(x)
	print*, 'velocity min =', minval(v)
	print*, 'velocity max =', maxval(v)
endif

! call link_list
hoc=0
pp=1
do
	if (pp>npnode) exit
	ic=floor(x(:,pp)/4)+1
	if (minval(ic)<1 .or. maxval(ic)>nc) stop 'particle out of bound'
	ll(pp)=hoc(ic(1),ic(2),ic(3))
	hoc(ic(1),ic(2),ic(3))=pp
	pp=pp+1
enddo

if (head) then
	print*, 'created linked list, pp-1 =', pp-1
endif

! call checkpoint
itile=0
do k=1,nnt
do j=1,nnt
do i=1,nnt
	itile=itile+1
	vmax(:,i,j,k)=maxval(abs(v(:,nptile*(itile-1)+1:nptile*itile)),2)
enddo
enddo
enddo

v_i2r=vmax*v_resolution
v_r2i=1./v_i2r

open(unit=10,file='zip0.dat',status='replace',access='stream')
open(unit=11,file='zip1.dat',status='replace',access='stream')
open(unit=12,file='zip2.dat',status='replace',access='stream')
open(unit=13,file='zip3.dat',status='replace',access='stream')

write(10) npnode,v_i2r
rhoc=0
do k=1,nc
do j=1,nc
do i=1,nc
	l=hoc(i,j,k)
	do while (l>0)
		rhoc(i,j,k)=rhoc(i,j,k)+1
		zip0_i1=int(mod(x(:,l)/4,1.)*256,kind=1)
		zip1_i2(1:3)=nint(v(1:3,l)*v_r2i(1:3),kind=2)
		write(10) zip0_i1
		write(11) zip1_i2
		l=ll(l)
	enddo
enddo
enddo
enddo
write(12) int(min(rhoc,255),kind=1)
write(13) pack(rhoc,rhoc>=255)
close(10); close(11); close(12); close(13)

if (head) then
	print*, 'checkpoint done'
	print*, 'max np per cell =', maxval(rhoc)
	print*, 'number of peaks (np/cell>=255) =', count(rhoc>=255)
	print*, 'max v_i =', vmax
	print*, 'velocity resolution =', v_i2r
endif

endprogram
