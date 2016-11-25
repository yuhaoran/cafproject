integer,parameter :: ng=128 ! number of grid per dim
integer,parameter :: nplot=4 ! thickness of the column

integer i,j,k

real den(ng,ng,ng)
real def(ng,ng,ng)
real dspz(ng,ng,ng)
real proj(ng,ng)

! read in density (for delta, please +1)
open(11,file='den.dat',status='old',access='stream')
read(11) den
close(11)

! read in deformation potential
open(11,file='def.dat',status='old',access='stream')
read(11) def
close(11)

! get displacement field
do k=1,ng
do j=1,ng
do i=1,ng
  kp=mod(k,ng)+1
  jp=mod(j,ng)+1
  ip=mod(i,ng)+1
  dspz(i,j,k)=def(i,j,kp)+def(ip,j,kp)+def(i,jp,kp)+def(ip,jp,kp)-def(i,j,k)-def(ip,j,k)-def(i,jp,k)-def(ip,jp,k)
enddo
enddo
enddo
dspz=dspz/4

! get corrected column density
k1=ng/2
k2=k1+nplot-1
do j=1,ng
do i=1,ng
  jm=mod(j-2,ng)+1
  im=mod(i-2,ng)+1
  k1m=mod(k1-2,ng)+1
  kk1=k1-1+(dsp(3,im,jm,k1m)+dsp(3,i,jm,k1m)+dsp(3,im,j,k1m)+dsp(3,i,j,k1m))/4
  kk2=k2+(dsp(3,im,jm,k2)+dsp(3,i,jm,k2)+dsp(3,im,j,k2)+dsp(3,i,j,k2))/4
  nk1=floor(kk1)+1
  nk2=floor(kk2)+1
  r1=nk1-kk1
  r2=kk2-nk2+1
  proj(i,j)=r1*den(i,j,nk1)+sum(den(i,j,nk1+1:nk2-1))+r2*den(i,j,nk2)-sum(den(i,j,nk2:nk1))
enddo 
enddo
proj=proj/nplot

open(11,file='proj.dat',status='replace',access='stream')
write(11) proj
close(11)

print*, proj(1,1), proj(ng,ng)

end
