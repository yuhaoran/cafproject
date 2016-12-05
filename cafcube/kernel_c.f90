#define debug
!#define LRCKCORR

subroutine kernel_c
use variables
use penfft
implicit none
save
include 'fftw3.f'

character(*),parameter :: dir_kern='./kernels/'
integer hc1
integer il,ih,jl,jh,kl,kh,kx,ky,kz
integer i1,j1,k1
real rx,ry,rz,r,ck_table(3,4,4,4)
integer itemp(3)
real rtemp(2)

if (head) print*, 'coarse kernel initialization'

hc1=nc*nn/2+1
ck=0

! offsets for local volume

il=1+nc*(icx-1)
jl=1+nc*(icy-1)
kl=1+nc*(icz-1)

ih=nc*icx
jh=nc*icy
kh=nc*icz

! construct uncorrected force kernel

do k=kl,kh
  k1=k-kl+1
  if (k<nc/2+2) then
    rz=k-1
  else
    rz=k-1-nc
  endif
  rz=rz*ncell

  do j=jl,jh
    j1=j-jl+1
    if (j<nc/2+2) then
      ry=j-1
    else
      ry=j-1-nc
    endif
    ry=ry*ncell

    do i=il,ih
      i1=i-il+1
      if (i<nc/2+2) then
        rx=i-1
      else
        rx=i-1-nc
      endif
      rx=rx*ncell
      r=sqrt(rx**2+ry**2+rz**2)
      if (r==0.0) then
        ck(:,i1,j1,k1)=0.0
      else
        ck(1,i1,j1,k1)=-rx/r**3
        ck(2,i1,j1,k1)=-ry/r**3
        ck(3,i1,j1,k1)=-rz/r**3
      endif
    enddo

  enddo

enddo

! Read in kernel for 2-level matching

open(20,file=dir_kern//'wfxyzc.2.ascii',status='old')
do k=1,4
do j=1,4
do i=1,4
  read(20,'(3i4,3e16.8)') itemp(1:3),ck_table(:,i,j,k)
enddo
enddo
enddo
close(20)

! copy corrections to other octants of the kernel
! need nc > 8

if (nn==1) then

  ck(:,:4,:4,:4)=ck_table(:,:,:,:)
  do k=2,4
    ck(1:2,:4,:4,nc-k+2)=ck_table(1:2,:,:,k)
    ck(3,:4,:4,nc-k+2)=-ck_table(3,:,:,k)
  enddo
  do j=2,4
    ck(1,:4,nc-j+2,:4)=ck_table(1,:,j,:)
    ck(2,:4,nc-j+2,:4)=-ck_table(2,:,j,:)
    ck(3,:4,nc-j+2,:4)=ck_table(3,:,j,:)
  enddo
  do k=2,4
  do j=2,4
    ck(1,:4,nc-j+2,nc-k+2)=ck_table(1,:,j,k)
    ck(2:3,:4,nc-j+2,nc-k+2)=-ck_table(2:3,:,j,k)
  enddo
  enddo
  do i=2,4
    ck(1,nc-i+2,:4,:4)=-ck_table(1,i,:,:)
    ck(2:3,nc-i+2,:4,:4)=ck_table(2:3,i,:,:)
  enddo
  do k=2,4
  do i=2,4
    ck(1,nc-i+2,:4,nc-k+2)=-ck_table(1,i,:,k)
    ck(2,nc-i+2,:4,nc-k+2)=ck_table(2,i,:,k)
    ck(3,nc-i+2,:4,nc-k+2)=-ck_table(3,i,:,k)
  enddo
  enddo
  do j=2,4
  do i=2,4
    ck(1:2,nc-i+2,nc-j+2,:4)=-ck_table(1:2,i,j,:)
    ck(3,nc-i+2,nc-j+2,:4)=ck_table(3,i,j,:)
  enddo
  enddo
  do k=2,4
  do j=2,4
  do i=2,4
    ck(1:3,nc-i+2,nc-j+2,nc-k+2)=-ck_table(:,i,j,k)
  enddo
  enddo
  enddo

else

  if (icx==1 .and. icy==1 .and. icz==1) then
    ck(:,:4,:4,:4)=ck_table

  elseif (icx==1 .and. icy==1 .and. icz==nn) then
    do k=2,4
      ck(1:2,:4,:4,nc-k+2)=ck_table(1:2,:,:,k)
      ck(3,:4,:4,nc-k+2)=-ck_table(3,:,:,k)
    enddo
  elseif (icx==1 .and. icy==nn .and. icz==1) then
    do j=2,4
      ck(1,:4,nc-j+2,:4)=ck_table(1,:,j,:)
      ck(2,:4,nc-j+2,:4)=-ck_table(2,:,j,:)
      ck(3,:4,nc-j+2,:4)=ck_table(3,:,j,:)
    enddo
  elseif (icx==nn .and. icy==1 .and. icz==1) then
    do i=2,4
      ck(1,nc-i+2,:4,:4)=-ck_table(1,i,:,:)
      ck(2:3,nc-i+2,:4,:4)=ck_table(2:3,i,:,:)
    enddo


  elseif (icx==1 .and. icy==nn .and. icz==nn) then
    do k=2,4
    do j=2,4
      ck(1,:4,nc-j+2,nc-k+2)=ck_table(1,:,j,k)
      ck(2:3,:4,nc-j+2,nc-k+2)=-ck_table(2:3,:,j,k)
    enddo
    enddo
  elseif (icx==nn .and. icy==1 .and. icz==nn) then
    do k=2,4
    do i=2,4
      ck(1,nc-i+2,:4,nc-k+2)=-ck_table(1,i,:,k)
      ck(2,nc-i+2,:4,nc-k+2)=ck_table(2,i,:,k)
      ck(3,nc-i+2,:4,nc-k+2)=-ck_table(3,i,:,k)
    enddo
    enddo
  elseif (icx==nn .and. icy==nn .and. icz==1) then
    do j=2,4
    do i=2,4
      ck(1:2,nc-i+2,nc-j+2,:4)=-ck_table(1:2,i,j,:)
      ck(3,nc-i+2,nc-j+2,:4)=ck_table(3,i,j,:)
    enddo
    enddo

  elseif (icx==nn .and. icy==nn .and. icz==nn) then
    do k=2,4
    do j=2,4
    do i=2,4
      ck(1:3,nc-i+2,nc-j+2,nc-k+2)=-ck_table(:,i,j,k)
    enddo
    enddo
    enddo
  endif

endif
if (head) print*, 'finished octants'

#ifdef LRCKCORR

#else
r3=ck(1,:,:,:)
call fft_cube2pencil
call trans_zxy2xyz
!kern_c(1,:,:,:)=imag(cx)
do k=1,npen
do j=1,nc
do i=1,nc*nn/2+1
  ! or: kern_c(1,i,j,k)=rxlong(2*i,j,k)
  kern_c(1,i,j,k)=imag(cx(i,j,k))
enddo
enddo
enddo
! or: kern_c(1,:,:)=rxlong(2::2,:,:)
! or: kern_c(1,:,:)=imag(cx)

r3=ck(2,:,:,:)
call fft_cube2pencil
call trans_zxy2xyz
do k=1,npen 
do j=1,nc 
do i=1,nc*nn/2+1
  kern_c(2,i,j,k)=imag(cx(i,j,k))
enddo
enddo
enddo

r3=ck(3,:,:,:)
call fft_cube2pencil
call trans_zxy2xyz
do k=1,npen
do j=1,nc
do i=1,nc*nn/2+1
  kern_c(3,i,j,k)=imag(cx(i,j,k))
enddo
enddo
enddo
if (head) print*, 'kernel_c done'

#endif

!open(21,file='kern_c.dat',access='stream',status='replace')
!write(21) kern_c
!close(21)
!print*, 'sum of kern_c', sum(kern_c)

endsubroutine kernel_c
