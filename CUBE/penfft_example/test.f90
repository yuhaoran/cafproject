implicit none

integer(8) planx,plany,planz,iplanx,iplany,iplanz
integer,parameter :: NULL=0

integer,parameter :: nn=1
integer,parameter :: nc=512 ! nc/image/dim, in physical volume, >=24
integer,parameter :: npen=nc/nn ! nc /dim in shorter side of the pencil, for pencil decomposition

!! MPI images !!
integer,parameter :: rank=0                         ! MPI_rank
integer,parameter :: icz=rank/(nn**2)+1             ! image_z
integer,parameter :: icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
integer,parameter :: icx=mod(rank,nn)+1             ! image_x

real r3(nc,nc,nc)[*]
real rxlong(nc*nn+2,nc,npen)[*]
complex cx(nc*nn/2+1,nc,npen)[*]
real crho_c(nc*nn+2,nc,npen)
complex cy(npen,nn,nn,nc/2+1,npen)[*]
complex cz(npen,nn,nn,nc/2+1,npen)[*]
complex cyxz(npen*nn,npen*nn/2+1,npen)
complex cyyxz(npen,nn,npen*nn/2+1,npen)
equivalence(cyxz,cyyxz)

print*, 'call fft plan'
call create_penfft_plan
print*, 'initialize den'
r3=1.23
print*, 'call ftran'
call fft_cube2pencil
print*, 'transpose to xyz'
!call trans_zxy2xyz
print*, 'transpose to zxy'
!call trans_xyz2zxy
print*, 'call btran'
r3=0
!call ifft_pencil2cube

print*, r3(1:4,1,1)

call destroy_penfft_plan


contains

subroutine fft_cube2pencil
!use variables
implicit none
save
integer m,m1,m2,m3,i0,i1,i2

m1=icx
m2=icy
m3=icz

m=num_images()

do i0=1,nn ! cube->x
  i1=mod(m2+i0-2,nn)+1
  !print*, i1,m1,m3,image1d(i1,m1,m3)
  rxlong(npen*nn*(i1-1)+1:npen*nn*i1,:,:)=r3(:,:,npen*(m2-1)+1:npen*m2)[image1d(i1,m1,m3)]
enddo
!print*, 'got here'
sync all

call sfftw_execute(planx)
sync all

!print*,'called x'

cx=cmplx(rxlong(::2,:,:),rxlong(2::2,:,:))

sync all

do i0=1,nn ! x->y
  i1=mod(m1+i0-2,nn)+1
  !cyxz=trans12(cx(npen*nn/2*(m1-1)+1:npen*nn/2*m1+1,:,:)[image1d(i1,m2,m3)],npen*nn/2+1,npen*nn,npen)
  ! for m1/=nn, two more redundant layers
  cy(:,:,i1, :, :)=cyyxz
enddo

sync all

call sfftw_execute(plany)

sync all

!print*, 'called y'

do i0=1,nn**2 ! y->z
  i1=mod(i0+m-2,nn)+1
  i2=mod((i0+m-2)/nn,nn)+1
!print*, i1,i2,image1d(m1,i1,i2)
  cz(:,i1,i2,:,:)=trans13(cy(:,m2,m3, :, :)[image1d(m1,i1,i2)],npen,npen*nn/2+1,npen)
!print*,'here'
enddo


sync all

call sfftw_execute(planz)

!print*, 'called z'

sync all
endsubroutine fft_cube2pencil






subroutine ifft_pencil2cube
!use variables
implicit none
save
integer m,m1,m2,m3,i0,i1,i2

m1=icx
m2=icy
m3=icz

m=num_images()

call sfftw_execute(iplanz)

sync all

do i0=1,nn**2
  i1=mod(i0+m-2,nn)+1
  i2=mod((i0+m-2)/nn,nn)+1
  cy(:,i1,i2,:,:)=trans13(cz(:,m2,m3, :, :)[image1d(m1,i1,i2)],npen,npen*nn/2+1,npen)
enddo

sync all

call sfftw_execute(iplany)

sync all

do i0=1,nn
  i1=mod(m1+i0-2,nn)+1
  cyyxz=cy(:,:,m1, :, :)[image1d(i1,m2,m3)]
  cx(npen*nn/2*(i1-1)+1:npen*nn/2*i1+1,:,:)=trans12(cyxz,npen*nn,npen*nn/2+1,npen)
enddo

sync all

rxlong(::2,:,:)=real(cx); rxlong(2::2,:,:)=imag(cx)

call sfftw_execute(iplanx)

do i0=1,nn
  i1=mod(m2+i0-2,nn)+1
  r3(:,:,npen*(i1-1)+1:npen*i1)=rxlong(npen*nn*(m1-1)+1:npen*nn*m1,:,:)[image1d(m2,i1,m3)]
enddo

sync all

r3=r3/(nc*nn)**3

endsubroutine ifft_pencil2cube




subroutine trans_zxy2xyz
!use variables
implicit none
save
integer m,m1,m2,m3,i0,i1,i2

m1=icx
m2=icy
m3=icz

m=num_images()

do i0=1,nn**2
  i1=mod(i0+m-2,nn)+1
  i2=mod((i0+m-2)/nn,nn)+1
  cy(:,i1,i2,:,:)=trans13(cz(:,m2,m3, :, :)[image1d(m1,i1,i2)],npen,npen*nn/2+1,npen)
enddo

sync all

do i0=1,nn
  i1=mod(m1+i0-2,nn)+1
  cyyxz=cy(:,:,m1, :, :)[image1d(i1,m2,m3)]
  cx(npen*nn/2*(i1-1)+1:npen*nn/2*i1+1,:,:)=trans12(cyxz,npen*nn,npen*nn/2+1,npen)
enddo

rxlong(::2,:,:)=real(cx); rxlong(2::2,:,:)=imag(cx)

endsubroutine trans_zxy2xyz

subroutine trans_xyz2zxy
!use variables
implicit none
save
integer m,m1,m2,m3,i0,i1,i2

m1=icx
m2=icy
m3=icz

m=num_images()

do i0=1,nn
  i1=mod(m1+i0-2,nn)+1
  cyxz=trans12(cx(npen*nn/2*(m1-1)+1:npen*nn/2*m1+1,:,:)[image1d(i1,m2,m3)],npen*nn/2+1,npen*nn,npen)
  ! for m1/=nn, two more redundant layers
  cy(:,:,i1, :, :)=cyyxz
enddo

sync all

do i0=1,nn**2
  i1=mod(i0+m-2,nn)+1
  i2=mod((i0+m-2)/nn,nn)+1
  cz(:,i1,i2,:,:)=trans13(cy(:,m2,m3, :, :)[image1d(m1,i1,i2)],npen,npen*nn/2+1,npen)
enddo

sync all

endsubroutine trans_xyz2zxy

subroutine create_penfft_plan
!  use variables
  implicit none
  save
  include 'fftw3.f'
  call sfftw_plan_many_dft_r2c(planx,1,nc*nn,nc*npen,rxlong,NULL,1,nc*nn+2,rxlong,NULL,1,nc*nn/2+1,FFTW_MEASURE)
  call sfftw_plan_many_dft_c2r(iplanx,1,nc*nn,nc*npen,rxlong,NULL,1,nc*nn/2+1,rxlong,NULL,1,nc*nn+2,FFTW_MEASURE)
  call sfftw_plan_many_dft(plany,1,nc*nn,(nc/2+1)*npen,cy,NULL,1,nc*nn,cy,NULL,1,nc*nn,FFTW_FORWARD,FFTW_MEASURE)
  call sfftw_plan_many_dft(iplany,1,nc*nn,(nc/2+1)*npen,cy,NULL,1,nc*nn,cy,NULL,1,nc*nn,FFTW_BACKWARD,FFTW_MEASURE)
  call sfftw_plan_many_dft(planz,1,nc*nn,(nc/2+1)*npen,cz,NULL,1,nc*nn,cz,NULL,1,nc*nn,FFTW_FORWARD,FFTW_MEASURE)
  call sfftw_plan_many_dft(iplanz,1,nc*nn,(nc/2+1)*npen,cz,NULL,1,nc*nn,cz,NULL,1,nc*nn,FFTW_BACKWARD,FFTW_MEASURE)
endsubroutine create_penfft_plan

subroutine destroy_penfft_plan
!  use variables
  implicit none
  save
  include 'fftw3.f'
  call sfftw_destroy_plan(planx)
  call sfftw_destroy_plan(iplanx)
  call sfftw_destroy_plan(plany)
  call sfftw_destroy_plan(iplany)
  call sfftw_destroy_plan(planz)
  call sfftw_destroy_plan(iplanz)
endsubroutine destroy_penfft_plan

pure function trans12(cmatrix,m,n,h)
  implicit none
  integer(4) :: i
  integer(4), intent(in) :: m,n,h
  complex(4), intent(in) :: cmatrix(m,n,h)
  complex(4) :: trans12(n,m,h)
  do i=1,h
    trans12(:,:,i)=transpose(cmatrix(:,:,i))
  enddo
endfunction

pure function trans13(cmatrix,m,h,n)
implicit none
integer(4) :: i
integer(4), intent(in) :: m,n,h
complex(4), intent(in) :: cmatrix(m,h,n)
complex(4) :: trans13(n,h,m)
do i=1,h
  trans13(:,i,:)=transpose(cmatrix(:,i,:))
enddo
endfunction

function image1d(cx,cy,cz)
integer image1d,cx,cy,cz
image1d=cx+nn*(cy-1)+nn**2*(cz-1)
endfunction

!pure function trans23(cmatrix,h,m,n)
!implicit none
!integer(4) :: i
!integer(4), intent(in) :: m,n,h
!complex(4), intent(in) :: cmatrix(h,m,n)
!complex(4) :: trans23(h,n,m)
!do i=1,h
!  trans23(i,:,:)=transpose(cmatrix(i,:,:))
!enddo
!endfunction

end
