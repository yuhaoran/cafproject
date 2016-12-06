module penfft_fine
use penfft_config
implicit none
save

real,parameter :: pi=2*asin(1.)
integer,parameter :: nk_ny=ng*nn/2 ! Nyquist wave number
integer,parameter :: nbin=floor(4*log(nk_ny*sqrt(3.)/0.95)/log(2.))

! FFT plans
integer(8) planx_fine,plany_fine,planz_fine,iplanx_fine,iplany_fine,iplanz_fine

! pencil fft arrays
real cube(ng,ng,ng)[*]
integer,parameter :: NULL=0
real rxlong(ng*nn+2,ng,ngpen)[*]
complex cx(ng*nn/2+1,ng,ngpen)[*]
complex cy(ngpen,nn,nn,ng/2+1,ngpen)[*]
complex cz(ngpen,nn,nn,ng/2+1,ngpen)[*]
complex cyxz(ngpen*nn,ngpen*nn/2+1,ngpen)
complex cyyxz(ngpen,nn,ngpen*nn/2+1,ngpen)
equivalence(cyxz,cyyxz)



contains

subroutine create_penfft_fine_plan
  implicit none
  save
  include 'fftw3.f'
  if (head) then
    print*, 'create fftw plans for penfft'
    print*, 'ng =', ng
  endif
  call sfftw_plan_many_dft_r2c(planx_fine,1,ng*nn,ng*ngpen,rxlong,NULL,1,ng*nn+2,rxlong,NULL,1,ng*nn/2+1,FFTW_MEASURE)
  call sfftw_plan_many_dft_c2r(iplanx_fine,1,ng*nn,ng*ngpen,rxlong,NULL,1,ng*nn/2+1,rxlong,NULL,1,ng*nn+2,FFTW_MEASURE)
  call sfftw_plan_many_dft(plany_fine,1,ng*nn,(ng/2+1)*ngpen,cy,NULL,1,ng*nn,cy,NULL,1,ng*nn,FFTW_FORWARD,FFTW_MEASURE)
  call sfftw_plan_many_dft(iplany_fine,1,ng*nn,(ng/2+1)*ngpen,cy,NULL,1,ng*nn,cy,NULL,1,ng*nn,FFTW_BACKWARD,FFTW_MEASURE)
  call sfftw_plan_many_dft(planz_fine,1,ng*nn,(ng/2+1)*ngpen,cz,NULL,1,ng*nn,cz,NULL,1,ng*nn,FFTW_FORWARD,FFTW_MEASURE)
  call sfftw_plan_many_dft(iplanz_fine,1,ng*nn,(ng/2+1)*ngpen,cz,NULL,1,ng*nn,cz,NULL,1,ng*nn,FFTW_BACKWARD,FFTW_MEASURE)
endsubroutine

subroutine destroy_penfft_fine_plan
  implicit none
  save
  include 'fftw3.f'
  call sfftw_destroy_plan(planx_fine)
  call sfftw_destroy_plan(iplanx_fine)
  call sfftw_destroy_plan(plany_fine)
  call sfftw_destroy_plan(iplany_fine)
  call sfftw_destroy_plan(planz_fine)
  call sfftw_destroy_plan(iplanz_fine)
endsubroutine


subroutine fft_cube2pencil_fine
implicit none
save
integer m,m1,m2,m3,i0,i1,i2

m1=icx
m2=icy
m3=icz

m=num_images()

do i0=1,nn ! cube->x
  i1=mod(m2+i0-2,nn)+1
  rxlong(ngpen*nn*(i1-1)+1:ngpen*nn*i1,:,:)=cube(:,:,ngpen*(m2-1)+1:ngpen*m2)[image1d(i1,m1,m3)]
enddo

sync all

call sfftw_execute(planx_fine)

sync all

cx=cmplx(rxlong(::2,:,:),rxlong(2::2,:,:))

sync all

do i0=1,nn ! x->y
  i1=mod(m1+i0-2,nn)+1
  cyxz=trans12(cx(ngpen*nn/2*(m1-1)+1:ngpen*nn/2*m1+1,:,:)[image1d(i1,m2,m3)],ngpen*nn/2+1,ngpen*nn,ngpen)
  ! for m1/=nn, two more redundant layers
  cy(:,:,i1, :, :)=cyyxz
enddo

sync all

call sfftw_execute(plany_fine)

sync all


do i0=1,nn**2 ! y->z
  i1=mod(i0+m-2,nn)+1
  i2=mod((i0+m-2)/nn,nn)+1
  cz(:,i1,i2,:,:)=trans13(cy(:,m2,m3, :, :)[image1d(m1,i1,i2)],ngpen,ngpen*nn/2+1,ngpen)
enddo

sync all

call sfftw_execute(planz_fine)

sync all

endsubroutine fft_cube2pencil_fine

subroutine ifft_pencil2cube_fine
implicit none
save
integer m,m1,m2,m3,i0,i1,i2

m1=icx
m2=icy
m3=icz

m=num_images()

call sfftw_execute(iplanz_fine)

sync all

do i0=1,nn**2
  i1=mod(i0+m-2,nn)+1
  i2=mod((i0+m-2)/nn,nn)+1
  cy(:,i1,i2,:,:)=trans13(cz(:,m2,m3, :, :)[image1d(m1,i1,i2)],ngpen,ngpen*nn/2+1,ngpen)
enddo

sync all

call sfftw_execute(iplany_fine)

sync all

do i0=1,nn
  i1=mod(m1+i0-2,nn)+1
  cyyxz=cy(:,:,m1, :, :)[image1d(i1,m2,m3)]
  cx(ngpen*nn/2*(i1-1)+1:ngpen*nn/2*i1+1,:,:)=trans12(cyxz,ngpen*nn,ngpen*nn/2+1,ngpen)
enddo

sync all

rxlong(::2,:,:)=real(cx); rxlong(2::2,:,:)=imag(cx)

call sfftw_execute(iplanx_fine)

do i0=1,nn
  i1=mod(m2+i0-2,nn)+1
  cube(:,:,ngpen*(i1-1)+1:ngpen*i1)=rxlong(ngpen*nn*(m1-1)+1:ngpen*nn*m1,:,:)[image1d(m2,i1,m3)]
enddo

cube=cube/(ng*nn)**3

sync all

endsubroutine ifft_pencil2cube_fine

subroutine trans_zxy2xyz_fine
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
  cy(:,i1,i2,:,:)=trans13(cz(:,m2,m3, :, :)[image1d(m1,i1,i2)],ngpen,ngpen*nn/2+1,ngpen)
enddo

sync all

do i0=1,nn
  i1=mod(m1+i0-2,nn)+1
  cyyxz=cy(:,:,m1, :, :)[image1d(i1,m2,m3)]
  cx(ngpen*nn/2*(i1-1)+1:ngpen*nn/2*i1+1,:,:)=trans12(cyxz,ngpen*nn,ngpen*nn/2+1,ngpen)
enddo

rxlong(::2,:,:)=real(cx); rxlong(2::2,:,:)=imag(cx)

sync all
endsubroutine trans_zxy2xyz_fine

subroutine trans_xyz2zxy_fine
implicit none
save
integer m,m1,m2,m3,i0,i1,i2

m1=icx
m2=icy
m3=icz

m=num_images()

do i0=1,nn
  i1=mod(m1+i0-2,nn)+1
  cyxz=trans12(cx(ngpen*nn/2*(m1-1)+1:ngpen*nn/2*m1+1,:,:)[image1d(i1,m2,m3)],ngpen*nn/2+1,ngpen*nn,ngpen)
  ! for m1/=nn, two more redundant layers
  cy(:,:,i1, :, :)=cyyxz
enddo

sync all

do i0=1,nn**2
  i1=mod(i0+m-2,nn)+1
  i2=mod((i0+m-2)/nn,nn)+1
  cz(:,i1,i2,:,:)=trans13(cy(:,m2,m3, :, :)[image1d(m1,i1,i2)],ngpen,ngpen*nn/2+1,ngpen)
enddo

sync all

endsubroutine trans_xyz2zxy_fine

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


endmodule 
