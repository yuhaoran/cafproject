! to be optimized - transpose, image1d, ctransfer more than one slab
program test_penfft
  use iso_fortran_env, only : int64
  implicit none

  integer(8) planx,plany,planz,iplanx,iplany,iplanz
  integer,parameter :: NULL=0
  real,parameter :: pi=3.14159

  integer,parameter :: nn=1
  integer,parameter :: nc=768 ! nc/image/dim, in physical volume, >=24
  integer,parameter :: npen=nc/nn ! nc /dim in shorter side of the pencil, for pencil decomposition

  !! MPI images !!
  integer rank,icx,icy,icz

  real temp_r,temp_theta
  integer i,j,k,l,seedsize
  integer(int64) time64
  integer,allocatable :: iseed(:)

  real        r3(nc,nc,nc),r0(nc,nc)
  complex     c3(nc/2,nc,nc)
  real        rxyz(nc*nn+2  ,nc,npen)
  complex     cxyz(nc*nn/2+1,nc,npen)
  complex     cyyyxz(npen,nn,nn,nc/2+1,npen)
  complex     cyyxz(nc,     nn,nc/2+1,npen)
  complex     czzzxy(npen,nn,nn,nc/2+1,npen)

  complex ctransfer1(nc/2,nc,nn)[*]
  complex ctransfer2(nc,nc/2+1,nn)[*]
  complex ctransfer3(npen,npen,nn,nn)[*]
  complex ctransfer4(nc/2+1,nc,nn)[*]

  equivalence(cyyyxz,cyyxz,r3,c3)
  equivalence(czzzxy,rxyz,cxyz)

  rank=this_image()-1            ! MPI_rank
  icz=rank/(nn**2)+1             ! image_z
  icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
  icx=mod(rank,nn)+1             ! image_x
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !print*, 'ilp64 test: default integer'
  !i=2147483647
  !i=i+1
  !print*,i
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (this_image()==1) print*, 'call fft plan'
  call create_penfft_plan

  !call system('hostname'); sync all

  if (this_image()==1) print*, 'initialize random number'
  call random_seed(size=seedsize)
  !print*,'min seedsize =', seedsize
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
  call system_clock(time64)
  do i=1,seedsize
    time64=time64*this_image()
    iseed(i)=lcg(time64)
    !print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
  enddo
  call random_seed(put=iseed)
  call random_number(r3)

  deallocate(iseed)

  if (this_image()==1) print*,'Box-Muller transform'
  do k=1,nc
  do j=1,nc
  do i=1,nc,2
    temp_theta=2*pi*r3(i,j,k)
    temp_r=sqrt(-2*log(1-r3(i+1,j,k)))
    r3(i,j,k)=temp_r*cos(temp_theta)
    r3(i+1,j,k)=temp_r*sin(temp_theta)
  enddo
  enddo
  enddo
  !r3=this_image()
  r0=r3(:,:,nc)
  print*, r3(1,1,1),r3(nc,nc,nc)!,nc*nc*nc

  sync all

  if (this_image()==1) print*, 'call ftran'
  call fft_cube2pencil
  if (this_image()==1) print*, 'transpose to xyz'
  call trans_zxy2xyz
  if (this_image()==1) print*, 'transpose to zxy'
  call trans_xyz2zxy
  if (this_image()==1) print*, 'call btran'
  call ifft_pencil2cube
  if (this_image()==1) print*, 'destroy plan'
  call destroy_penfft_plan
  sync all
  print*, 'precision on image',this_image(),maxval(abs(r3(:,:,nc)-r0))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    !new cube->x loop
    do l=1,npen ! loop over cells in z, extract slabs
      ctransfer1(:,:,1:nn)=c3(:,:,l::npen) ! nn slabs of c3 copied to ctransfer1
      sync all
      do i1=1,nn ! loop over parts in x, get slabs from each y node
        ! i1=mod()
        cxyz(nc*(i1-1)/2+1:nc*i1/2,:,l)=ctransfer1(:,:,m2)[image1d(i1,m1,m3)]
      enddo
      sync all
    enddo
!if (this_image()==8) print*, cxyz
!stop

    ! cxyz(nc*nn/2+1,:,:) are left zeros.
    call sfftw_execute(planx)

    !new x->y loop
    do l=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer2(:,:,i1)=transpose(cxyz(nc/2*(i1-1)+1:nc/2*i1+1,:,l))
      enddo
      sync all
      do i1=1,nn
        cyyxz(:,i1,:,l)=ctransfer2(:,:,m1)[image1d(i1,m2,m3)]
      enddo
      sync all
    enddo
    ! for m1/=nn, there are extra elements in x+
!if (this_image()==8) print*, cyyxz
!stop

    call sfftw_execute(plany)

    !new y->z loop
    do l=1,nc/2+1 ! loop over slices in x direction
      do i2=1,nn
      do i1=1,nn
        ctransfer3(:,:,i1,i2)=transpose(cyyyxz(:,i1,i2,l,:))
      enddo
      enddo
      sync all
      do i2=1,nn
      do i1=1,nn
        czzzxy(:,i1,i2,l,:)=ctransfer3(:,:,m2,m3)[image1d(m1,i1,i2)]
      enddo
      enddo
      sync all
    enddo
    ! for m1/=nn, there are extra elements in x+
!if (this_image()==8) print*, czzzxy
!stop

    call sfftw_execute(planz)

  endsubroutine fft_cube2pencil

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    !new z->y loop
!if (this_image()==7) print*, czzzxy
!stop
    do l=1,nc/2+1 ! loop over slices in x direction
      do i2=1,nn
      do i1=1,nn
        ctransfer3(:,:,i1,i2)=transpose(czzzxy(:,i1,i2,l,:))
      enddo
      enddo
      sync all
      do i2=1,nn
      do i1=1,nn
        cyyyxz(:,i1,i2,l,:)=ctransfer3(:,:,m2,m3)[image1d(m1,i1,i2)]
      enddo
      enddo
      sync all
    enddo
!if (this_image()==1) print*, cyyxz
!stop
    call sfftw_execute(iplany)

    !new y->x loop
    do l=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer4(:,:,i1)=transpose(cyyxz(:,i1,:,l))
      enddo
      sync all
      do i1=1,nn
        cxyz(nc/2*(i1-1)+1:nc/2*i1+1,:,l)=ctransfer4(:,:,m1)[image1d(i1,m2,m3)]
      enddo
      sync all
    enddo
!if (this_image()==8) print*, cxyz
!stop
    call sfftw_execute(iplanx)

    !new x->cube loop
!sync all
    do l=1,npen
      do i1=1,nn
        ctransfer1(:,:,i1)=cxyz(nc*(i1-1)/2+1:nc*i1/2,:,l)
      enddo
!if (this_image()==5) print*,ctransfer1
!stop
      sync all
      do i1=1,nn
        c3(:,:,l+(i1-1)*npen)=ctransfer1(:,:,m1)[image1d(m2,i1,m3)]
      enddo
      sync all
    enddo
!if (this_image()==8) print*, c3
!stop

    r3=r3/(nc*nn)/(nc*nn)/(nc*nn)
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

  endsubroutine trans_xyz2zxy

  subroutine create_penfft_plan
  !  use variables
    implicit none
    save
    include 'fftw3.f'
    call sfftw_plan_many_dft_r2c(planx,1,nc*nn,nc*npen,cxyz,NULL,1,nc*nn+2,cxyz,NULL,1,nc*nn/2+1,FFTW_MEASURE)
    call sfftw_plan_many_dft_c2r(iplanx,1,nc*nn,nc*npen,cxyz,NULL,1,nc*nn/2+1,cxyz,NULL,1,nc*nn+2,FFTW_MEASURE)
    call sfftw_plan_many_dft(plany,1,nc*nn,(nc/2+1)*npen,cyyxz,NULL,1,nc*nn,cyyxz,NULL,1,nc*nn,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplany,1,nc*nn,(nc/2+1)*npen,cyyxz,NULL,1,nc*nn,cyyxz,NULL,1,nc*nn,FFTW_BACKWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(planz,1,nc*nn,(nc/2+1)*npen,czzzxy,NULL,1,nc*nn,czzzxy,NULL,1,nc*nn,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplanz,1,nc*nn,(nc/2+1)*npen,czzzxy,NULL,1,nc*nn,czzzxy,NULL,1,nc*nn,FFTW_BACKWARD,FFTW_MEASURE)
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

  function image1d(cx,cy,cz)
  integer image1d,cx,cy,cz
  image1d=cx+nn*(cy-1)+nn**2*(cz-1)
  endfunction

  function lcg(s) !// Linear congruential generator
    implicit none
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  endfunction lcg
end
