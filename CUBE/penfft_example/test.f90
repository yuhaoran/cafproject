! to be optimized - check algorithm, transpose, image1d,
! equivalence, ctransfer more than one slab, parallel
program test
  use iso_fortran_env, only : int64
  implicit none

  integer(8) planx,plany,planz,iplanx,iplany,iplanz
  integer,parameter :: NULL=0
  real,parameter :: pi=3.14159

  integer,parameter :: nn=1
  integer,parameter :: nc=128 ! nc/image/dim, in physical volume, >=24
  integer,parameter :: npen=nc/nn ! nc /dim in shorter side of the pencil, for pencil decomposition

  !! MPI images !!
  integer,parameter :: rank=0                         ! MPI_rank
  integer,parameter :: icz=rank/(nn**2)+1             ! image_z
  integer,parameter :: icy=(rank-nn**2*(icz-1))/nn+1  ! image_y
  integer,parameter :: icx=mod(rank,nn)+1             ! image_x
  real temp_r,temp_theta
  integer i,j,k,l,seedsize
  integer(int64) time64
  integer,allocatable :: iseed(:)

  real        r3(nc,nc,nc)
  complex     c3(nc/2,nc,nc)
  equivalence(r3,c3)

  real        rx(nc*nn+2  ,nc,npen)
  complex     cx(nc*nn/2+1,nc,npen)
  equivalence(rx,cx)

  complex     cyyyxz(npen,nn,nn,nc/2+1,npen)
  complex     cyyxz(nc,     nn,nc/2+1,npen)
  equivalence(cyyyxz,cyyxz)

  complex     cz(npen,nn,nn,nc/2+1,npen)

  complex ctransfer1(nc/2,nc,nn)[*]
  complex ctransfer2(nc,nc/2+1,nn)[*]
  complex ctransfer3(npen,npen,nn,nn)[*]
  complex ctransfer4(nc/2+1,nc,nn)[*]

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*, 'ilp64 test: default integer'
  i=2147483647
  i=i+1
  print*,i

  print*, 'call fft plan'
  call create_penfft_plan


  print*, 'initialize r3'
  call random_seed(size=seedsize)
  !print*,'min seedsize =', seedsize
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
  call system_clock(time64)
  do i=1,seedsize
    iseed(i)=lcg(time64)
    !print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
  enddo
  call random_seed(put=iseed)
  call random_number(r3)

  deallocate(iseed)

  print*,'Box-Muller transform'
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
  print*, r3(1:10,1,1)


  print*, 'call ftran'
  call fft_cube2pencil
  print*, 'transpose to xyz'
  call trans_zxy2xyz
  print*, 'transpose to zxy'
  call trans_xyz2zxy
  print*, 'call btran'
  call ifft_pencil2cube

  print*, r3(1:10,1,1)
  call destroy_penfft_plan

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
      do i1=1,nn ! loop over parts in x, get slabs from each y node
        ! i1=mod()
        cx(nc*(i1-1)/2+1:nc*i1/2,:,l)=ctransfer1(:,:,m2)[image1d(i1,m1,m3)]
      enddo
    enddo
    ! cx(nc*nn/2+1,:,:) are left zeros.

    call sfftw_execute(planx)

    !new x->y loop
    do l=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer2(:,:,i1)=transpose(cx(nc/2*(i1-1)+1:nc/2*i1+1,:,l))
      enddo
      sync all
      do i1=1,nn
        cyyxz(:,i1,:,l)=ctransfer2(:,:,m1)[image1d(i1,m2,m3)]
      enddo
    enddo

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
        cz(:,i1,i2,l,:)=ctransfer3(:,:,m2,m3)[image1d(m1,i1,i2)]
      enddo
      enddo
    enddo

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
    do l=1,nc/2+1 ! loop over slices in x direction
      do i2=1,nn
      do i1=1,nn
        ctransfer3(:,:,i1,i2)=transpose(cz(:,i1,i2,l,:))
      enddo
      enddo
      sync all
      do i2=1,nn
      do i1=1,nn
        cyyyxz(:,i1,i2,l,:)=ctransfer3(:,:,m2,m3)[image1d(m1,i1,i2)]
      enddo
      enddo
    enddo

    call sfftw_execute(iplany)

    !new y->x loop
    do l=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer4(:,:,i1)=transpose(cyyxz(:,i1,:,l))
      enddo
      sync all
      do i1=1,nn
        cx(nc/2*(i1-1)+1:nc/2*i1+1,:,l)=ctransfer4(:,:,m1)[image1d(i1,m2,m3)]
      enddo
    enddo

    call sfftw_execute(iplanx)

    !new x->cube loop
    do l=1,npen
      do i1=1,nn
        ctransfer1(:,:,i1)=cx(nc*(i1-1)/2+1:nc*i1/2,:,l)
      enddo
      sync all
      do i1=1,nn
        c3(:,:,l+(i1-1)*npen)=ctransfer1(:,:,m2)[image1d(i1,m1,m3)]
      enddo
    enddo

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
    call sfftw_plan_many_dft_r2c(planx,1,nc*nn,nc*npen,cx,NULL,1,nc*nn+2,cx,NULL,1,nc*nn/2+1,FFTW_MEASURE)
    call sfftw_plan_many_dft_c2r(iplanx,1,nc*nn,nc*npen,cx,NULL,1,nc*nn/2+1,cx,NULL,1,nc*nn+2,FFTW_MEASURE)
    call sfftw_plan_many_dft(plany,1,nc*nn,(nc/2+1)*npen,cyyxz,NULL,1,nc*nn,cyyxz,NULL,1,nc*nn,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplany,1,nc*nn,(nc/2+1)*npen,cyyxz,NULL,1,nc*nn,cyyxz,NULL,1,nc*nn,FFTW_BACKWARD,FFTW_MEASURE)
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
