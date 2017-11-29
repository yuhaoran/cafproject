! read qspace Jacobian matrix and compute the det and inv
!
program detJmatrix
  use parameters
  !use pencil_fft
  !use powerspectrum
  implicit none
  save
  integer i,j,k,i_dim,j_dim
  real cube33(3,3,ng,ng,ng),cube(ng,ng,ng),a3(3),a33(3,3),cubeinv(ng,ng,ng),rgrad(3,ng,ng,ng)

  call geometry
  print*, 'det Jmatrix'
  print*, 'ng=',ng

  print*, 'checkpoint at:'
  open(16,file='../main/redshifts.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    print*, z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  print*,''

  !call create_penfft_plan
  print*, 'ng=',ng

  do cur_checkpoint= 1,n_checkpoint
    cube=0; cube33=0
    print*,output_name('Jmatrix')
    open(10,file=output_name('Jmatrix'),status='old',access='stream')
    do i_dim=1,3
    do j_dim=1,3
      read(10) cube33(j_dim,i_dim,:,:,:)
    enddo
    enddo
    close(10)

    do k=1,ng
    do j=1,ng
    do i=1,ng
      a33=cube33(:,:,i,j,k)
      !determinant
      cube(i,j,k)=(a33(1,1)*a33(2,2)*a33(3,3) - a33(1,1)*a33(2,3)*a33(3,2)&
                - a33(1,2)*a33(2,1)*a33(3,3) + a33(1,2)*a33(2,3)*a33(3,1)&
                + a33(1,3)*a33(2,1)*a33(3,2) - a33(1,3)*a33(2,2)*a33(3,1))
      !inverse
      a33=matinv3(a33)
      cube33(:,:,i,j,k)=a33
      !det(inv)
      cubeinv(i,j,k)=(a33(1,1)*a33(2,2)*a33(3,3) - a33(1,1)*a33(2,3)*a33(3,2)&
                - a33(1,2)*a33(2,1)*a33(3,3) + a33(1,2)*a33(2,3)*a33(3,1)&
                + a33(1,3)*a33(2,1)*a33(3,2) - a33(1,3)*a33(2,2)*a33(3,1))
    enddo
    enddo
    enddo
    open(10,file=output_name('detJ'),status='replace',access='stream')
    write(10) cube
    close(10)
    open(10,file=output_name('invJ'),status='replace',access='stream')
    write(10) cube33
    close(10)
    open(10,file=output_name('detinvJ'),status='replace',access='stream')
    write(10) cubeinv
    close(10)

    open(10,file=output_name('gradphiq'),status='old',access='stream')
      read(10) rgrad(1,:,:,:)
      read(10) rgrad(2,:,:,:)
      read(10) rgrad(3,:,:,:)
    close(10)

    do k=1,ng
    do j=1,ng
    do i=1,ng
      a3=rgrad(:,i,j,k)
      a33=cube33(:,:,i,j,k)
      a3=matmul(a3,a33)
      rgrad(:,i,j,k)=a3
    enddo
    enddo
    enddo
    
    print*, output_name('gradphiqinvJ')
    open(10,file=output_name('gradphiqinvJ'),status='replace',access='stream')
      write(10) rgrad(1,:,:,:)
      write(10) rgrad(2,:,:,:)
      write(10) rgrad(3,:,:,:)
    close(10)

  enddo


  !call destroy_penfft_plan

contains

  pure function matinv3(A) result(B)
    real, intent(in) :: A(3,3)
    real             :: B(3,3)
    real             :: detinv
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function

end
