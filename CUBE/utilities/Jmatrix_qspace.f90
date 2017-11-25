! read qspace fields and compute the divergence
!
program Jmatrix_qspace
  use pencil_fft
  use powerspectrum
  implicit none
  save
  integer i_redshift,i_dim
  real cube3(3,ng,ng,ng),cube(ng,ng,ng),delta_L(ng,ng,ng)
  real xi(10,nbin)[*]

  call geometry
  print*, 'Jmatrix @ q-space'
  print*, 'checkpoint at:'
  open(16,file='../main/redshifts.txt',status='old')
  do i_redshift=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i_redshift)
    print*, z_checkpoint(i_redshift)
  enddo
  71 n_checkpoint=i_redshift-1
  close(16)
  print*,''

  call create_penfft_plan
  print*, 'ng=',ng

  do cur_checkpoint= 1,n_checkpoint
    cube=0; cube3=0
    print*,output_name('dsp')
    open(10,file=output_name('dsp'),status='old',access='stream')
    open(11,file=output_name('Jmatrix'),status='replace',access='stream')
    do i_dim=1,3
      print*,'analyzing dimension i=',i_dim
      read(10) cube(:,:,:)
      call get_gradient(cube,cube3)
      cube3(i_dim,:,:,:)=cube3(i_dim,:,:,:)+1
      write(11) cube3(1,:,:,:)
      write(11) cube3(2,:,:,:)
      write(11) cube3(3,:,:,:)
    enddo
    close(10)
    close(11)
  enddo
  call destroy_penfft_plan

contains

  subroutine get_gradient(scalarfield,rgrad)

    integer j_dim,dim_1,dim_2,dim_3
    integer i,j,k,kg,jg,ig
    real kx(3),kr
    real scalarfield(ng,ng,ng),rgrad(3,ng,ng,ng)
    complex cgrad(3,ng*nn/2+1,ng,npen)
    complex pdim, ekx(3)

    r3=scalarfield
    cgrad=0
    call pencil_fft_forward
    do k=1,npen
    do j=1,ng
    do i=1,ng*nn/2+1
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*ng+j
      ig=i
      kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
      kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
      ekx=exp(2*pi*(0,1)*kx/ng)
      do j_dim=1,3
        dim_1=j_dim
        dim_2=mod(dim_1,3)+1
        dim_3=mod(dim_2,3)+1
        pdim=(ekx(dim_1)-1)*(ekx(dim_2)+1)*(ekx(dim_3)+1)/4
        cgrad(j_dim,i,j,k)=cgrad(j_dim,i,j,k)+pdim*cxyz(i,j,k)
      enddo
    enddo
    enddo
    enddo

    do j_dim=1,3
      cxyz=cgrad(j_dim,:,:,:)
      call pencil_fft_backward
      rgrad(j_dim,:,:,:)=r3
    enddo
  endsubroutine

  pure function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real, intent(in) :: A(3,3)   !! Matrix
    real             :: B(3,3)   !! Inverse matrix
    real             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
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
