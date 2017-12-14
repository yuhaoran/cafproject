! compile with
! gfortran lower_resolution_delta_L.f90
implicit none

integer,parameter :: ng=256 ! original resolution
integer,parameter :: n=128 ! lower resolution
integer i,j,k

real delta_L(ng,ng,ng)
real delta_new(n,n,n)

open(11,file='../output/universe3/image1/delta_L_1.bin',status='old',access='stream')
read(11) delta_L
close(11)

do k=1,n
do j=1,n
do i=1,n
  delta_new(i,j,k) = 0.125 * sum( delta_L(i*2-1:i*2, j*2-1:j*2, k*2-1:k*2) )
enddo
enddo
enddo

open(11,file='../output/universe4/image1/delta_L_1.bin',status='replace',access='stream')
write(11) delta_new
close(11)

end
