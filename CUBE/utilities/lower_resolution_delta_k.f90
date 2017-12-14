implicit none

integer,parameter :: ng=256 ! original resolution
integer,parameter :: n=128 ! lower resolution

complex cxyz(ng/2+1,ng,ng)
complex cnew(n/2+1,n,n)

open(11,file='../output/universe1/image1/delta_k_1.bin',status='old',access='stream')
read(11) cxyz
close(11)



open(11,file='../output/universe2/image1/delta_k_1.bin',status='old',access='stream')
read(11) cnew
close(11)

end
