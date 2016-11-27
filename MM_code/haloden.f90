implicit none
integer,parameter :: np=223479 ! total number of particles
integer,parameter :: nf=4608 ! number of fine cells per dim
integer,parameter :: ng=128 ! number of grids per dim in output resolution

integer i,j,k,ip,ig,pp
integer i_s,j_s,k_s,ii_s,jj_s,kk_s,r
integer wp(np),pp_s
integer hoc_g(ng,ng,ng),ll_p(np) ! linked list
integer hoc_p(np),ll_g(ng**3) ! inverse linked list
real xp(3,np),den(ng,ng,ng),den_ngp(ng,ng,ng)
real r2,r2min,gpos(3),hpos(3),dpos(3)

logical search

! read in particle positions
open(11,file='xhalo.dat',status='old',access='stream')
read(11) xp
close(11)

xp=xp*ng/nf
wp=1 ! weight -- number of grids assigned to the particle

! assign halos to grids
do ip=1,np
  i=floor(xp(1,ip))+1
  j=floor(xp(2,ip))+1
  k=floor(xp(3,ip))+1
  den_ngp(i,j,k)=den(i,j,k)+1
  ll_p(ip)=hoc_g(i,j,k)
  hoc_g(i,j,k)=ip
enddo

! assign non-empty grids to halos
do k=1,ng
do j=1,ng
do i=1,ng
  ig=ng**2*(k-1)+ng*(j-1)+i
  pp=hoc_g(i,j,k)
  do
    if (pp==0) exit
    hoc_p(pp)=ig
    pp=ll_p(pp)
  enddo
enddo
enddo
enddo

! assign empty grids to halos
do k=1,ng
do j=1,ng
do i=1,ng
  gpos=(/i,j,k/)-0.5
  ig=ng**2*(k-1)+ng*(j-1)+i
  if (hoc_g(i,j,k)==0) then
    r=0
    search=.true.
    r2min=1000*ng**2
    do ! search larger layer
      if (.not.search) exit
      r=r+1
      do k_s=k-r,k+r
      do j_s=j-r,j+r
      do i_s=i-r,i+r
        kk_s=modulo(k_s-1,ng)+1
        jj_s=modulo(j_s-1,ng)+1
        ii_s=modulo(i_s-1,ng)+1
        pp=hoc_g(ii_s,jj_s,kk_s)
        do
          if (pp==0) exit ! find next cell
          search=.false. ! found particle
          hpos=xp(:,pp)
          dpos=hpos-gpos
          dpos=modulo(dpos+ng/2,real(ng))-ng/2
          r2=sum(dpos**2)
          if (r2<r2min) then
            r2min=r2
            pp_s=pp
          endif
          pp=ll_p(pp)
        enddo
      enddo
      enddo
      enddo
    enddo ! finished searching current layer
    ! found nearest particle, exit searching loop
    ll_g(ig)=hoc_p(pp_s) ! let current grid point to pp_s's chain
    hoc_p(pp_s)=ig ! replace pp_s's hoc to be current grid
    wp(pp_s)=wp(pp_s)+1
  endif
enddo
enddo
enddo

! assign density
do ip=1,np
  ig=hoc_p(ip)
  do
    if (ig==0) exit
    k=(ig-1)/ng**2+1
    j=(ig-1-ng**2*(k-1))/ng+1
    i=modulo(ig-1,ng)+1
    den(i,j,k)=den(i,j,k)+1./real(wp(ip));
    ig=ll_g(ig);
  enddo
enddo

print*, 'np =',np
print*, 'sum(den) =',sum(den*1.d0)
print*, 'minval(den) =',minval(den)

open(11,file='den_interp.dat',status='replace',access='stream')
write(11) den
close(11)
print*, 'wrote into file den_interp.dat'

end
