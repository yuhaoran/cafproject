use parameters
use penfft_config
implicit none

integer,parameter :: nmax_redshift=100

integer nplocal
integer,parameter :: npnode=nf**3 ! only true for this project
real,parameter :: density_buffer=1.2
integer,parameter :: npmax=npnode*density_buffer
integer i,j,k,l,np,nlast,ind,ip,dx,dxy,kg,mg,jg,ig,i0,j0,k0,pp,itx,ity,itz,idx,imove,nf_shake,ibin
integer i_s,j_s,k_s,ii_s,jj_s,kk_s,r
integer wp(npmax),pp_s

integer cur_checkpoint, n_checkpoint
integer rhoc(nt,nt,nt,nnt,nnt,nnt)

integer(izipx) x(3,npmax)
integer(2)   pid(4,npmax)
real xp(3,npmax),pos1(3),gpos(3),dpos(3),r2,r2min
real den(ng,ng,ng)

integer hoc_g(ng,ng,ng),ll_p(npmax) ! linked list
integer hoc_p(npmax),ll_g(ng**3) ! inverse linked list

real z_checkpoint(nmax_redshift)
character (200) :: fn0,fn2
logical search


print*, 'Voronoi density field interpolation'
print*, 'ng=',ng
print*, 'npmax=',npmax
do cur_checkpoint= n_checkpoint,n_checkpoint
  fn0='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip0_'//image2str(this_image()-1)//'.dat'
  fn2='.'//opath//'/node'//image2str(this_image()-1)//'/'//z2str(z_checkpoint(cur_checkpoint)) &
  //'zip2_'//image2str(this_image()-1)//'.dat'

  if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
  
  open(12,file=fn2,status='old',action='read',access='stream')
  read(12) sim
  ! check zip format and read rhoc
  if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
    print*, 'zip format incompatable'
    close(12)
    stop
  endif
  read(12) rhoc
  close(12)
  call print_header(sim)
  !print*, sum(rhoc), sim%nplocal
  nplocal=sim%nplocal
  if (head) print*, 'nplocal =',nplocal

  open(10,file=fn0,status='old',action='read',access='stream')
  read(10) x(:,:nplocal)
  close(10)

  ! convert zip into xv
  nlast=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
  if (head) print*, 'CIC interpolation on tile'
  do k=1,nt
  do j=1,nt
  do i=1,nt
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np
      ip=nlast+l
      pos1=nt*((/itx,ity,itz/)-1) + (/i,j,k/)-1 + ((x(:,ip)+ishift)+rshift)*x_resolution
      pos1=real(ng)*((/icx,icy,icz/)-1) + pos1*real(ng)/real(nc)
      xp(:,ip)=pos1 ! convert to real coordinates and record in xp
    enddo
    nlast=nlast+np
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  print*,ip,nplocal

  ! assign particles to grids
  do ip=1,nplocal
    i=floor(xp(1,ip))+1
    j=floor(xp(2,ip))+1
    k=floor(xp(3,ip))+1
    ll_p(ip)=hoc_g(i,j,k)
    hoc_g(i,j,k)=ip
  enddo

  ! assign non-empty grids to particles
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

  wp=1 ! weight -- number of grids assigned to the particle
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
            pos1=xp(:,pp)
            dpos=pos1-gpos
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
  den=0
  do ip=1,nplocal
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

  print*, 'nplocal =',nplocal
  print*, 'sum(den) =',sum(den*1.d0)
  print*, 'minval(den) =',minval(den)
  print*, 'convert to density contrast'
  den=den/(sum(den*1.d0)/ng/ng/ng)-1

  print*, 'mean(den) =',sum(den*1.d0)/ng/ng/ng
  print*, 'maxval(den) =',maxval(den)
  print*, 'minval(den) =',minval(den)
  
  open(11,file='.'//opath//'node'//image2str(this_image()-1)//'/delta_voronoi.dat',status='replace',access='stream')
  write(11) den
  close(11)
  print*, 'done'
enddo
end
