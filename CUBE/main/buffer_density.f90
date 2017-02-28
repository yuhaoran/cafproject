!#define debug
!#define redundant
subroutine buffer_density
use variables
implicit none
save
!integer ntemp
integer nshift,nlen,nlast,ifrom

if (head) print*, 'buffer_density'

! sync buffer regions in rhoc for each tile
!x
rhoc(:0,:,1:nt,1,:,:)=rhoc(nt-ncb+1:nt,:,1:nt,nnt,:,:)[image1d(inx,icy,icz)]
rhoc(:0,:,1:nt,2:,:,:)=rhoc(nt-ncb+1:nt,:,1:nt,:nnt-1,:,:)
rhoc(nt+1:,:,1:nt,nnt,:,:)=rhoc(1:ncb,:,1:nt,1,:,:)[image1d(ipx,icy,icz)]
rhoc(nt+1:,:,1:nt,:nnt-1,:,:)=rhoc(1:ncb,:,1:nt,2:,:,:)
sync all
!print*, 'sync x done', sum(rhoc)

!y
rhoc(:,:0,:,:,1,:)=rhoc(:,nt-ncb+1:nt,:,:,nnt,:)[image1d(icx,iny,icz)]
rhoc(:,:0,:,:,2:,:)=rhoc(:,nt-ncb+1:nt,:,:,1:nnt-1,:)
rhoc(:,nt+1:,:,:,nnt,:)=rhoc(:,1:ncb,:,:,1,:)[image1d(icx,ipy,icz)]
rhoc(:,nt+1:,:,:,:nnt-1,:)=rhoc(:,1:ncb,:,1:,2:,:)
sync all
!print*, 'sync y done', sum(rhoc)

!z
rhoc(:,:,:0,:,:,1)=rhoc(:,:,nt-ncb+1:nt,:,:,nnt)[image1d(icx,icy,inz)]
rhoc(:,:,:0,:,:,2:)=rhoc(:,:,nt-ncb+1:nt,:,:,:nnt-1)
rhoc(:,:,nt+1:,:,:,nnt)=rhoc(:,:,1:ncb,:,:,1)[image1d(icx,icy,ipz)]
rhoc(:,:,nt+1:,:,:,:nnt-1)=rhoc(:,:,1:ncb,:,:,2:)
sync all
!print*, 'sync z done', sum(rhoc)

if (sum(rhoc)>npmax) then
  print*, 'npmax should be set larger'
  print*, sum(rhoc), npmax
  stop
endif

!print*, 'x before shift =', sum(x*unit8)
!print*, 'v before shift =', sum(v*unit8)

! move xv to the end, leave enough space for buffer
!print*,sum(rhoc),npmax
nshift=npmax-nplocal
x(:,nshift+1:npmax)=x(:,1:nplocal)
v(:,nshift+1:npmax)=v(:,1:nplocal)
!x(:,:nshift)=0 ! redundant
!v(:,:nshift)=0 ! redundant

#ifdef PID
pid(:,nshift+1:npmax)=pid(:,1:nplocal)
!pid(:,:nshift)=0
#endif

!print*, 'x after shift=', sum(x*unit8)

#ifdef debug
  x=0
  x(1,npmax-nplocal+1:npmax)=1 ! only count particles
#endif

cum=cumsum6(rhoc)
!print*, 'cumsum of all particles =', cum(nt+ncb,nt+ncb,nt+ncb,nnt,nnt,nnt)
! create extended zip

! redistribute local particles -----------------

! nshift is the offset between local (backed up at the end) and extented
ifrom=nshift ! ifrom is the index of shifted particles, contiuous

do itz=1,nnt
do ity=1,nnt
do itx=1,nnt ! loop over tiles
  do iz=1,nt ! loop over z slab
  do iy=1,nt ! loog over y slot
    ! nlast is the last particle's index
    nlast=cum(nt,iy,iz,itx,ity,itz)
    ! nlen is the number of particle in this x-slot
    nlen=nlast-cum(0,iy,iz,itx,ity,itz)
    x(:,nlast-nlen+1:nlast)=x(:,ifrom+1:ifrom+nlen)
    v(:,nlast-nlen+1:nlast)=v(:,ifrom+1:ifrom+nlen)
#ifdef PID
      pid(:,nlast-nlen+1:nlast)=pid(:,ifrom+1:ifrom+nlen)
#endif
    ifrom=ifrom+nlen
  enddo
  enddo
enddo
enddo
enddo
!print*, 'redistributed local particles', sum(x*unit8)
sync all
endsubroutine buffer_density
