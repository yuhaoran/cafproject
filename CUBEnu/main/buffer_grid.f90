module buffer_grid_subroutines

contains

subroutine buffer_grid
  use variables
  use neutrinos
  implicit none
  save

  call buffer_np(rhoc)
  call buffer_np(rhoc_nu)
  call buffer_vc(vfield)
  call buffer_vc(vfield_nu)
  call redistribute_cdm()
  call redistribute_nu()

endsubroutine

subroutine buffer_np(rhoc)
  use parameters
  implicit none
  save

  integer(4) rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  integer(8),parameter :: unit8=1
  integer(8) itest
  if (head) print*, 'buffer_np'
  !x
  rhoc(:0,:,:,1,:,:)=rhoc(nt-ncb+1:nt,:,:,nnt,:,:)[image1d(inx,icy,icz)]
  rhoc(:0,:,:,2:,:,:)=rhoc(nt-ncb+1:nt,:,:,:nnt-1,:,:)
  rhoc(nt+1:,:,:,nnt,:,:)=rhoc(1:ncb,:,:,1,:,:)[image1d(ipx,icy,icz)]
  rhoc(nt+1:,:,:,:nnt-1,:,:)=rhoc(1:ncb,:,:,2:,:,:)
  sync all
  !y
  rhoc(:,:0,:,:,1,:)=rhoc(:,nt-ncb+1:nt,:,:,nnt,:)[image1d(icx,iny,icz)]
  rhoc(:,:0,:,:,2:,:)=rhoc(:,nt-ncb+1:nt,:,:,1:nnt-1,:)
  rhoc(:,nt+1:,:,:,nnt,:)=rhoc(:,1:ncb,:,:,1,:)[image1d(icx,ipy,icz)]
  rhoc(:,nt+1:,:,:,:nnt-1,:)=rhoc(:,1:ncb,:,1:,2:,:)
  sync all
  !z
  rhoc(:,:,:0,:,:,1)=rhoc(:,:,nt-ncb+1:nt,:,:,nnt)[image1d(icx,icy,inz)]
  rhoc(:,:,:0,:,:,2:)=rhoc(:,:,nt-ncb+1:nt,:,:,:nnt-1)
  rhoc(:,:,nt+1:,:,:,nnt)=rhoc(:,:,1:ncb,:,:,1)[image1d(icx,icy,ipz)]
  rhoc(:,:,nt+1:,:,:,:nnt-1)=rhoc(:,:,1:ncb,:,:,2:)
  sync all
endsubroutine

subroutine buffer_vc(vfield)
  use parameters
  implicit none
  save

  real(4) vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
  ! the following variables are introduced because
  ! gcc only allows <= 7 ranks in arrays
  real(4) vtransx(3,ncb,nt+2*ncb,nt+2*ncb,nnt,nnt)[*]
  real(4) vtransy(3,nt+2*ncb,ncb,nt+2*ncb,nnt,nnt)[*]
  real(4) vtransz(3,nt+2*ncb,nt+2*ncb,ncb,nnt,nnt)[*]
  if (head) print*, 'buffer_vc'
  !x
  vtransx=vfield(:,nt-ncb+1:nt,:,:,nnt,:,:)
  sync all
  vfield(:,:0,:,:,1,:,:)=vtransx(:,:,:,:,:,:)[image1d(inx,icy,icz)]
  sync all
  vfield(:,:0,:,:,2:,:,:)=vfield(:,nt-ncb+1:nt,:,:,:nnt-1,:,:)

  vtransx=vfield(:,1:ncb,:,:,1,:,:)
  sync all
  vfield(:,nt+1:,:,:,nnt,:,:)=vtransx(:,:,:,:,:,:)[image1d(ipx,icy,icz)]
  sync all
  vfield(:,nt+1:,:,:,:nnt-1,:,:)=vfield(:,1:ncb,:,:,2:,:,:)
  sync all

  !y
  vtransy=vfield(:,:,nt-ncb+1:nt,:,:,nnt,:)
  sync all
  vfield(:,:,:0,:,:,1,:)=vtransy(:,:,:,:,:,:)[image1d(icx,iny,icz)]
  sync all
  vfield(:,:,:0,:,:,2:,:)=vfield(:,:,nt-ncb+1:nt,:,:,1:nnt-1,:)

  vtransy=vfield(:,:,1:ncb,:,:,1,:)
  sync all
  vfield(:,:,nt+1:,:,:,nnt,:)=vtransy(:,:,:,:,:,:)[image1d(icx,ipy,icz)]
  sync all
  vfield(:,:,nt+1:,:,:,:nnt-1,:)=vfield(:,:,1:ncb,:,1:,2:,:)
  sync all

  !z
  vtransz=vfield(:,:,:,nt-ncb+1:nt,:,:,nnt)
  sync all
  vfield(:,:,:,:0,:,:,1)=vtransz(:,:,:,:,:,:)[image1d(icx,icy,inz)]
  sync all
  vfield(:,:,:,:0,:,:,2:)=vfield(:,:,:,nt-ncb+1:nt,:,:,:nnt-1)

  vtransz=vfield(:,:,:,1:ncb,:,:,1)
  sync all
  vfield(:,:,:,nt+1:,:,:,nnt)=vtransz(:,:,:,:,:,:)[image1d(icx,icy,ipz)]
  sync all
  vfield(:,:,:,nt+1:,:,:,:nnt-1)=vfield(:,:,:,1:ncb,:,:,2:)
  sync all
endsubroutine

subroutine redistribute_cdm()
  use variables
  implicit none
  save
  integer(8) nshift,nlen,nlast,ifrom,checkxp0,checkxp1

  if (head) print*, 'redistribute_cdm'
  ! check
  overhead_image=sum(rhoc*unit8)/real(np_image_max,8)
  sync all
  do i=1,nn**3
    overhead_image=max(overhead_image,overhead_image[i])
  enddo
  if (head) then
    print*, '  image overhead',overhead_image*100,'% full'
    print*, '  comsumed image_buffer =',overhead_image*image_buffer,'/',image_buffer
  endif
  sync all

  if (overhead_image>1d0) then
    print*, '  error: too many particles in this image+buffer'
    print*, '  ',sum(rhoc*unit8),np_image_max
    print*, '  on',this_image()
    print*, '  please set image_buffer larger'
    stop
  endif

  ! shift to right
  checkxp0=sum(xp*int(1,kind=8))
  nshift=np_image_max-nplocal
  xp(:,nshift+1:np_image_max)=xp(:,1:nplocal)
  xp(:,1:nshift)=0
  vp(:,nshift+1:np_image_max)=vp(:,1:nplocal)
  vp(:,1:nshift)=0
# ifdef PID
    pid(nshift+1:np_image_max)=pid(1:nplocal)
    pid(1:nshift)=0
# endif
  checkxp1=sum(xp*int(1,kind=8))
  print*, '  ',np_image_max,nplocal,nshift
  if (checkxp0/=checkxp1) then
    print*, '  error in shifting right',image,checkxp0,checkxp1
    stop
  endif

  ! shift back
  cum=cumsum6(rhoc)
  ifrom=nshift

  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt ! loop over tiles
    do iz=1,nt ! loop over z slab
    do iy=1,nt ! loog over y slot
      ! nlast is the last particle's index
      nlast=cum(nt,iy,iz,itx,ity,itz)
      ! nlen is the number of particle in this x-slot
      nlen=nlast-cum(0,iy,iz,itx,ity,itz)
      xp(:,nlast-nlen+1:nlast)=xp(:,ifrom+1:ifrom+nlen)
      vp(:,nlast-nlen+1:nlast)=vp(:,ifrom+1:ifrom+nlen)
      xp(:,ifrom+1:ifrom+nlen)=0
      vp(:,ifrom+1:ifrom+nlen)=0
#     ifdef PID
        pid(nlast-nlen+1:nlast)=pid(ifrom+1:ifrom+nlen)
        pid(ifrom+1:ifrom+nlen)=0
#     endif
      ifrom=ifrom+nlen
    enddo
    enddo
  enddo
  enddo
  enddo

  checkxp1=sum(xp*int(1,kind=8))
  if (checkxp0/=checkxp1) then
    print*, '  error in shifting back',image,checkxp0,checkxp1
    stop
  endif
  sync all
endsubroutine

subroutine redistribute_nu
  use variables
  use neutrinos
  implicit none
  save
  integer(8) nshift,nlen,nlast,ifrom,checkxp0,checkxp1

  if (head) print*, 'redistribute_nu'
  ! check
  overhead_image=sum(rhoc_nu*unit8)/real(np_image_max_nu,8)
  sync all
  do i=1,nn**3
    overhead_image=max(overhead_image,overhead_image[i])
  enddo
  if (head) then
    print*, '  image overhead',overhead_image*100,'% full'
    print*, '  comsumed image_buffer_nu =',overhead_image*image_buffer_nu,'/',image_buffer_nu
  endif
  sync all

  if (overhead_image>1d0) then
    print*, '  error: too many particles in this image+buffer'
    print*, '  ',sum(rhoc_nu*unit8),np_image_max_nu
    print*, '  on',this_image()
    print*, '  please set image_buffer_nu larger'
    stop
  endif

  ! shift to right
  checkxp0=sum(xp*int(1,kind=8))
  nshift=np_image_max_nu-nplocal_nu
  xp_nu(:,nshift+1:np_image_max_nu)=xp_nu(:,1:nplocal_nu)
  xp_nu(:,1:nshift)=0
  vp_nu(:,nshift+1:np_image_max_nu)=vp_nu(:,1:nplocal_nu)
  vp_nu(:,1:nshift)=0
# ifdef EID
    pid_nu(nshift+1:np_image_max_nu)=pid(1:nplocal_nu)
    pid_nu(1:nshift)=0
# endif
  checkxp1=sum(xp*int(1,kind=8))
  print*, '  ',np_image_max_nu,nplocal_nu,nshift
  if (checkxp0/=checkxp1) then
    print*, '  error in shifting right',image,checkxp0,checkxp1
    stop
  endif

  ! shift back
  cum_nu=cumsum6(rhoc_nu)
  ifrom=nshift

  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt ! loop over tiles
    do iz=1,nt ! loop over z slab
    do iy=1,nt ! loog over y slot
      ! nlast is the last particle's index
      nlast=cum_nu(nt,iy,iz,itx,ity,itz)
      ! nlen is the number of particle in this x-slot
      nlen=nlast-cum_nu(0,iy,iz,itx,ity,itz)
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,ifrom+1:ifrom+nlen)
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,ifrom+1:ifrom+nlen)
      xp_nu(:,ifrom+1:ifrom+nlen)=0
      vp_nu(:,ifrom+1:ifrom+nlen)=0
#     ifdef EID
        pid_nu(nlast-nlen+1:nlast)=pid(ifrom+1:ifrom+nlen)
        pid_nu(ifrom+1:ifrom+nlen)=0
#     endif
      ifrom=ifrom+nlen
    enddo
    enddo
  enddo
  enddo
  enddo

  checkxp1=sum(xp*int(1,kind=8))
  if (checkxp0/=checkxp1) then
    print*, '  error in shifting back',image,checkxp0,checkxp1
    stop
  endif
  sync all

endsubroutine

endmodule
