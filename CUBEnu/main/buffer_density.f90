subroutine buffer_density
  use variables
  implicit none
  save
  integer(8) nshift,nlen,nlast,ifrom,checkxp0,checkxp1

  if (head) print*, 'buffer_density'
  overhead_image=0
  ! sync buffer regions in rhoc for each tile
  !x
  rhoc(:0,:,:,1,:,:)=rhoc(nt-ncb+1:nt,:,:,nnt,:,:)[image1d(inx,icy,icz)]
  rhoc(:0,:,:,2:,:,:)=rhoc(nt-ncb+1:nt,:,:,:nnt-1,:,:)
  rhoc(nt+1:,:,:,nnt,:,:)=rhoc(1:ncb,:,:,1,:,:)[image1d(ipx,icy,icz)]
  rhoc(nt+1:,:,:,:nnt-1,:,:)=rhoc(1:ncb,:,:,2:,:,:)

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
  !print*, 'sync x done', sum(rhoc)

  !y
  rhoc(:,:0,:,:,1,:)=rhoc(:,nt-ncb+1:nt,:,:,nnt,:)[image1d(icx,iny,icz)]
  rhoc(:,:0,:,:,2:,:)=rhoc(:,nt-ncb+1:nt,:,:,1:nnt-1,:)
  rhoc(:,nt+1:,:,:,nnt,:)=rhoc(:,1:ncb,:,:,1,:)[image1d(icx,ipy,icz)]
  rhoc(:,nt+1:,:,:,:nnt-1,:)=rhoc(:,1:ncb,:,1:,2:,:)
  sync all

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
  !print*, 'sync y done', sum(rhoc)

  !z
  rhoc(:,:,:0,:,:,1)=rhoc(:,:,nt-ncb+1:nt,:,:,nnt)[image1d(icx,icy,inz)]
  rhoc(:,:,:0,:,:,2:)=rhoc(:,:,nt-ncb+1:nt,:,:,:nnt-1)
  rhoc(:,:,nt+1:,:,:,nnt)=rhoc(:,:,1:ncb,:,:,1)[image1d(icx,icy,ipz)]
  rhoc(:,:,nt+1:,:,:,:nnt-1)=rhoc(:,:,1:ncb,:,:,2:)
  sync all

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
  !print*, 'sync z done', sum(rhoc)

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

  

  nshift=np_image_max-nplocal

  checkxp0=sum(xp*int(1,kind=8))

  xp(:,nshift+1:np_image_max)=xp(:,1:nplocal)
  vp(:,nshift+1:np_image_max)=vp(:,1:nplocal)
  xp(:,1:np_image_max-nplocal)=0
  vp(:,1:np_image_max-nplocal)=0

  checkxp1=sum(xp*int(1,kind=8))
  if (checkxp0/=checkxp1) then
    print*, 'wrong1',image,checkxp0,checkxp1
  endif

# ifdef PID
    pid(nshift+1:np_image_max)=pid(1:nplocal)
    pid(1:np_image_max-nplocal)=0
# endif

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
    print*, 'wrong2',image,checkxp0,checkxp1
  endif
  !print*, 'redistributed local particles', sum(x*unit8)
  sync all
endsubroutine buffer_density
