module buffer_grid_subroutines

contains

subroutine buffer_grid
  use variables
  use neutrinos
  implicit none
  save

  call buffer_np(rhoc)
# ifdef NEUTRINOS
    if (neutrino_flag) call buffer_np(rhoc_nu)
# endif

  call buffer_vc(vfield)
# ifdef NEUTRINOS
    if (neutrino_flag) call buffer_vc(vfield_nu)
# endif

  call redistribute_cdm()
# ifdef NEUTRINOS
    if (neutrino_flag) call redistribute_nu()
# endif

endsubroutine

subroutine buffer_np(rhoc)
  use omp_lib
  use parameters
  implicit none
  save

  integer(4) rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  integer(8),parameter :: unit8=1
  integer(8) itest
  if (head) then
    print*, 'buffer_np'
  endif
# ifdef FORCETEST
    print*,'sum of rhoc =',sum(rhoc)
# endif
  !x
  !$omp workshare
  rhoc(:0,:,:,1,:,:)=rhoc(nt-ncb+1:nt,:,:,nnt,:,:)[image1d(inx,icy,icz)]
  rhoc(:0,:,:,2:,:,:)=rhoc(nt-ncb+1:nt,:,:,:nnt-1,:,:)
  rhoc(nt+1:,:,:,nnt,:,:)=rhoc(1:ncb,:,:,1,:,:)[image1d(ipx,icy,icz)]
  rhoc(nt+1:,:,:,:nnt-1,:,:)=rhoc(1:ncb,:,:,2:,:,:)
  !$omp endworkshare
  sync all
  !y
  !$omp workshare
  rhoc(:,:0,:,:,1,:)=rhoc(:,nt-ncb+1:nt,:,:,nnt,:)[image1d(icx,iny,icz)]
  rhoc(:,:0,:,:,2:,:)=rhoc(:,nt-ncb+1:nt,:,:,1:nnt-1,:)
  rhoc(:,nt+1:,:,:,nnt,:)=rhoc(:,1:ncb,:,:,1,:)[image1d(icx,ipy,icz)]
  rhoc(:,nt+1:,:,:,:nnt-1,:)=rhoc(:,1:ncb,:,1:,2:,:)
  !$omp endworkshare
  sync all
  !z
  !$omp workshare
  rhoc(:,:,:0,:,:,1)=rhoc(:,:,nt-ncb+1:nt,:,:,nnt)[image1d(icx,icy,inz)]
  rhoc(:,:,:0,:,:,2:)=rhoc(:,:,nt-ncb+1:nt,:,:,:nnt-1)
  rhoc(:,:,nt+1:,:,:,nnt)=rhoc(:,:,1:ncb,:,:,1)[image1d(icx,icy,ipz)]
  rhoc(:,:,nt+1:,:,:,:nnt-1)=rhoc(:,:,1:ncb,:,:,2:)
  !$omp endworkshare
# ifdef FORCETEST
    print*,'sum of rhoc =',sum(rhoc)
# endif
  sync all
endsubroutine

subroutine buffer_vc(vfield1)
  use variables
  implicit none
  save

  real(4) vfield1(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
  ! the following variables are introduced because
  ! gcc only allows <= 7 ranks in arrays
  !real(4) vtransx(3,ncb,nt+2*ncb,nt+2*ncb,nnt,nnt)[*]
  !real(4) vtransy(3,nt+2*ncb,ncb,nt+2*ncb,nnt,nnt)[*]
  !real(4) vtransz(3,nt+2*ncb,nt+2*ncb,ncb,nnt,nnt)[*]
  if (head) print*, 'buffer_vc'
  !x
  vtransx=vfield1(:,nt-ncb+1:nt,:,:,nnt,:,:)
  sync all
  vfield1(:,:0,:,:,1,:,:)=vtransx(:,:,:,:,:,:)[image1d(inx,icy,icz)]
  sync all
  vfield1(:,:0,:,:,2:,:,:)=vfield1(:,nt-ncb+1:nt,:,:,:nnt-1,:,:)

  vtransx=vfield1(:,1:ncb,:,:,1,:,:)
  sync all
  vfield1(:,nt+1:,:,:,nnt,:,:)=vtransx(:,:,:,:,:,:)[image1d(ipx,icy,icz)]
  sync all
  vfield1(:,nt+1:,:,:,:nnt-1,:,:)=vfield1(:,1:ncb,:,:,2:,:,:)
  sync all

  !y
  vtransy=vfield1(:,:,nt-ncb+1:nt,:,:,nnt,:)
  sync all
  vfield1(:,:,:0,:,:,1,:)=vtransy(:,:,:,:,:,:)[image1d(icx,iny,icz)]
  sync all
  vfield1(:,:,:0,:,:,2:,:)=vfield1(:,:,nt-ncb+1:nt,:,:,1:nnt-1,:)

  vtransy=vfield1(:,:,1:ncb,:,:,1,:)
  sync all
  vfield1(:,:,nt+1:,:,:,nnt,:)=vtransy(:,:,:,:,:,:)[image1d(icx,ipy,icz)]
  sync all
  vfield1(:,:,nt+1:,:,:,:nnt-1,:)=vfield1(:,:,1:ncb,:,1:,2:,:)
  sync all

  !z
  vtransz=vfield1(:,:,:,nt-ncb+1:nt,:,:,nnt)
  sync all
  vfield1(:,:,:,:0,:,:,1)=vtransz(:,:,:,:,:,:)[image1d(icx,icy,inz)]
  sync all
  vfield1(:,:,:,:0,:,:,2:)=vfield1(:,:,:,nt-ncb+1:nt,:,:,:nnt-1)

  vtransz=vfield1(:,:,:,1:ncb,:,:,1)
  sync all
  vfield1(:,:,:,nt+1:,:,:,nnt)=vtransz(:,:,:,:,:,:)[image1d(icx,icy,ipz)]
  sync all
  vfield1(:,:,:,nt+1:,:,:,:nnt-1)=vfield1(:,:,:,1:ncb,:,:,2:)
  sync all
endsubroutine

subroutine redistribute_cdm()
  use variables
  use neutrinos
  implicit none
  save
  integer(8) nshift,nlen,ifrom,checkxp0,checkxp1

  if (head) then
    print*,''
    print*, 'redistribute_cdm'
    call system_clock(t1,t_rate)
  endif
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
  nshift=np_image_max-sim%nplocal
  !$omp parallelsections default(shared)
  !$omp section
    xp(:,nshift+1:np_image_max)=xp(:,1:sim%nplocal)
    xp(:,1:nshift)=0
  !$omp section
    vp(:,nshift+1:np_image_max)=vp(:,1:sim%nplocal)
    vp(:,1:nshift)=0
# ifdef PID
    !$omp section
    pid(nshift+1:np_image_max)=pid(1:sim%nplocal)
    pid(1:nshift)=0
# endif
  !$omp endparallelsections
  checkxp1=sum(xp*int(1,kind=8))

  if (checkxp0/=checkxp1) then
    print*, '  error in shifting right',image,checkxp0,checkxp1
    stop
  endif
  ! shift back
  call spine_image(rhoc,idx_b_l,idx_b_r,ppl0,pplr,pprl,ppr0,ppl,ppr)

  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do iz=1,nt
    do iy=1,nt
      xp(:,ppl0(iy,iz,itx,ity,itz)+1:ppr0(iy,iz,itx,ity,itz))=xp(:,nshift+ppl(iy,iz,itx,ity,itz)+1:nshift+ppr(iy,iz,itx,ity,itz))
      vp(:,ppl0(iy,iz,itx,ity,itz)+1:ppr0(iy,iz,itx,ity,itz))=vp(:,nshift+ppl(iy,iz,itx,ity,itz)+1:nshift+ppr(iy,iz,itx,ity,itz))
#       ifdef PID
        pid(ppl0(iy,iz,itx,ity,itz)+1:ppr0(iy,iz,itx,ity,itz))=pid(nshift+ppl(iy,iz,itx,ity,itz)+1:nshift+ppr(iy,iz,itx,ity,itz))
#       endif
    enddo
    enddo
  enddo
  enddo
  enddo
  !xp(:,idx_b_r(nt+ncb,nt+ncb,nnt,nnt,nnt):)=0

  !checkxp1=sum(xp*int(1,kind=8))
  !if (checkxp0/=checkxp1) then
    !print*, '  error in shifting back',image,checkxp0,checkxp1
    !stop
  !endif
  if (head) then
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, ''
  endif
  sync all
endsubroutine

#ifdef NEUTRINOS
subroutine redistribute_nu
  use variables
  use neutrinos
  implicit none
  save
  integer(8) nshift,nlen,ifrom,checkxp0,checkxp1

  if (head) then
    print*,''
    print*, 'redistribute_nu'
    call system_clock(t1,t_rate)
  endif
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
  nshift=np_image_max_nu-sim%nplocal_nu
  !$omp parallelsections
  !$omp section
    xp_nu(:,nshift+1:np_image_max_nu)=xp_nu(:,1:sim%nplocal_nu)
    xp_nu(:,1:nshift)=0
  !$omp section
    vp_nu(:,nshift+1:np_image_max_nu)=vp_nu(:,1:sim%nplocal_nu)
    vp_nu(:,1:nshift)=0
# ifdef EID
  !$omp section
    pid_nu(nshift+1:np_image_max_nu)=pid_nu(1:sim%nplocal_nu)
    pid_nu(1:nshift)=0
# endif
  !$omp endparallelsections
  checkxp1=sum(xp*int(1,kind=8))
  !print*, '  ',np_image_max_nu,sim%nplocal_nu,nshift
  if (checkxp0/=checkxp1) then
    print*, '  error in shifting right',image,checkxp0,checkxp1
    stop
  endif

  ! shift back
  cum_nu=cumsum6(rhoc_nu)
  ifrom=nshift
  !$omp parallelsections default(shared) &
  !$omp& private(itz,ity,itx,iz,iy,nlast,nlen,ifrom)
  !$omp section
  ifrom=nshift
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(nt,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(0,iy,iz,itx,ity,itz)
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,ifrom+1:ifrom+nlen)
      xp_nu(:,ifrom+1:ifrom+nlen)=0
      ifrom=ifrom+nlen
    enddo
    enddo
  enddo
  enddo
  enddo
  !$omp section
  ifrom=nshift
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(nt,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(0,iy,iz,itx,ity,itz)
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,ifrom+1:ifrom+nlen)
      vp_nu(:,ifrom+1:ifrom+nlen)=0
      ifrom=ifrom+nlen
    enddo
    enddo
  enddo
  enddo
  enddo
#ifdef EID
  !$omp section
  ifrom=nshift
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(nt,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(0,iy,iz,itx,ity,itz)
      pid_nu(nlast-nlen+1:nlast)=pid_nu(ifrom+1:ifrom+nlen)
      pid_nu(ifrom+1:ifrom+nlen)=0
      ifrom=ifrom+nlen
    enddo
    enddo
  enddo
  enddo
  enddo
#endif
  !$omp endparallelsections

  checkxp1=sum(xp*int(1,kind=8))
  if (checkxp0/=checkxp1) then
    print*, '  error in shifting back',image,checkxp0,checkxp1
    stop
  endif
  if (head) then
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, ''
  endif
  sync all

endsubroutine
#endif

endmodule
