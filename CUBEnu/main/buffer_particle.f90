module buffer_particle_subroutines

contains

subroutine buffer_x
  use variables
  use neutrinos
  implicit none
  save
  call buffer_xvp(.true.,.false.)
#ifdef NEUTRINOS
  if (neutrino_flag) call buffer_xp_nu
#endif
endsubroutine

subroutine buffer_v
  use variables
  use neutrinos
  implicit none
  save
  call buffer_xvp(.false.,.true.)
#ifdef NEUTRINOS
  if (neutrino_flag) call buffer_vp_nu
#endif
endsubroutine

subroutine buffer_xvp(dox,dov)
  use variables
  implicit none
  save
  logical dox,dov
  !integer(8) nshift,nlen,ifrom,mlast
  integer(8),dimension(nt,nt,nnt,nnt) :: myzl,myzr
  integer(8),dimension(1-ncb:nt+ncb,nnt,nnt) :: mzl,mzr
  integer(8),dimension(nnt,nnt) :: ml,mr

  if (head) print*, 'buffer_xp'
# ifdef FORCETEST
    print*,'xp =',xp(:,1:20)
# endif
  call system_clock(t1,t_rate)
  !cum5=spine6(rhoc)

  ! buffer x direction
  ! sync x- buffer with node on the left
  myzl=pprl(:,:,nnt,:,:)[image1d(inx,icy,icz)]
  myzr=ppr0(:,:,nnt,:,:)[image1d(inx,icy,icz)]
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,1 ! do only tile_x=1
    !$omp paralleldo default(shared) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
      xp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(inx,icy,icz)]
      if (dov) &
      vp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
      vp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(inx,icy,icz)]
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x- buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=2,nnt ! skip tile_x=1
    !$omp paralleldo default(shared) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
      xp(:,pprl(iy,iz,itx-1,ity,itz)+1:ppr0(iy,iz,itx-1,ity,itz))
      if (dov) &
      vp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
      vp(:,pprl(iy,iz,itx-1,ity,itz)+1:ppr0(iy,iz,itx-1,ity,itz))
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! sync x+ buffer with node on the right
  myzl=ppl0(:,:,1,:,:)[image1d(ipx,icy,icz)]
  myzr=pplr(:,:,1,:,:)[image1d(ipx,icy,icz)]
  do itz=1,nnt
  do ity=1,nnt
  do itx=nnt,nnt ! do only tile_x=nnt
    !$omp paralleldo default(shared) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
      xp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(ipx,icy,icz)]
      if (dov) &
      vp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
      vp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(ipx,icy,icz)]
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x+ buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt-1 ! skip tile_x=nnt
    !$omp paralleldo default(shared) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
      xp(:,ppl0(iy,iz,itx+1,ity,itz)+1:pplr(iy,iz,itx+1,ity,itz))
      if (dov) &
      vp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
      vp(:,ppl0(iy,iz,itx+1,ity,itz)+1:pplr(iy,iz,itx+1,ity,itz))
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! buffer y direction
  ! sync y-
  mzl=idx_b_l(nt-ncb+1,:,:,nnt,:)[image1d(icx,iny,icz)]
  mzr=idx_b_r(nt,:,:,nnt,:)[image1d(icx,iny,icz)]
  do itz=1,nnt
  do ity=1,1 ! do only ity=1
  do itx=1,nnt
    !$omp paralleldo default(shared) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
      xp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,iny,icz)]
      if (dov) &
      vp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
      vp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,iny,icz)]
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y-
  do itz=1,nnt
  do ity=2,nnt ! skip ity=1
  do itx=1,nnt
    !$omp paralleldo default(shared) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
      xp(:,idx_b_l(nt-ncb+1,iz,itx,ity-1,itz)+1:idx_b_r(nt,iz,itx,ity-1,itz))
      if (dov) &
      vp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
      vp(:,idx_b_l(nt-ncb+1,iz,itx,ity-1,itz)+1:idx_b_r(nt,iz,itx,ity-1,itz))
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! sync y+
  mzl=idx_b_l(1,:,:,1,:)[image1d(icx,ipy,icz)]
  mzr=idx_b_r(ncb,:,:,1,:)[image1d(icx,ipy,icz)]
  do itz=1,nnt
  do ity=nnt,nnt ! do only ity=nnt
  do itx=1,nnt
    !$omp paralleldo default(shared) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
      xp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))
      if (dov) &
      vp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
      vp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y+
  do itz=1,nnt
  do ity=1,nnt-1 ! skip ity=nnt
  do itx=1,nnt
    !$omp paralleldo default(shared) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
      xp(:,idx_b_l(1,iz,itx,ity+1,itz)+1:idx_b_r(ncb,iz,itx,ity+1,itz))
      if (dov) &
      vp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
      vp(:,idx_b_l(1,iz,itx,ity+1,itz)+1:idx_b_r(ncb,iz,itx,ity+1,itz))
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! buffer z direction
  ! sync z-
  ml=idx_b_l(1-ncb,nt-ncb+1,:,:,nnt)[image1d(icx,icy,inz)]
  mr=idx_b_r(nt+ncb,nt,:,:,nnt)[image1d(icx,icy,inz)]
  do itz=1,1 ! do only itz=1
  !$omp paralleldo default(shared) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
    xp(:,ml(itx,ity)+1:mr(itx,ity))
    if (dov) &
    vp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
    vp(:,ml(itx,ity)+1:mr(itx,ity))
  enddo
  enddo
  !$omp endparalleldo
  enddo
  sync all

  ! redistribute z-
  do itz=2,nnt ! skip itz=1
  !$omp paralleldo default(shared) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
    xp(:,idx_b_l(1-ncb,nt-ncb+1,itx,ity,itz-1)+1:idx_b_r(nt+ncb,nt,itx,ity,itz-1))
    if (dov) &
    vp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
    vp(:,idx_b_l(1-ncb,nt-ncb+1,itx,ity,itz-1)+1:idx_b_r(nt+ncb,nt,itx,ity,itz-1))
  enddo
  enddo
  !$omp endparalleldo
  enddo
  sync all

  ! sync z+
  ml=idx_b_l(1-ncb,1,:,:,1)[image1d(icx,icy,ipz)]
  mr=idx_b_r(nt+ncb,ncb,:,:,1)[image1d(icx,icy,ipz)]
  do itz=nnt,nnt ! do only itz=nnt
  !$omp paralleldo default(shared) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
    xp(:,ml(itx,ity)+1:mr(itx,ity))
    if (dov) &
    vp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
    vp(:,ml(itx,ity)+1:mr(itx,ity))
  enddo
  enddo
  !$omp endparalleldo
  enddo
  sync all

  ! redistribute z+
  do itz=1,nnt-1 ! skip itz=nnt
  !$omp paralleldo default(shared) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
    xp(:,idx_b_l(1-ncb,1,itx,ity,itz+1)+1:idx_b_r(nt+ncb,ncb,itx,ity,itz+1))
    if (dov) &
    vp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
    vp(:,idx_b_l(1-ncb,1,itx,ity,itz+1)+1:idx_b_r(nt+ncb,ncb,itx,ity,itz+1))
  enddo
  enddo
  !$omp endparalleldo
  enddo
  sync all
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
# ifdef FORCETEST
    print*,'xp =',xp(:,1:20)
# endif
endsubroutine buffer_xvp

#ifdef NEUTRINOS
subroutine buffer_xp_nu
  use variables
  use neutrinos
  implicit none
  save
  integer(8) nshift,nlen,ifrom,mlast

  if (head) print*, 'buffer_xp_nu'
  ! buffer x direction
  ! sync x- buffer with node on the left
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,1 ! do only tile_x=1
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(0,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,iy,iz,itx,ity,itz)+rhoc_nu(1-ncb,iy,iz,itx,ity,itz)
      mlast=cum_nu(nt,iy,iz,nnt,ity,itz)[image1d(inx,icy,icz)]
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x- buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=2,nnt ! skip tile_x=1
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(0,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,iy,iz,itx,ity,itz)+rhoc_nu(1-ncb,iy,iz,itx,ity,itz)
      mlast=cum_nu(nt,iy,iz,itx-1,ity,itz)
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! sync x+ buffer with node on the right
  do itz=1,nnt
  do ity=1,nnt
  do itx=nnt,nnt ! do only tile_x=nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(nt+ncb,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt,iy,iz,itx,ity,itz)
      mlast=cum_nu(ncb,iy,iz,1,ity,itz)[image1d(ipx,icy,icz)]
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x+ buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt-1 ! skip tile_x=nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(nt+ncb,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt,iy,iz,itx,ity,itz)
      mlast=cum_nu(ncb,iy,iz,itx+1,ity,itz)
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! buffer y direction
  ! sync y-
  do itz=1,nnt
  do ity=1,1 ! do only ity=1
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,0,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,nt,iz,itx,nnt,itz)[image1d(icx,iny,icz)]
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y-
  do itz=1,nnt
  do ity=2,nnt ! skip ity=1
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,0,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,nt,iz,itx,ity-1,itz)
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)
    enddo
  enddo
  enddo
  enddo
  sync all

  ! sync y+
  do itz=1,nnt
  do ity=nnt,nnt ! do only ity=nnt
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,nt+ncb,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt+ncb,nt,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,ncb,iz,itx,1,itz)[image1d(icx,ipy,icz)]
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y+
  do itz=1,nnt
  do ity=1,nnt-1 ! skip ity=nnt
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,nt+ncb,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt+ncb,nt,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,ncb,iz,itx,ity+1,itz)
      xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)
    enddo
  enddo
  enddo
  enddo
  sync all

  ! buffer z direction
  ! sync z-
  do itz=1,1 ! do only itz=1
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,0,itx,ity,itz)
    nlen=nlast-cum_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,nt,itx,ity,nnt)[image1d(icx,icy,inz)]
    xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
  enddo
  enddo
  enddo
  sync all

  ! redistribute z-
  do itz=2,nnt ! skip itz=1
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,0,itx,ity,itz)
    nlen=nlast-cum_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,nt,itx,ity,itz-1)
    xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)
  enddo
  enddo
  enddo
  sync all

  ! sync z+
  do itz=nnt,nnt ! do only itz=nnt
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
    nlen=nlast-cum_nu(nt+ncb,nt+ncb,nt,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,ncb,itx,ity,1)[image1d(icx,icy,ipz)]
    xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
  enddo
  enddo
  enddo
  sync all

  ! redistribute z+
  do itz=1,nnt-1 ! skip itz=nnt
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
    nlen=nlast-cum_nu(nt+ncb,nt+ncb,nt,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,ncb,itx,ity,itz+1)
    xp_nu(:,nlast-nlen+1:nlast)=xp_nu(:,mlast-nlen+1:mlast)
  enddo
  enddo
  enddo
  sync all
endsubroutine buffer_xp_nu
#endif

#ifdef OLD
subroutine buffer_vp
  use variables
  implicit none
  save
  integer(8) nshift,nlen,ifrom,mlast

# ifdef PID
    if (head) print*, 'buffer_vp (vp & pid)'
# else
    if (head) print*, 'buffer_vp'
# endif
  ! buffer x direction
  ! sync x- buffer with node on the left
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,1 ! do only tile_x=1
    do iz=1,nt
    do iy=1,nt
      nlast=cum(0,iy,iz,itx,ity,itz)
      nlen=nlast-cum(1-ncb,iy,iz,itx,ity,itz)+rhoc(1-ncb,iy,iz,itx,ity,itz)
      mlast=cum(nt,iy,iz,nnt,ity,itz)[image1d(inx,icy,icz)]
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x- buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=2,nnt ! skip tile_x=1
    do iz=1,nt
    do iy=1,nt
      nlast=cum(0,iy,iz,itx,ity,itz)
      nlen=nlast-cum(1-ncb,iy,iz,itx,ity,itz)+rhoc(1-ncb,iy,iz,itx,ity,itz)
      mlast=cum(nt,iy,iz,itx-1,ity,itz)
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! sync x+ buffer with node on the right
  do itz=1,nnt
  do ity=1,nnt
  do itx=nnt,nnt ! do only tile_x=nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum(nt+ncb,iy,iz,itx,ity,itz)
      nlen=nlast-cum(nt,iy,iz,itx,ity,itz)
      mlast=cum(ncb,iy,iz,1,ity,itz)[image1d(ipx,icy,icz)]
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x+ buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt-1 ! skip tile_x=nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum(nt+ncb,iy,iz,itx,ity,itz)
      nlen=nlast-cum(nt,iy,iz,itx,ity,itz)
      mlast=cum(ncb,iy,iz,itx+1,ity,itz)
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! buffer y direction
  ! sync y-
  do itz=1,nnt
  do ity=1,1 ! do only ity=1
  do itx=1,nnt
    do iz=1,nt
      nlast=cum(nt+ncb,0,iz,itx,ity,itz)
      nlen=nlast-cum(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc(1-ncb,1-ncb,iz,itx,ity,itz)
      mlast=cum(nt+ncb,nt,iz,itx,nnt,itz)[image1d(icx,iny,icz)]
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y-
  do itz=1,nnt
  do ity=2,nnt ! skip ity=1
  do itx=1,nnt
    do iz=1,nt
      nlast=cum(nt+ncb,0,iz,itx,ity,itz)
      nlen=nlast-cum(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc(1-ncb,1-ncb,iz,itx,ity,itz)
      mlast=cum(nt+ncb,nt,iz,itx,ity-1,itz)
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
    enddo
  enddo
  enddo
  enddo
  sync all

  ! sync y+
  do itz=1,nnt
  do ity=nnt,nnt ! do only ity=nnt
  do itx=1,nnt
    do iz=1,nt
      nlast=cum(nt+ncb,nt+ncb,iz,itx,ity,itz)
      nlen=nlast-cum(nt+ncb,nt,iz,itx,ity,itz)
      mlast=cum(nt+ncb,ncb,iz,itx,1,itz)[image1d(icx,ipy,icz)]
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y+
  do itz=1,nnt
  do ity=1,nnt-1 ! skip ity=nnt
  do itx=1,nnt
    do iz=1,nt
      nlast=cum(nt+ncb,nt+ncb,iz,itx,ity,itz)
      nlen=nlast-cum(nt+ncb,nt,iz,itx,ity,itz)
      mlast=cum(nt+ncb,ncb,iz,itx,ity+1,itz)
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
    enddo
  enddo
  enddo
  enddo
  sync all

  ! buffer z direction
  ! sync z-
  do itz=1,1 ! do only itz=1
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum(nt+ncb,nt+ncb,0,itx,ity,itz)
    nlen=nlast-cum(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc(1-ncb,1-ncb,1-ncb,itx,ity,itz)
    mlast=cum(nt+ncb,nt+ncb,nt,itx,ity,nnt)[image1d(icx,icy,inz)]
# ifdef PID
    pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
# endif
    vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
  enddo
  enddo
  enddo
  sync all

  ! redistribute z-
  do itz=2,nnt ! skip itz=1
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum(nt+ncb,nt+ncb,0,itx,ity,itz)
    nlen=nlast-cum(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc(1-ncb,1-ncb,1-ncb,itx,ity,itz)
    mlast=cum(nt+ncb,nt+ncb,nt,itx,ity,itz-1)
# ifdef PID
    pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
# endif
    vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
  enddo
  enddo
  enddo
  sync all

  ! sync z+
  do itz=nnt,nnt ! do only itz=nnt
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
    nlen=nlast-cum(nt+ncb,nt+ncb,nt,itx,ity,itz)
    mlast=cum(nt+ncb,nt+ncb,ncb,itx,ity,1)[image1d(icx,icy,ipz)]
# ifdef PID
    pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
# endif
    vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
  enddo
  enddo
  enddo
  sync all

  ! redistribute z+
  do itz=1,nnt-1 ! skip itz=nnt
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
    nlen=nlast-cum(nt+ncb,nt+ncb,nt,itx,ity,itz)
    mlast=cum(nt+ncb,nt+ncb,ncb,itx,ity,itz+1)
#   ifdef PID
    pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
    vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
  enddo
  enddo
  enddo
  sync all
endsubroutine buffer_vp
#endif

#ifdef NEUTRINOS
subroutine buffer_vp_nu
  use variables
  use neutrinos
  implicit none
  save
  integer(8) nshift,nlen,ifrom,mlast

# ifdef EID
    if (head) print*, 'buffer_vp_nu (vp_nu & pid_nu)'
# else
    if (head) print*, 'buffer_vp_nu'
# endif
  ! buffer x direction
  ! sync x- buffer with node on the left
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,1 ! do only tile_x=1
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(0,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,iy,iz,itx,ity,itz)+rhoc_nu(1-ncb,iy,iz,itx,ity,itz)
      mlast=cum_nu(nt,iy,iz,nnt,ity,itz)[image1d(inx,icy,icz)]
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x- buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=2,nnt ! skip tile_x=1
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(0,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,iy,iz,itx,ity,itz)+rhoc_nu(1-ncb,iy,iz,itx,ity,itz)
      mlast=cum_nu(nt,iy,iz,itx-1,ity,itz)
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! sync x+ buffer with node on the right
  do itz=1,nnt
  do ity=1,nnt
  do itx=nnt,nnt ! do only tile_x=nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(nt+ncb,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt,iy,iz,itx,ity,itz)
      mlast=cum_nu(ncb,iy,iz,1,ity,itz)[image1d(ipx,icy,icz)]
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x+ buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt-1 ! skip tile_x=nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum_nu(nt+ncb,iy,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt,iy,iz,itx,ity,itz)
      mlast=cum_nu(ncb,iy,iz,itx+1,ity,itz)
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)
    enddo
    enddo
  enddo
  enddo
  enddo
  sync all

  ! buffer y direction
  ! sync y-
  do itz=1,nnt
  do ity=1,1 ! do only ity=1
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,0,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,nt,iz,itx,nnt,itz)[image1d(icx,iny,icz)]
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y-
  do itz=1,nnt
  do ity=2,nnt ! skip ity=1
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,0,iz,itx,ity,itz)
      nlen=nlast-cum_nu(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,nt,iz,itx,ity-1,itz)
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)
    enddo
  enddo
  enddo
  enddo
  sync all

  ! sync y+
  do itz=1,nnt
  do ity=nnt,nnt ! do only ity=nnt
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,nt+ncb,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt+ncb,nt,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,ncb,iz,itx,1,itz)[image1d(icx,ipy,icz)]
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
    enddo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y+
  do itz=1,nnt
  do ity=1,nnt-1 ! skip ity=nnt
  do itx=1,nnt
    do iz=1,nt
      nlast=cum_nu(nt+ncb,nt+ncb,iz,itx,ity,itz)
      nlen=nlast-cum_nu(nt+ncb,nt,iz,itx,ity,itz)
      mlast=cum_nu(nt+ncb,ncb,iz,itx,ity+1,itz)
#   ifdef PID
      pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)
#   endif
      vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)
    enddo
  enddo
  enddo
  enddo
  sync all

  ! buffer z direction
  ! sync z-
  do itz=1,1 ! do only itz=1
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,0,itx,ity,itz)
    nlen=nlast-cum_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,nt,itx,ity,nnt)[image1d(icx,icy,inz)]
# ifdef PID
    pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
# endif
    vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
  enddo
  enddo
  enddo
  sync all

  ! redistribute z-
  do itz=2,nnt ! skip itz=1
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,0,itx,ity,itz)
    nlen=nlast-cum_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc_nu(1-ncb,1-ncb,1-ncb,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,nt,itx,ity,itz-1)
# ifdef PID
    pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)
# endif
    vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)
  enddo
  enddo
  enddo
  sync all

  ! sync z+
  do itz=nnt,nnt ! do only itz=nnt
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
    nlen=nlast-cum_nu(nt+ncb,nt+ncb,nt,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,ncb,itx,ity,1)[image1d(icx,icy,ipz)]
# ifdef PID
    pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
# endif
    vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
  enddo
  enddo
  enddo
  sync all

  ! redistribute z+
  do itz=1,nnt-1 ! skip itz=nnt
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum_nu(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
    nlen=nlast-cum_nu(nt+ncb,nt+ncb,nt,itx,ity,itz)
    mlast=cum_nu(nt+ncb,nt+ncb,ncb,itx,ity,itz+1)
#   ifdef PID
    pid_nu(nlast-nlen+1:nlast)=pid_nu(mlast-nlen+1:mlast)
#   endif
    vp_nu(:,nlast-nlen+1:nlast)=vp_nu(:,mlast-nlen+1:mlast)
  enddo
  enddo
  enddo
  sync all
endsubroutine buffer_vp_nu
#endif

endmodule
