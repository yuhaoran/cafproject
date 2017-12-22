subroutine buffer_xp
  use variables
  implicit none
  save
  integer(8) nshift,nlen,nlast,ifrom,mlast

  if (head) print*, 'buffer_xp'
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
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo

  ! redistribute local x- buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=2,nnt ! skip tile_x=1
    do iz=1,nt
    do iy=1,nt
      nlast=cum(0,iy,iz,itx,ity,itz)
      nlen=nlast-cum(1-ncb,iy,iz,itx,ity,itz)+rhoc(1-ncb,iy,iz,itx,ity,itz)
      mlast=cum(nt,iy,iz,itx-1,ity,itz)
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
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
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
    enddo
    enddo
  enddo
  enddo
  enddo

  ! redistribute local x+ buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt-1 ! skip tile_x=nnt
    do iz=1,nt
    do iy=1,nt
      nlast=cum(nt+ncb,iy,iz,itx,ity,itz)
      nlen=nlast-cum(nt,iy,iz,itx,ity,itz)
      mlast=cum(ncb,iy,iz,itx+1,ity,itz)
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
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
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
    enddo
  enddo
  enddo
  enddo

  ! redistribute y-
  do itz=1,nnt
  do ity=2,nnt ! skip ity=1
  do itx=1,nnt
    do iz=1,nt
      nlast=cum(nt+ncb,0,iz,itx,ity,itz)
      nlen=nlast-cum(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc(1-ncb,1-ncb,iz,itx,ity,itz)
      mlast=cum(nt+ncb,nt,iz,itx,ity-1,itz)
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
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
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
    enddo
  enddo
  enddo
  enddo

  ! redistribute y+
  do itz=1,nnt
  do ity=1,nnt-1 ! skip ity=nnt
  do itx=1,nnt
    do iz=1,nt
      nlast=cum(nt+ncb,nt+ncb,iz,itx,ity,itz)
      nlen=nlast-cum(nt+ncb,nt,iz,itx,ity,itz)
      mlast=cum(nt+ncb,ncb,iz,itx,ity+1,itz)
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
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
    xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
  enddo
  enddo
  enddo

  ! redistribute z-
  do itz=2,nnt ! skip itz=1
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum(nt+ncb,nt+ncb,0,itx,ity,itz)
    nlen=nlast-cum(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc(1-ncb,1-ncb,1-ncb,itx,ity,itz)
    mlast=cum(nt+ncb,nt+ncb,nt,itx,ity,itz-1)
    xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
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
    xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
  enddo
  enddo
  enddo

  ! redistribute z+
  do itz=1,nnt-1 ! skip itz=nnt
  do ity=1,nnt
  do itx=1,nnt
    nlast=cum(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
    nlen=nlast-cum(nt+ncb,nt+ncb,nt,itx,ity,itz)
    mlast=cum(nt+ncb,nt+ncb,ncb,itx,ity,itz+1)
    xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
  enddo
  enddo
  enddo
  sync all

endsubroutine buffer_xp

subroutine buffer_xp_nu
  use variables
  use neutrinos
  implicit none
  save
  integer(8) nshift,nlen,nlast,ifrom,mlast

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
