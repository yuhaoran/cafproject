subroutine update_particle
  use variables
  implicit none
  save

  ! x(1:3,:)
  integer(8) ileft,iright,nlast,nlen,idx
  integer(8) g(3),np ! index of grid
  integer(4) rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! double buffer tile
  integer(8) cume(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(4) rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! count writing
  integer(8) ii,jj,kk
  integer(8) iii(3)
  integer(8) pr1,pr2,pr3,pr4

  if (head) print*,'update_particle'
  dt_mid=(dt_old+dt)/2


  overhead_tile=0
  iright=0
  do itz=1,nnt ! loop over tile
  do ity=1,nnt
  do itx=1,nnt
  !if (head) print*,'  tile',int(itx,1),int(ity,1),int(itz,1)
    rhoce=0
    rholocal=0
    x_new=0
    v_new=0
#   ifdef PID
    pid_new=0
#   endif

    !if (head) print*,'    density loop'
    do k=1-ncb,nt+ncb ! loop over coarse grid
    do j=1-ncb,nt+ncb
    do i=1-ncb,nt+ncb
      nlast=cum(i,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        xq=((/i,j,k/)-1d0) + (int(x(:,ip)+ishift,izipx)+rshift)*x_resolution
        vreal=tan((pi*real(v(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/sigma_vi_old)
        deltax=(dt_mid*vreal)/ncell
        g=ceiling(xq+deltax)
        rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1 ! update mesh
        !print*,x(:,ip)
        !print*,xq,deltax,g;stop
      enddo
    enddo
    enddo
    enddo

    !if (head) print*,'    cubesum3'
    cume=cumsum3(rhoce)
    overhead_tile=max(overhead_tile,cume(nt+2*ncb,nt+2*ncb,nt+2*ncb)/real(np_tile_max))

    if (cume(nt+2*ncb,nt+2*ncb,nt+2*ncb)>np_tile_max) then
      print*, '  error: too many particles in this tile+buffer'
      print*, '  ',cume(nt+2*ncb,nt+2*ncb,nt+2*ncb),'>',np_tile_max
      print*, '  on',this_image(), itx,ity,itz
      print*, '  please set tile_buffer larger'
      stop
    endif

    ! create a new x and v for this local tile

  !if (head) print*,'    xvnew loop'
    do k=1-ncb,nt+ncb ! update particle
    do j=1-ncb,nt+ncb
    do i=1-ncb,nt+ncb
      nlast=cum(i,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        xq=((/i,j,k/)-1d0) + (int(x(:,ip)+ishift,izipx)+rshift)*x_resolution
        vreal=tan((pi*real(v(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/sigma_vi_old)
        deltax=(dt_mid*vreal)/ncell
        g=ceiling(xq+deltax)
        rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
        idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3)) ! index for writing

#       ifdef debug
          x_new(:,idx)=x_new(:,idx)+1
#       else
        x_new(:,idx)=x(:,ip)+nint(dt_mid*vreal/(x_resolution*ncell))
  !print*, nlast,np,l,ip
  !print*, x(:,ip)
  !print*, idx
  !print*, x_new(:,idx)
  !stop
#endif
        v_new(:,idx)=v(:,ip)
  !print*, v(:,ip)
  !print*, v_new(:,idx)
  !stop
#ifdef PID
        pid_new(idx)=pid(ip)
#endif
        ! update mesh
        ! rhoce_new
      enddo
    enddo
    enddo
    enddo

#ifdef debug
    print*, sum(abs(x_new(1,1:nptile_old)-1)*1)
#endif
    ! delete particles

    sync all

  !if (head) print*,'    delete_particle loop'
    do k=1,nt
    do j=1,nt
      ileft=iright+1
      nlast=cume(nt,j,k)
      nlen=nlast-cume(0,j,k)
      iright=ileft+nlen-1
      x(:,ileft:iright)=x_new(:,nlast-nlen+1:nlast)
      v(:,ileft:iright)=v_new(:,nlast-nlen+1:nlast)
  !print*,ileft,iright,nlast,nlen
  !print*, x(:,ileft:iright)
  !print*, x_new(:,nlast-nlen+1:nlast)
  !stop

#ifdef PID
      pid(ileft:iright)=pid_new(nlast-nlen+1:nlast)
#endif
    enddo
    enddo

    ! update rhoc
    rhoc(1:nt,1:nt,1:nt,itx,ity,itz)=rhoce(1:nt,1:nt,1:nt)
    !print*, sum(rhoce(1:nt,1:nt,1:nt)), sum(rhoce), cume(nt+2*ncb,nt+2*ncb,nt+2*ncb), sum(rhoc(1:nt,1:nt,1:nt,itx,ity,itz))
  enddo
  enddo
  enddo

  nplocal=iright
  print*,x(:,1),v(:,1)
  do i=1,nn**3
    overhead_tile=max(overhead_tile,overhead_tile[i])
  enddo

  if (head) then
    print*, '  tile overhead',overhead_tile*100,'% full'
    print*, '  comsumed tile_buffer =',overhead_tile*tile_buffer,'/',tile_buffer
    print*, '  clean buffer of rhoc'
  endif

  ! clean up buffer region of rhoc

  rhoc(:0,:,:,:,:,:)=0
  rhoc(nt+1:,:,:,:,:,:)=0
  rhoc(:,:0,:,:,:,:)=0
  rhoc(:,nt+1:,:,:,:,:)=0
  rhoc(:,:,:0,:,:,:)=0
  rhoc(:,:,nt+1:,:,:,:)=0

  ! clean up x, v beyond nplocal
  !x(:,nplocal+1:)=0
  !v(:,nplocal+1:)=0
#ifdef PID
  !pid(:,nplocal+1:)=0
#endif
  !print*,'  check xv of first particle',x(:,1),v(:,1)

  sync all

  if (head) then
    npcheck=0
    do i=1,nn**3
      npcheck=npcheck+nplocal[i]
    enddo
    print*, '  npcheck,npglobal=', npcheck,npglobal
  endif

  sync all

endsubroutine update_particle
