module update_particle

contains

subroutine update_x
  use variables
  use neutrinos
  implicit none
  save
  call update_xp
  call update_xp_nu
endsubroutine

subroutine update_xp()
  use variables
  implicit none
  save

  integer(8) ileft,iright,nlen,idx
  integer(8) g(3) ! index of grid
  integer(4) rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! double buffer tile
  real(4) vfield_new(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  real(8),parameter :: weight_v=0.1 ! how previous-step vfield is mostly weighted
  integer(8) cume(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(4) rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! count writing
  integer(8) checkv0,checkv1

  if (head) print*,'update_xp'
  call system_clock(t1,t_rate)
  dt_mid=(dt_old+dt)/2
  overhead_tile=0
  iright=0
  do itz=1,nnt ! loop over tile
  do ity=1,nnt
  do itx=1,nnt
  !if (head) print*,'  tile',int(itx,1),int(ity,1),int(itz,1)
    rhoce=0
    rholocal=0
    ! for empty coarse cells, use previous-step vfield
    vfield_new=0
    vfield_new(:,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb)=vfield(:,:,:,:,itx,ity,itz)*weight_v
    xp_new=0
    vp_new=0
#   ifdef PID
    pid_new=0
#   endif
    !if (head) print*,'    density loop'
    !$omp paralleldo &
    !$omp& default(shared) &
    !$omp& private(k,j,i,nlast,np,l,ip,xq,vreal,deltax,g) 
!    !$omp& reduction(+:rhoce,vfield_new)
    do k=1-ncb,nt+ncb ! loop over coarse grid
    do j=1-ncb,nt+ncb
    do i=1-ncb,nt+ncb
      nlast=cum(i,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
        vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
        vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
        deltax=(dt_mid*vreal)/ncell
        g=ceiling(xq+deltax)
        rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1 ! update mesh
        vfield_new(:,g(1),g(2),g(3))=vfield_new(:,g(1),g(2),g(3))+vreal
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    ! vfield_new is kept the same as previous-step if the grid is empty.
    ! vfield_new is at most weight_v weighted by previous-step for non-empty grids.
    ! this also gets rid of if statements.
    vfield_new(1,:,:,:)=vfield_new(1,:,:,:)/(rhoce+weight_v)
    vfield_new(2,:,:,:)=vfield_new(2,:,:,:)/(rhoce+weight_v)
    vfield_new(3,:,:,:)=vfield_new(3,:,:,:)/(rhoce+weight_v)
    cume=cumsum3(rhoce)
    overhead_tile=max(overhead_tile,cume(nt+2*ncb,nt+2*ncb,nt+2*ncb)/real(np_tile_max))
    if (cume(nt+2*ncb,nt+2*ncb,nt+2*ncb)>np_tile_max) then
      print*, '  error: too many particles in this tile+buffer'
      print*, '  ',cume(nt+2*ncb,nt+2*ncb,nt+2*ncb),'>',np_tile_max
      print*, '  on',this_image(), itx,ity,itz
      print*, '  please set tile_buffer larger'
      stop
    endif
    ! create x_new v_new for this local tile
    !if (head) print*,'    xvnew loop'
    !$omp paralleldo &
    !$omp& default(shared) &
    !$omp& private(k,j,i,nlast,np,l,ip,xq,vreal,deltax,g,idx)
!    !$omp& reduction(+:rholocal)
    do k=1-ncb,nt+ncb ! update particle
    do j=1-ncb,nt+ncb
    do i=1-ncb,nt+ncb
      nlast=cum(i,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
        vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
        vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
        deltax=(dt_mid*vreal)/ncell
        g=ceiling(xq+deltax)
        rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
        idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3))
        xp_new(:,idx)=xp(:,ip)+nint(dt_mid*vreal/(x_resolution*ncell))
        vreal=vreal-vfield_new(:,g(1),g(2),g(3))
        vp_new(:,idx)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
#       ifdef PID
          pid_new(idx)=pid(ip)
#       endif
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    sync all
    ! delete particles
    !if (head) print*,'    delete_particle loop'
    do k=1,nt
    do j=1,nt
      ileft=iright+1
      nlast=cume(nt,j,k)
      nlen=nlast-cume(0,j,k)
      iright=ileft+nlen-1
      xp(:,ileft:iright)=xp_new(:,nlast-nlen+1:nlast)
      vp(:,ileft:iright)=vp_new(:,nlast-nlen+1:nlast)
#     ifdef PID
        pid(ileft:iright)=pid_new(nlast-nlen+1:nlast)
#     endif
    enddo
    enddo
    ! update rhoc
    rhoc(:,:,:,itx,ity,itz)=0
    rhoc(1:nt,1:nt,1:nt,itx,ity,itz)=rhoce(1:nt,1:nt,1:nt)
    vfield(:,:,:,:,itx,ity,itz)=0
    vfield(:,1:nt,1:nt,1:nt,itx,ity,itz)=vfield_new(:,1:nt,1:nt,1:nt)
  enddo
  enddo
  enddo ! end looping over tiles
  nplocal=iright
  xp(:,nplocal+1:)=0
  vp(:,nplocal+1:)=0
# ifdef PID
    pid(nplocal+1:)=0
# endif
  !nplocal=sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
  !print*, iright,sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
  !stop
  sync all

  ! calculate std of the velocity field
  cum=cumsum6(rhoc)
  std_vsim=0; std_vsim_c=0; std_vsim_res=0
  !$omp paralleldo &
  !$omp& default(shared) &
  !$omp& private(itz,ity,itx,k,j,i,nlast,np,l,ip,vreal) &
  !$omp& reduction(+:std_vsim_c,std_vsim,std_vsim_res)
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do k=1,nt
    do j=1,nt
    do i=1,nt
      std_vsim_c=std_vsim_c+sum(vfield(:,i,j,k,itx,ity,itz)**2)
      nlast=cum(i,j,k,itx,ity,itz)
      np=rhoc(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        std_vsim_res=std_vsim_res+sum(vreal**2)
        vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
        std_vsim=std_vsim+sum(vreal**2)
      enddo
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  !$omp endparalleldo
  sync all

  ! co_sum
  if (head) then
    do i=2,nn**3
      std_vsim_c=std_vsim_c+std_vsim_c[i]
      std_vsim_res=std_vsim_res+std_vsim_res[i]
      std_vsim=std_vsim+std_vsim[i]
    enddo
  endif
  sync all

  ! broadcast
  std_vsim_c=std_vsim_c[1]
  std_vsim_res=std_vsim_res[1]
  std_vsim=std_vsim[1]
  sync all

  ! divide
  std_vsim=sqrt(std_vsim/npglobal)
  std_vsim_c=sqrt(std_vsim_c/nc/nc/nc/nn/nn/nn)
  std_vsim_res=sqrt(std_vsim_res/npglobal)

  ! set sigma_vi_new according to particle statistics
  sigma_vi_new=std_vsim_res/sqrt(3.)
  if (head) then
    print*,'  std_vsim    ',std_vsim*sim%vsim2phys,'km/s'
    print*,'  std_vsim_c  ',std_vsim_c*sim%vsim2phys,'km/s'
    print*,'  std_vsim_res',std_vsim_res*sim%vsim2phys,'km/s'
    print*,'  sigma_vi    ',sigma_vi,'(simulation unit)'
    print*,'  sigma_vi_new',sigma_vi_new,'(simulation unit)'
    write(77) a-da,real((/std_vsim,std_vsim_c,std_vsim_res/))
  endif
  sync all

  do i=1,nn**3
    overhead_tile=max(overhead_tile,overhead_tile[i])
  enddo

  if (head) then
    print*, '  tile overhead',overhead_tile*100,'% full'
    print*, '  comsumed tile_buffer =',overhead_tile*tile_buffer,'/',tile_buffer
    print*, '  clean buffer of rhoc'
  endif

  if (head) then
    npcheck=0
    do i=1,nn**3
      npcheck=npcheck+nplocal[i]
    enddo
    print*, '  npcheck,npglobal=', npcheck,npglobal
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, ''
  endif
  sync all
endsubroutine update_xp







subroutine update_xp_nu()
  use variables
  use neutrinos
  implicit none
  save

  integer(8) ileft,iright,nlen,idx
  integer(8) g(3) ! index of grid
  integer(4) rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! double buffer tile
  real(4) vfield_new(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  real(8),parameter :: weight_v=0.1 ! how previous-step vfield is mostly weighted
  integer(8) cume(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(4) rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! count writing
  integer(8) checkv0,checkv1

  if (head) print*,'update_vp_nu'
  call system_clock(t1,t_rate)
  dt_mid=(dt_old+dt)/2
  overhead_tile=0
  iright=0
  do itz=1,nnt ! loop over tile
  do ity=1,nnt
  do itx=1,nnt
    rhoce=0
    rholocal=0
    vfield_new=0
    vfield_new(:,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb)=vfield_nu(:,:,:,:,itx,ity,itz)*weight_v
    xp_new_nu=0
    vp_new_nu=0
#   ifdef EID
    pid_new_nu=0
#   endif
    !if (head) print*,'    density loop'
    !$omp paralleldo &
    !$omp& default(shared) &
    !$omp& private(k,j,i,nlast,np,l,ip,xq,vreal,deltax,g)
!    !$omp& reduction(+:rhoce,vfield_new)
    do k=1-ncb,nt+ncb ! loop over coarse grid
    do j=1-ncb,nt+ncb
    do i=1-ncb,nt+ncb
      nlast=cum_nu(i,j,k,itx,ity,itz)
      np=rhoc_nu(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        xq=((/i,j,k/)-1d0) + (int(xp_nu(:,ip)+ishift_nu,izipx_nu)+rshift_nu)*x_resolution_nu
        vreal=tan((pi*real(vp_nu(:,ip)))/real(nvbin_nu-1)) / (sqrt(pi/2)/(sigma_vi_nu*vrel_boost))
        vreal=vreal+vfield_nu(:,i,j,k,itx,ity,itz)
        deltax=(dt_mid*vreal)/ncell
        g=ceiling(xq+deltax)
        rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1 ! update mesh
        vfield_new(:,g(1),g(2),g(3))=vfield_new(:,g(1),g(2),g(3))+vreal
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    ! vfield_new is kept the same as previous-step if the grid is empty.
    ! vfield_new is at most weight_v weighted by previous-step for non-empty grids.
    ! this also gets rid of if statements.
    vfield_new(1,:,:,:)=vfield_new(1,:,:,:)/(rhoce+weight_v)
    vfield_new(2,:,:,:)=vfield_new(2,:,:,:)/(rhoce+weight_v)
    vfield_new(3,:,:,:)=vfield_new(3,:,:,:)/(rhoce+weight_v)
    cume=cumsum3(rhoce)
    overhead_tile=max(overhead_tile,cume(nt+2*ncb,nt+2*ncb,nt+2*ncb)/real(np_tile_max_nu))
    if (cume(nt+2*ncb,nt+2*ncb,nt+2*ncb)>np_tile_max_nu) then
      print*, '  error: too many particles in this tile+buffer'
      print*, '  ',cume(nt+2*ncb,nt+2*ncb,nt+2*ncb),'>',np_tile_max_nu
      print*, '  on',this_image(), itx,ity,itz
      print*, '  please set tile_buffer_nu larger'
      stop
    endif
    ! create x_new v_new for this local tile
    !if (head) print*,'    xvnew loop'
    !$omp paralleldo &
    !$omp& default(shared) &
    !$omp& private(k,j,i,nlast,np,l,ip,xq,vreal,deltax,g,idx)
!    !$omp& reduction(+:rholocal)
    do k=1-ncb,nt+ncb ! update particle
    do j=1-ncb,nt+ncb
    do i=1-ncb,nt+ncb
      nlast=cum_nu(i,j,k,itx,ity,itz)
      np=rhoc_nu(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        xq=((/i,j,k/)-1d0) + (int(xp_nu(:,ip)+ishift_nu,izipx_nu)+rshift_nu)*x_resolution_nu
        vreal=tan((pi*real(vp_nu(:,ip)))/real(nvbin_nu-1)) / (sqrt(pi/2)/(sigma_vi_nu*vrel_boost))
        vreal=vreal+vfield_nu(:,i,j,k,itx,ity,itz)
        deltax=(dt_mid*vreal)/ncell
        g=ceiling(xq+deltax)
        rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
        idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3))
        xp_new_nu(:,idx)=xp_nu(:,ip)+nint(dt_mid*vreal/(x_resolution_nu*ncell))
        vreal=vreal-vfield_new(:,g(1),g(2),g(3))
        vp_new_nu(:,idx)=nint(real(nvbin_nu-1)*atan(sqrt(pi/2)/(sigma_vi_nu*vrel_boost)*vreal)/pi,kind=izipv_nu)
#       ifdef EID
          pid_new_nu(idx)=pid_nu(ip)
#       endif
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    sync all
    ! delete particles
    !if (head) print*,'    delete_particle loop'
    do k=1,nt
    do j=1,nt
      ileft=iright+1
      nlast=cume(nt,j,k)
      nlen=nlast-cume(0,j,k)
      iright=ileft+nlen-1
      xp_nu(:,ileft:iright)=xp_new_nu(:,nlast-nlen+1:nlast)
      vp_nu(:,ileft:iright)=vp_new_nu(:,nlast-nlen+1:nlast)
#     ifdef EID
        pid_nu(ileft:iright)=pid_new_nu(nlast-nlen+1:nlast)
#     endif
    enddo
    enddo
    ! update rhoc
    rhoc_nu(:,:,:,itx,ity,itz)=0
    rhoc_nu(1:nt,1:nt,1:nt,itx,ity,itz)=rhoce(1:nt,1:nt,1:nt)
    vfield_nu(:,:,:,:,itx,ity,itz)=0
    vfield_nu(:,1:nt,1:nt,1:nt,itx,ity,itz)=vfield_new(:,1:nt,1:nt,1:nt)
  enddo
  enddo
  enddo ! end looping over tiles
  nplocal_nu=iright
  xp_nu(:,nplocal_nu+1:)=0
  vp_nu(:,nplocal_nu+1:)=0
# ifdef PID
    pid_nu(nplocal_nu+1:)=0
# endif
  !nplocal_nu=sum(rhoc_nu(1:nt,1:nt,1:nt,:,:,:))
  !print*, iright,nplocal_nu
  !stop
  sync all

  ! calculate std of the velocity field
  cum_nu=cumsum6(rhoc_nu)
  std_vsim_nu=0; std_vsim_c_nu=0; std_vsim_res_nu=0
  !$omp paralleldo &
  !$omp& default(shared) &
  !$omp& private(itz,ity,itx,k,j,i,nlast,np,l,ip,vreal) &
  !$omp& reduction(+:std_vsim_c_nu,std_vsim_nu,std_vsim_res_nu)
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do k=1,nt
    do j=1,nt
    do i=1,nt
      std_vsim_c_nu=std_vsim_c_nu+sum(vfield_nu(:,i,j,k,itx,ity,itz)**2)
      nlast=cum_nu(i,j,k,itx,ity,itz)
      np=rhoc_nu(i,j,k,itx,ity,itz)
      do l=1,np
        ip=nlast-np+l
        vreal=tan(pi*real(vp_nu(:,ip))/real(nvbin_nu-1))/(sqrt(pi/2)/(sigma_vi_nu*vrel_boost))
        std_vsim_res_nu=std_vsim_res_nu+sum(vreal**2)
        vreal=vreal+vfield_nu(:,i,j,k,itx,ity,itz)
        std_vsim_nu=std_vsim_nu+sum(vreal**2)
      enddo
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  !$omp endparalleldo
  sync all

  ! co_sum
  if (head) then
    do i=2,nn**3
      std_vsim_c_nu=std_vsim_c_nu+std_vsim_c_nu[i]
      std_vsim_res_nu=std_vsim_res_nu+std_vsim_res_nu[i]
      std_vsim_nu=std_vsim_nu+std_vsim_nu[i]
    enddo
  endif
  sync all

  ! broadcast
  std_vsim_c_nu=std_vsim_c_nu[1]
  std_vsim_res_nu=std_vsim_res_nu[1]
  std_vsim_nu=std_vsim_nu[1]
  sync all

  ! divide
  std_vsim_nu=sqrt(std_vsim_nu/npglobal_nu)
  std_vsim_c_nu=sqrt(std_vsim_c_nu/nc/nc/nc/nn/nn/nn)
  std_vsim_res_nu=sqrt(std_vsim_res_nu/npglobal_nu)

  ! set sigma_vi_new according to particle statistics
  sigma_vi_new_nu=std_vsim_res_nu/sqrt(3.)
  if (head) then
    print*,'  std_vsim_nu    ',std_vsim_nu*sim%vsim2phys,'km/s'
    print*,'  std_vsim_c_nu  ',std_vsim_c_nu*sim%vsim2phys,'km/s'
    print*,'  std_vsim_res_nu',std_vsim_res_nu*sim%vsim2phys,'km/s'
    print*,'  sigma_vi_nu    ',sigma_vi_nu,'(simulation unit)'
    print*,'  sigma_vi_new_nu',sigma_vi_new_nu,'(simulation unit)'
    write(77) a-da,real((/std_vsim_nu,std_vsim_c_nu,std_vsim_res_nu/))
  endif
  sync all

  do i=1,nn**3
    overhead_tile=max(overhead_tile,overhead_tile[i])
  enddo

  if (head) then
    print*, '  tile overhead (nu)',overhead_tile*100,'% full'
    print*, '  comsumed tile_buffer_nu =',overhead_tile*tile_buffer_nu,'/',tile_buffer_nu
    print*, '  clean buffer of rhoc_nu'
  endif

  if (head) then
    npcheck=0
    do i=1,nn**3
      npcheck=npcheck+nplocal_nu[i]
    enddo
    print*, '  npcheck,npglobal=', npcheck,npglobal_nu
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*,''
    sync all
  endif

endsubroutine update_xp_nu

endmodule
