!#define debug
!#define redundant
!#define randomv
subroutine update_particle
use variables
implicit none
save

! x(1:3,:)
integer ileft,iright,nlast,nlen,idx
integer g(3),np ! index of grid
integer rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! double buffer tile
integer cume(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
integer rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! count writing
real vrand(3), v_th
integer ii,jj,kk
integer(8) irand(3)
integer iii(3)
integer pr1,pr2,pr3,pr4

if (head) print*,'update_particle'

#ifdef randomv
  v_th=1.5/real(2**(izipx*8))*ncell*0.5

  pr1=71
  pr2=73
  pr3=79
  pr4=83
  iii=(/1,3,7/)
#endif
dt_mid=(dt_old+dt)/2


overhead_tile=0
iright=0
do itz=1,nnt ! loop over tile
do ity=1,nnt
do itx=1,nnt
  rhoce=0
  rholocal=0
  x_new=0
  v_new=0
# ifdef PID
  pid_new=0
# endif

  do k=1-ncb,nt+ncb ! loop over coarse grid
  do j=1-ncb,nt+ncb
  do i=1-ncb,nt+ncb
    nlast=cum(i,j,k,itx,ity,itz)
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np
      ip=nlast-np+l
#ifdef randomv
      ii=i+(itx-1)*nt
      jj=j+(ity-1)*nt
      kk=k+(itz-1)*nt
      irand=(abs(ii-jj)+1)*(abs(jj-kk)+1)*ii*jj*kk*(nc-ii+1)*(nc-jj+1)*(nc-kk+1)+iii+its+l
      irand=mod(irand,pr1)
      irand=mod(irand**2+ii*pr1,pr2)
      irand=mod(irand**2+jj*pr2,pr3)
      irand=mod(irand**2+kk*pr3,pr4)
      vrand=(irand+0.5)/pr4-0.5
      vrand=vrand*(1/(1+int(abs(dt_mid*v(:,ip)*v_i2r/v_th))))
#else
      vrand=0
#endif
      !print*, 'x0', ( 4.*nt*((/itx,ity,itz/)-1) + 4.*((/i,j,k/)-1) + 4*((x(:,ip)+int(-128,1))+128.5)/256 )
      !print*, 'dx', dt_mid*v(:,ip)*v_i2r(:)
      !g=ceiling( &
      !  ( (/i,j,k/)-1 + ((x(:,ip)+int(-128,1))+128.5)/256. ) &
      !  + dt_mid*v(:,ip)*v_i2r(:)/4. &
      !  )
      xq=(/i,j,k/)-1d0 + ((x(:,ip)+ishift)+rshift)*x_resolution
      deltax=dt_mid*v(:,ip)*v_i2r/4+(x_resolution*ncell)*vrand
      g=ceiling(xq+deltax)
      !g=ceiling(((/i,j,k/)-1+(x(:,ip)+ishift+rshift)*x_resolution)+dt_mid*v(:,ip)*v_i2r/4 &
      !          +(x_resolution*ncell)*vrand)
      !print*, 'particle',ip,itx,ity,itz,g
      !print*, dt_mid*v(:,ip)*v_i2r(:)/4.
      !print*, g,vrand
      rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1 ! update mesh
    enddo
  enddo
  enddo
  enddo

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

  do k=1-ncb,nt+ncb ! update particle
  do j=1-ncb,nt+ncb
  do i=1-ncb,nt+ncb
    nlast=cum(i,j,k,itx,ity,itz)
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np
      ip=nlast-np+l
#ifdef randomv
      ii=i+(itx-1)*nt
      jj=j+(ity-1)*nt
      kk=k+(itz-1)*nt
      irand=(abs(ii-jj)+1)*(abs(jj-kk)+1)*ii*jj*kk*(nc-ii+1)*(nc-jj+1)*(nc-kk+1)+iii+its+l
      irand=mod(irand,pr1)
      irand=mod(irand**2+ii*pr1,pr2)
      irand=mod(irand**2+jj*pr2,pr3)
      irand=mod(irand**2+kk*pr3,pr4)
      vrand=(irand+0.5)/pr4-0.5
      vrand=vrand*(1/(1+int(abs(dt_mid*v(:,ip)*v_i2r/v_th))))
#else
      vrand=0
#endif
      xq=(/i,j,k/)-1d0 + ((x(:,ip)+ishift)+rshift)*x_resolution
      deltax=dt_mid*v(:,ip)*v_i2r/4+(x_resolution*ncell)*vrand
      g=ceiling(xq+deltax)
      !g=ceiling(((/i,j,k/)-1+(x(:,ip)+ishift+rshift)*x_resolution)+dt_mid*v(:,ip)*v_i2r(:)/4 &
      !          +(x_resolution*ncell)*vrand)
      rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
      idx=cume(g(1),g(2),g(3))-rhoce(g(1),g(2),g(3))+rholocal(g(1),g(2),g(3)) ! index for writing

#ifdef debug
      x_new(:,idx)=x_new(:,idx)+1
#else
      x_new(:,idx)=x(:,ip)+nint(dt_mid*v(:,ip)*v_i2r(:)/(x_resolution*ncell) &
                                + vrand)
#endif
			v_new(:,idx)=v(:,ip)
#ifdef PID
      pid_new(:,idx)=pid(:,ip)
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

  do k=1,nt
  do j=1,nt
    ileft=iright+1
    nlast=cume(nt,j,k)
    nlen=nlast-cume(0,j,k)
    iright=ileft+nlen-1
    x(:,ileft:iright)=x_new(:,nlast-nlen+1:nlast)
    v(:,ileft:iright)=v_new(:,nlast-nlen+1:nlast)
#ifdef PID
    pid(:,ileft:iright)=pid_new(:,nlast-nlen+1:nlast)
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
x(:,nplocal+1:)=0
v(:,nplocal+1:)=0
#ifdef PID
pid(:,nplocal+1:)=0
#endif

sync all

if (head) then
  npcheck=0
  do i=1,nn**3
    npcheck=npcheck+nplocal[i]
  enddo
  print*, '  npcheck,npglobal=', npcheck,npglobal
endif

sync all

!print*, 'nplocal, v_i2r =', iright, v_i2r
!print*, 'sum of x', sum(x(:,:nplocal)*unit8)
!print*, 'sum of v', sum(v(:,:nplocal)*unit8)
!print*, 'sum of rhoc', sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
!print*, x(:,1), v(:,1)
!print*, 'update_particle done'

endsubroutine update_particle
