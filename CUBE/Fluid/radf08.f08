  integer, parameter :: n=64,nx=n,ny=n,nz=n,neq=4,nnx=1,nny=1
  real, dimension(neq,nx,ny,nz) :: u[nnx,nny,*]
  real dt[*]
  character*5 fn5
  integer myi(1),myi3(3)

  myi=this_image(dt)
  write(fn5,'(i5.5)') myi
  nnz=num_images()/(nnx*nny)
  write(*,*) 'dim=',nnx,nny,nnz
  call init
  cfl=0.8
  open(10,file='u0.'//fn5//'.dat',access='stream')
  open(11,file='udiv.'//fn5//'.dat',access='stream')
  open(12,file='ucrl.'//fn5//'.dat',access='stream')
  t=0
  do it=1,4000
     do i=1,n/40+1
        call boundary
        dt=timestep()
        call co_min(dt)
        call sweepx
        call sweepy
        call sweepz
        call sweepz
        call sweepy
        call sweepx
        t=t+dt
     end do
     write(*,*) t
     write(10) u(1,:,ny/2,:)
     call writecurl
     if (mod(it,10) .eq. 1) call pspect
!     write(11) u(2,:,:,:)
  end do
  
contains

  subroutine co_min(x)
    implicit none
    real x[*]
    call co_min1(x,1)
  end subroutine co_min
  recursive subroutine co_min1(x,istride)
    implicit none
    real x[*],xtmp
    integer istride,itarget,myimage(1),idir
    myimage=this_image(x)
    if (istride >= num_images()) return
    idir=2*mod(myimage(1)/istride,2)-1
    itarget=mod(myimage(1)+2*num_images()+idir*istride-1,num_images())+1
    xtmp=max(x,x[itarget])
    sync all
    x=xtmp
    call co_min1(x,2*istride)
    return    
  end subroutine co_min1
  subroutine boundary
    implicit none
    integer,dimension(3):: mycoord,mm,mp
    integer, parameter :: nbuf=6
    
    mycoord=this_image(u)
    mp=mod(mycoord,(/nnx,nny,nnz/))+1
    mm=mod(mycoord+(/nnx,nny,nnz/)-2,(/nnx,nny,nnz/))+1
    u(:,:,:,nz-nbuf+1:)=u(:,:,:,nbuf+1:2*nbuf)[mycoord(1),mycoord(2),mp(3)]
    u(:,:,:,:nbuf)=u(:,:,:,nz-2*nbuf+1:nz-nbuf)[mycoord(1),mycoord(2),mm(3)]
    sync all
    u(:,:,ny-nbuf+1:,:)=u(:,:,nbuf+1:2*nbuf,:)[mycoord(1),mp(2),mycoord(3)]
    u(:,:,:nbuf,:)=u(:,:,ny-2*nbuf+1:ny-nbuf,:)[mycoord(1),mm(2),mycoord(3)]
    sync all
    u(:,nx-nbuf+1:,:,:)=u(:,nbuf+1:2*nbuf,:,:)[mp(1),mycoord(2),mycoord(3)]
    u(:,:nbuf,:,:)=u(:,nx-2*nbuf+1:nx-nbuf,:,:)[mm(1),mycoord(2),mycoord(3)]
    sync all    
  end subroutine boundary
  
  subroutine writecurl
  j=ny/2
  jp=j+1
  do k=1,nz
    kp=mod(k,nz)+1
    do i=1,nx
       ip=mod(i,nx)+1
       divx=(u(2,ip,j,k)-u(2,i,j,k)+u(2,ip,jp,k)-u(2,i,jp,k) + &
           u(2,ip,j,kp)-u(2,i,j,kp)+u(2,ip,jp,kp)-u(2,i,jp,kp) )/4
       divy=(u(3,i,jp,k)-u(3,i,j,k)+u(3,ip,jp,k)-u(3,ip,j,k) + &
           u(3,i,jp,kp)-u(3,i,j,kp)+u(3,ip,jp,kp)-u(3,ip,j,kp) )/4
       divz=(u(4,i,j,kp)-u(4,i,j,k)+u(4,ip,j,kp)-u(4,ip,j,k) + &
           u(4,i,jp,kp)-u(4,i,jp,k)+u(4,ip,jp,kp)-u(4,ip,jp,k) )/4
       write(11) divx+divy+divz
       vzx=(u(4,ip,j,k)-u(4,i,j,k)+u(4,ip,jp,k)-u(4,i,jp,k) + &
           u(4,ip,j,kp)-u(4,i,j,kp)+u(4,ip,jp,kp)-u(4,i,jp,kp) )/4
       vyx=(u(3,ip,j,k)-u(3,i,j,k)+u(3,ip,jp,k)-u(3,i,jp,k) + &
           u(3,ip,j,kp)-u(3,i,j,kp)+u(3,ip,jp,kp)-u(3,i,jp,kp) )/4
       vxy=(u(2,i,jp,k)-u(2,i,j,k)+u(2,ip,jp,k)-u(2,ip,j,k) + &
           u(2,i,jp,kp)-u(2,i,j,kp)+u(2,ip,jp,kp)-u(2,ip,j,kp) )/4
       vzy=(u(4,i,jp,k)-u(4,i,j,k)+u(4,ip,jp,k)-u(4,ip,j,k) + &
           u(4,i,jp,kp)-u(4,i,j,kp)+u(4,ip,jp,kp)-u(4,ip,j,kp) )/4
       vxz=(u(2,i,j,kp)-u(2,i,j,k)+u(2,ip,j,kp)-u(2,ip,j,k) + &
           u(2,i,jp,kp)-u(2,i,jp,k)+u(2,ip,jp,kp)-u(2,ip,jp,k) )/4
       vyz=(u(3,i,j,kp)-u(3,i,j,k)+u(3,ip,j,kp)-u(3,ip,j,k) + &
           u(3,i,jp,kp)-u(3,i,jp,k)+u(3,ip,jp,kp)-u(3,ip,jp,k) )/4
       curlx=vzy-vyz
       curly=vxz-vzx
       curlz=vyx-vxy
       write(12) curlx**2+curly**2+curlz**2
    enddo
  enddo
  end subroutine writecurl


  subroutine init
    implicit none
    integer,parameter :: nred=2
    real rhored(nx/nred,ny/nred,nz/nred)
    pi=4*atan(1.)

    open(100,file='uinit.'//fn5//'.dat',access='stream')
    read(100)  rhored
    u=0
    do i=1,nred
       do j=1,nred
          do k=1,nred
             u(1,i::nred,j::nred,k::nred)=rhored
          end do
       end do
    end do
    return
  end subroutine init

  real function timestep()
!$omp workshare
    timestep= cfl/maxval(1./sqrt(3.)+sqrt(sum(u(2:4,:,:,:)**2,1)))
!$omp end workshare
    return
  end function timestep
    
  subroutine sweepx
    implicit none
    integer j,k
    real u1d(neq,nx)

    !$omp parallel do default(shared) private(j,k,u1d)
    do k=1,nz
       do j=1,ny
          u1d=u(:,:,j,k)
          call relaxing(u1d,nx)
          u(:,:,j,k)=u1d
       end do
    end do
    return
  end subroutine sweepx

  subroutine sweepy
    implicit none
    integer i,k
    real u1d(neq,ny)

    !$omp parallel do default(shared) private(i,k,u1d)
    do k=1,nz
       do i=1,nx
          u1d((/1,3,2,4/),:)=u(:,i,:,k)
          call relaxing(u1d,ny)
          u(:,i,:,k) = u1d((/1,3,2,4/),:)
       end do
    end do
    return
  end subroutine sweepy
  subroutine sweepz
    implicit none
    integer i,j
    real u1d(neq,nz)
    !$omp parallel do default(shared) private(i,j,u1d)
    do j=1,ny
       do i=1,nx
          u1d((/1,4,3,2/),:)=u(:,i,j,:)
          call relaxing(u1d,nz)
          u(:,i,j,:)=u1d((/1,4,3,2/),:)
       end do
    end do
    return
  end subroutine sweepz
  
  recursive subroutine relaxing(u,n)
    implicit none
    integer n
    real cmax
    real, dimension(n) :: c
    real, dimension(neq,n) :: u,u1,w,fu,fr,fl,dfl,dfr

    cmax=averageflux(u,w,c,n)
    fr=(u*spread(c,1,neq)+w)/2
    fl=cshift(u*spread(c,1,neq)-w,1,2)/2
    fu=fr-fl
    u1=u-(fu-cshift(fu,-1,2))*dt/2

    cmax=averageflux(u1,w,c,n)
    fr=(u1*spread(c,1,neq)+w)/2
    dfl=(fr-cshift(fr,-1,2))/2
    dfr=cshift(dfl,1,2)
    call minmod(fr,dfl,dfr,n)

    fl=cshift(u1*spread(c,1,neq)-w,1,2)/2
    dfl=(cshift(fl,-1,2)-fl)/2
    dfr=cshift(dfl,1,2)
    call minmod(fl,dfl,dfr,n)

    fu=fr-fl
    u=u-(fu-cshift(fu,-1,2))*dt
    return
  end subroutine relaxing
  recursive real function averageflux(u,w,c,n)
    implicit none
    integer n
    real, dimension(n)::rho,rhog2,v2,v,w(neq,n),c(n)
    real u(neq,n)

    ! relativistic ideal fluid
    ! stress-energy tensor T=(\rho+p)u^\mu u^\nu + p g^\mu\nu
    ! metric convention g=diag(-1,1,1,1)
    ! Lorenz factor 1/\gamma^2=1-\beta^2
    ! 4-velocity u0=\gamma,u1=\gamma\beta, u_\mu u^\mu=-1
    ! equation of state: p=rho/3
    ! conservation law: T^00,0 + T^01,1=0 (continuity)
    ! T^01,0+T^11,1 = 0 (momentum)
    ! define: U0=T^00=(\rho+p)\gamma^2-p, U1=T^01 = (\rho+p)\gamma^2\beta
    ! simplify: U0=\rho(4\gamma^2-1)/3,U1=4\rho\gamma^2\beta/3
    ! U0=\rho(4/(1-\beta^2)-1)/3,U1=4\rho\beta/3/(1-\beta^2)
    ! U1/U0=4\beta/3+O(\beta^3)
    ! T11=\rho(4\gamma^2\beta^2+1)/3=U1^2/U0+p+O(beta^3)
    
    
    
    v2=3*sum(u(2:4,:)**2,1)/4
    rho=2*sqrt(u(1,:)**2-v2)-u(1,:)
    rhog2=rho/4+3*u(1,:)/4
    v=u(2,:)/rhog2
    c=(1/sqrt(3.)+abs(v))
    averageflux=maxval(c)
    w(1,:)=u(2,:)
    w(2,:)=3*u(2,:)**2/(3*u(1,:)+rho)+rho/3
    w(3,:)=3*u(2,:)*u(3,:)/(3*u(1,:)+rho)
    w(4,:)=3*u(2,:)*u(4,:)/(3*u(1,:)+rho)
    return

  end function averageflux
  
  subroutine vanleer(f,a,b,n)
    implicit none
    integer n
    real, dimension(neq,n) :: f,a,b,c
    
    c=a*b
    where(c>0) f=f+2*c/(a+b)
    return
  end subroutine vanleer
  recursive subroutine minmod(f,a,b,n)
    implicit none
    integer n
    real, dimension(neq,n) :: f,a,b

    f=f+(sign(1.,a)+sign(1.,b))*min(abs(a),abs(b))/2.
    return
  end subroutine minmod

  
end program
