  integer, parameter :: n=1000,neq=2
  real, dimension(neq,n) :: u

  call init

  cfl=0.8
  open(10,file='u0.dat',access='stream')
  open(11,file='u1.dat',access='stream')
  t=0
  do nt=1,1000
     do i=1,n/50
        dt=relaxing(u)
        t=t+dt
     end do
     write(*,*) t
     write(10) u(1,:)
     write(11) u(2,:)
  end do
  
contains
  
  subroutine init
    pi=4*atan(1.)
    u(1,:)=1+0.2*(/(sin(2*pi*(i-0.5)/n),i=1,n)/)
    u(2,:)=0
  end subroutine init
  
  real function relaxing(u)
    implicit none
    real cmax
    real, dimension(n) :: c
    real, dimension(neq,n) :: u,u1,w,fu,fr,fl,dfl,dfr

    cmax=averageflux(u,w,c)

    dt=cfl/cmax
    fr=(u*spread(c,1,neq)+w)/2
    fl=cshift(u*spread(c,1,neq)-w,1,2)/2
    fu=fr-fl
    u1=u-(fu-cshift(fu,-1,2))*dt/2

    cmax=averageflux(u1,w,c)
    fr=(u1*spread(c,1,neq)+w)/2
    dfl=(fr-cshift(fr,-1,2))/2
    dfr=cshift(dfl,1,2)
    call vanleer(fr,dfl,dfr)

    fl=cshift(u1*spread(c,1,neq)-w,1,2)/2
    dfl=(cshift(fl,-1,2)-fl)/2
    dfr=cshift(dfl,1,2)
    call vanleer(fl,dfl,dfr)

    fu=fr-fl
    u=u-(fu-cshift(fu,-1,2))*dt
    relaxing=dt
    return
  end function relaxing
  real function averageflux(u,w,c)
    implicit none
    integer i
    real, dimension(n)::p,v,w(neq,n),c(n)
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
    
    
    
    v=3*u(2,:)/u(1,:)/4
    c=(1/sqrt(3.)+abs(v))!*1.1
    averageflux=maxval(c)
    p=u(1,:)*(1-4*v**2/3)/3
    w(1,:)=u(2,:)
!    w(2,:)=u(2,:)**2/u(1,:)/2+u(1,:)/3
    w(2,:)=(5*u(1,:)-4*sqrt(u(1,:)**2-3*u(2,:)**2/4))/3
    
    return

  end function averageflux
  
  subroutine vanleer(f,a,b)
    implicit none
    real, dimension(neq,n) :: f,a,b,c
    
    c=a*b
    where(c>0) f=f+2*c/(a+b)
    return
  end subroutine vanleer
  subroutine minmod(f,a,b)
    implicit none
    real, dimension(neq,n) :: f,a,b

    f=f+(sign(1.,a)+sign(1.,b))*min(abs(a),abs(b))/2.
    return
  end subroutine minmod

  
end program
