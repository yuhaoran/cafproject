program THREE_POINT_CORRELATION_FAST
  implicit none

  integer,parameter :: ngal=1000
  integer,parameter :: nnode=2*ngal
  real(kind=8),parameter :: accuracy=0.5
  real(kind=8),parameter :: dstep=2**0.5,dmin=1E-1
    ! for testing purpose: due to small ngal
  !real(kind=8),parameter :: dstep=2**0.1,dmin=10.
  real(kind=8),parameter :: pi=3.14159265358979, &
                            e=2.718281828,eps2=1E-8
  integer,dimension(2,nnode)::inode
  integer,dimension(ngal) :: ll
  real(kind=8),dimension(ngal) :: rpos
  real(kind=8),dimension(7,3) :: data
  real(kind=8),dimension(3) :: gamma1,gamma2
  real(kind=8),dimension(7,nnode) :: xyr
  !xyr: 1-x,2-y,3-g1,4-g2,5-kappa,6-noise,7-max_radius_square
  real(kind=8),dimension(:,:,:,:),allocatable :: cor
  real(kind=8),dimension(:,:),allocatable :: cor2

  character(len=80),parameter :: fn_in= &
                    '/tmp/pet_fat/catalog_gau_sorted.dat'
  !                  '/mnt/r0/pet_fat/F14.dat'
  character(len=80),parameter :: fn_cor3fast= &
                    '/tmp/pet_fat/cor3fast_n1000_t0.5.raw'
  !                  '/mnt/r0/pet_fat/F14cor3fast.dat'

  integer :: n,i,j,k,b,c,q,p,np,p1,p2, &
             idstep,idstep1,idstep2,idstep3,iglobal,ndstep
  integer(kind=8) :: s,t
  real(kind=8) :: rimin,rimax,rmin,rmax,xmax,ymax, &
                  y0,r1,r2,r3,d1,d2,d3,tmpx,tmpy, &
                  t1,t21,t22,t3,t4
  real(kind=8) :: theta2,alpha2,phi2,gamma,kappa,w,weight

  open(10,file=fn_in,form="formatted")
  read(10,*) xyr(:6,:ngal)
  close(10)
  !call Generate
  !xyr(6,:ngal)=1.
  xyr(6,:ngal)=1./xyr(6,:ngal)**2

  xmax=maxval(xyr(1,:ngal))
  ymax=maxval(xyr(2,:ngal))
  rmax=sqrt(xmax**2+ymax**2)
  ndstep=log(rmax/dmin)/log(dstep)+1

  Allocate(cor(10,ndstep,ndstep,ndstep))
  !cor: 1-g1,2-g2,3-g3,...,8-g8,9-kappa,10-Ntriplets
  Allocate(cor2(4,ndstep))
  !cor2: 1-r,2-sum,3-difference,4-Npairs
  !length r: max -> min in order of index
  cor=0.
  cor2=0.
  !call dscale
  call CPU_time(t1)
  print*,'time1: initialization',t1

  ll=(/ (i,i=1,ngal) /)
  iglobal=ngal
  call buildtree(ll,ngal)
  xyr(7,:ngal)=1E-5
  !do i=1,nnode
  !  write(*,'(F5.0,3i5)'),xyr(6,i),i,inode(:,i)
  !enddo
  call CPU_time(t21)
  print*,'time2.1: buildtree',t21-t1

  p=0
  q=0
  c=0
  b=0
  s=0
  t=0

  iglobal=ngal+1
  call Subdivide(iglobal,iglobal,iglobal,0)
  call CPU_time(t22)
  print*,'time2.2: raw 3PCF computation',t22-t21

  do k=1,ndstep
  do j=1,ndstep
  do i=1,ndstep
    if(cor(10,i,j,k)>0.) then
      cor(:9,i,j,k)=cor(:9,i,j,k)/spread(cor(10,i,j,k),1,9)
      weight=weight+cor(10,i,j,k)
      cor(10,i,j,k)=1/sqrt(cor(10,i,j,k))
    endif
  enddo
  enddo
  enddo
  call CPU_time(t3)
  print*,'time3: normalize 3PCF',t3-t22

  open(12,file=fn_cor3fast,form='binary')
  write(12) cor(:,:,:,:)
  close(12)
  !do i=20,ndstep-30
  !  if(cor(10,i,70,70)>0.) then
  !    write(*,'(i2,4E12.4)') i,cor((/1,2,3,10/),i,70,70)
  !  endif
  !enddo
  call CPU_time(t4)
  print*,'time4: output 4D array cor',t4-t3

  print*,''
  print*,'CPU_time=',t4-t1
  print*,'accuracy=',accuracy
  print*,'ngal=',ngal
  print*,'rmax=',rmax
  print*,'ndstep=',ndstep !catalog_ran.dat: >=129; F14.dat: 113
  print*,'noise=',sqrt(1/weight)
  print*,'SUBDIVIDE=',t
  print*,'3PCF(triplets)=',s
  print*,'overlap=',c
  print*,'out of range=',b
  print*,''

contains
  recursive subroutine Subdivide(node1,node2,node3,counter)
    integer node1,node2,node3
    integer counter
    if(node1/=node2 .and. node3==node1) print*,'error'
    if(node1/=node2 .and. node3==node2) print*,'error'
t=t+1
  if(node1==node3) then
    if(node1>ngal) then
      call Subdivide(inode(1,node1),inode(1,node1),inode(1,node1),0)
      call Subdivide(inode(2,node1),inode(2,node1),inode(2,node1),0)
      call Subdivide(inode(1,node1),inode(1,node1),inode(2,node1),0)
      call Subdivide(inode(2,node1),inode(2,node1),inode(1,node1),0)
    endif
  elseif(node1==node2) then
    if(node1>ngal) then
      call Subdivide(inode(1,node1),inode(1,node1),node3,0)
      call Subdivide(inode(2,node1),inode(2,node1),node3,0)
      call Subdivide(inode(1,node1),inode(2,node1),node3,0)
    endif
  else
    data(:,1)=xyr(:,node1)
    data(:,2)=xyr(:,node2)
    data(:,3)=xyr(:,node3)
    d3=(data(1,1)-data(1,2))**2+(data(2,1)-data(2,2))**2
    d2=(data(1,1)-data(1,3))**2+(data(2,1)-data(2,3))**2
    d1=(data(1,2)-data(1,3))**2+(data(2,2)-data(2,3))**2
    r1=sqrt(d1)
    r2=sqrt(d2)
    r3=sqrt(d3)
    if(min(r1,r2,r3)<=0.) then
      print*,'overlapping nodes are', node1,node2,node3
      c=c+1  !overlap
      return
    endif
    idstep1=log(rmax/r1)/log(dstep)+1
    idstep2=log(rmax/r2)/log(dstep)+1
    idstep3=log(rmax/r3)/log(dstep)+1

      if( data(7,1)/min(d2,d3)<=accuracy**2 ) then
        if (counter==2) then
          if(min(idstep1,idstep2,idstep3)>=1 .and. &
             max(idstep1,idstep2,idstep3)<=ndstep) then
            call three_point(node1,node2,node3)
          else
             b=b+1 !out of range
          endif
          return
        else
          call Subdivide(node2,node3,node1,counter+1)
        endif
      else
        if (node1<=ngal) then
          if (counter==2) then
            if(min(idstep1,idstep2,idstep3)>=1 .and. &
               max(idstep1,idstep2,idstep3)<=ndstep) then
              call three_point(node1,node2,node3)
            else
              b=b+1 !out of range
            endif
            return
          else
            call Subdivide(node2,node3,node1,counter+1)
          endif
        else
          call Subdivide(inode(1,node1),node2,node3,0)
          call Subdivide(inode(2,node1),node2,node3,0)
        endif
      endif

  endif
  end subroutine Subdivide

  subroutine three_point(node1,node2,node3)
!    real ttt
    integer node1,node2,node3
    s=s+1
      data(:,1)=xyr(:,node1)
      data(:,2)=xyr(:,node2)
      data(:,3)=xyr(:,node3)
      tmpx=data(1,1)
      tmpy=data(2,1)
      np=0
      p1=0
      p2=0
      p=0
      do p=2,3
        if(data(1,p)<tmpx) then
           tmpx=data(1,p)
           p1=np
           np=p
        elseif(data(1,p)==tmpx .and. data(2,p)<tmpy) then
           tmpy=data(2,p)
           p1=np
           np=p
        else
           p2=p
        endif
      enddo
      if(np>1) then
        if(max(p1,p2)<1) write(*,*) p1,p2,np
        data=data(:,(/np,max(p1,p2),1/))
      endif
      if(abs(data(1,2)-data(1,1))<eps2 .or. &
         abs(data(1,3)-data(1,1))<eps2) then
        if(data(1,2)<data(1,3)) data=data(:,(/1,3,2/))
      elseif(data(1,2)/=data(1,1)) then
        y0=(data(2,2)-data(2,1))/ &
        (data(1,2)-data(1,1))*(data(1,3)-data(1,1))
        if((data(2,3)-data(2,1))<y0) data=data(:,(/1,3,2/))
      endif

      r3=sqrt((data(1,1)-data(1,2))**2+(data(2,1)-data(2,2))**2)
      r2=sqrt((data(1,1)-data(1,3))**2+(data(2,1)-data(2,3))**2)
      r1=sqrt((data(1,2)-data(1,3))**2+(data(2,2)-data(2,3))**2)

      idstep1=log(rmax/r1)/log(dstep)+1
      idstep2=log(rmax/r2)/log(dstep)+1
      idstep3=log(rmax/r3)/log(dstep)+1

      do i=1,3
        gamma1(i)=data(3,i)
        gamma2(i)=data(4,i)
        theta2=atan2(gamma2(i),gamma1(i))
        phi2=2*atan2( data(2,2)-data(2,1), data(1,2)-data(1,1) )
        alpha2=theta2-phi2
        gamma=sqrt( gamma1(i)**2 + gamma2(i)**2 )
        gamma1(i)=gamma*cos(alpha2)
        gamma2(i)=gamma*sin(alpha2)
      enddo
      kappa=data(5,1)*data(5,2)*data(5,3)
      w=data(6,1)*data(6,2)*data(6,3)

!      call random_number(ttt)
!      if (ttt<1E-5) write(*,*) min(sqrt(data(7,1))/min(r1,r2),sqrt(data(7,2))/min(r1,r3),sqrt(data(7,3))/min(r2,r3)), max(r1,r2,r3),min(r1,r2,r3)
      cor(1,idstep1,idstep2,idstep3) =cor(1,idstep1,idstep2,idstep3)+gamma1(1)*gamma1(2)*gamma1(3)*w
      cor(2,idstep1,idstep2,idstep3) =cor(2,idstep1,idstep2,idstep3)+gamma1(1)*gamma1(2)*gamma2(3)*w
      cor(3,idstep1,idstep2,idstep3) =cor(3,idstep1,idstep2,idstep3)+gamma1(1)*gamma2(2)*gamma1(3)*w
      cor(4,idstep1,idstep2,idstep3) =cor(4,idstep1,idstep2,idstep3)+gamma2(1)*gamma1(2)*gamma1(3)*w
      cor(5,idstep1,idstep2,idstep3) =cor(5,idstep1,idstep2,idstep3)+gamma1(1)*gamma2(2)*gamma2(3)*w
      cor(6,idstep1,idstep2,idstep3) =cor(6,idstep1,idstep2,idstep3)+gamma2(1)*gamma2(2)*gamma1(3)*w
      cor(7,idstep1,idstep2,idstep3) =cor(7,idstep1,idstep2,idstep3)+gamma2(1)*gamma1(2)*gamma2(3)*w
      cor(8,idstep1,idstep2,idstep3) =cor(8,idstep1,idstep2,idstep3)+gamma2(1)*gamma2(2)*gamma2(3)*w
      cor(9,idstep1,idstep2,idstep3) =cor(9,idstep1,idstep2,idstep3)+kappa*w
      cor(10,idstep1,idstep2,idstep3) =cor(10,idstep1,idstep2,idstep3)+w

  end subroutine three_point

  subroutine dscale
    rimin=rmax
    idstep=0
    do while(rimin>dmin) !i=1,ndstep
       idstep=idstep+1
       rimax=rimin
       rimin=rmax*dstep**(-idstep)
       cor2(1,idstep)=sqrt(rimax*rimin)
    enddo
    rmin=rimin
  end subroutine dscale

  recursive subroutine buildtree(ll,n)
    integer n1,n2,nl,nh,i,ig,n,itmp
    integer ll(n)
    integer ifar(1)
    real w1,w2,dx(2),cx(2),rm1,rm2
    real eps,rpivot

    iglobal=iglobal+1
    if (ngal .eq. 1) then
       xyr(:5,iglobal)=xyr(:5,1)
       xyr(7,iglobal)=0
       inode(:2,iglobal)=-1
       return
    end if

    xyr(6,iglobal)=sum(xyr(6,ll))
    xyr(:5,iglobal)=sum(xyr(:5,ll)*spread(xyr(6,ll),1,5),2) &
                    /spread(xyr(6,iglobal),1,5)
    cx=xyr(:2,iglobal)
    ifar=maxloc((cx(1)-xyr(1,ll))**2+(cx(2)-xyr(2,ll))**2)
    dx=cx-xyr(:2,ll(ifar(1)))
    if (sum(dx**2) .eq. 0) call random_number(dx)
    dx=dx/sqrt(sum(dx**2))
    rpos(ll)=dx(1)*xyr(1,ll)+dx(2)*xyr(2,ll)
    xyr(7,iglobal)=(cx(1)-xyr(1,ll(ifar(1))))**2 &
                   +(cx(2)-xyr(2,ll(ifar(1))))**2

    n1=1
    n2=n
    eps=1.e-8  ! needed to break identical galaxy positions
    do
       w1=((n+1)/2.-n1)/(n2-n1)
       w2=(n2-(n+1)/2.)/(n2-n1)
       if (w1 < .1/n) then
          w1=.1/(n2-n1)
          w2=1-w1
       endif
       if (w2 < .1/n) then
          w2=.1/(n2-n1)
          w1=1-w2
       endif
       rpivot=minval(rpos(ll(n1:n2)))*w1+maxval(rpos(ll(n1:n2)))*w2
       nl=n1
       nh=n2
       do
          do
             if (rpos(ll(nl))>rpivot+eps) exit
             nl=nl+1
             if (nl .eq. n2) exit
             ! happens if two positions are identical
          enddo
          do
             if (rpos(ll(nh))<rpivot-eps) exit
             nh=nh-1
             if (nh .eq. n1) exit
          enddo
          if (nh<=nl) exit
          itmp=ll(nh)
          ll(nh)=ll(nl)
          ll(nl)=itmp
       end do
! here, nl=nh+1
       if (2*nh-1>=n) then
          n2=nh
       else if (2*nl-1 <= n) then
          n1=nl
       else
          if (mod (n,2) .eq. 1) then
             n1=(n+1)/2
             n2=n2
          else
             n1=n/2
             n2=n/2+1
          end if
       end if
       if (n2-n1 < 1) exit
       if (n2-n1 < 2 .and. mod(n,2) .eq. 0) exit
    end do
! if n is even, n1=n/2, n2=n/2+1, else n1=n2=(n+1)/2
    if (mod(n,2) .eq. 1) then
       rm1=sum(rpos(ll(:n1-1)))/(n1-1)
       rm2=sum(rpos(ll(n1+1:)))/(n-n1)
       if (rpos(ll(n1))-rm1 < rm2-rpos(ll(n1))) then
          n2=n2+1
       else
          n1=n1-1
       end if
    end if
    ig=iglobal
!    if (mod(iglobal,100) .eq. 0) write(*,*) iglobal,n1,n2,n
!    if (iglobal <= ngal+20 .and. iglobal>=ngal+1) write(*,*) iglobal,n1,n2,n
    if (n1>1) then
       inode(1,ig)=iglobal+1
       call buildtree(ll(:n1),n1)
    else
       inode(1,ig)=ll(n1)
    end if
    if (n2<n) then
       inode(2,ig)=iglobal+1
       call buildtree(ll(n2:n),n-n2+1)
    else
       inode(2,ig)=ll(n2)
    end if
  end subroutine buildtree

end program
