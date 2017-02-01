!! module for power spectrum analysis
! uses pencil_fft however work for single image only
! cx1,cx2 can be memory optimized
! auto-power can be memory optimized
! check nexp frequently
#define linear_kbin

module powerspectrum
use pencil_fft

integer,parameter :: nk_ny=ng*nn/2 ! Nyquist wave number
#ifdef linear_kbin
  integer,parameter :: nbin=nint(nk_ny*sqrt(3.))
#else
  integer,parameter :: nbin=floor(4*log(nk_ny*sqrt(3.)/0.95)/log(2.))
#endif

contains

subroutine cross_power(xi,cube1,cube2)
implicit none

integer i,j,k,ig,jg,kg,ibin
real kr,kx(3),sincx,sincy,sincz,sinc,rbin

real cube1(ng,ng,ng),cube2(ng,ng,ng)
real xi(10,nbin)[*], xi_input(10,nbin)
real amp11,amp12,amp22
complex cx1(ng*nn/2+1,ng,npen),cx2(ng*nn/2+1,ng,npen)

real,parameter :: nexp=4.0 ! CIC kernel
  xi_input=xi
  xi=0

  r3=cube1
  call pencil_fft_forward
  cx1=cxyz

  r3=cube2
  call pencil_fft_forward
  cx2=cxyz

  xi=0
  do k=1,npen
  do j=1,ng
  do i=1,ng*nn/2+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ng+j
    ig=i
    kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    if ((ig==1.or.ig==ng*nn/2+1) .and. jg>ng*nn/2+1) cycle
    if ((ig==1.or.ig==ng*nn/2+1) .and. (jg==1.or.jg==ng*nn/2+1) .and. kg>ng*nn/2+1) cycle
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
    sincx=merge(1.0,sin(pi*kx(1)/ng/nn)/(pi*kx(1)/ng/nn),kx(1)==0.0)
    sincy=merge(1.0,sin(pi*kx(2)/ng/nn)/(pi*kx(2)/ng/nn),kx(2)==0.0)
    sincz=merge(1.0,sin(pi*kx(3)/ng/nn)/(pi*kx(3)/ng/nn),kx(3)==0.0)
    sinc=sincx*sincy*sincz
#ifdef linear_kbin
      ibin=nint(kr)
#else
      rbin=4.0/log(2.)*log(kr/0.95)
      ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
#endif
    xi(1,ibin)=xi(1,ibin)+1 ! number count
    xi(2,ibin)=xi(2,ibin)+kr ! k count
    amp11=real(cx1(i,j,k)*conjg(cx1(i,j,k)))/(ng**3)/(ng**3)/(sinc**0.0)*4*pi*kr**3 ! linear density field
    amp22=real(cx2(i,j,k)*conjg(cx2(i,j,k)))/(ng**3)/(ng**3)/(sinc**4.0)*4*pi*kr**3 ! CIC'ed density field
    amp12=real(cx1(i,j,k)*conjg(cx2(i,j,k)))/(ng**3)/(ng**3)/(sinc**2.0)*4*pi*kr**3 ! cross

    xi(3,ibin)=xi(3,ibin)+amp11 ! auto power 1
    xi(4,ibin)=xi(4,ibin)+amp22 ! auto power 2
    xi(5,ibin)=xi(5,ibin)+amp12 ! cross power
    xi(6,ibin)=xi(6,ibin)+1/sinc**2.0 ! kernel
    xi(7,ibin)=xi(7,ibin)+1/sinc**4.0 ! kernel

  enddo
  enddo
  enddo
  sync all

  if (head) then
    do i=2,nn**3
      xi=xi+xi(:,:)[i]
    enddo
    xi(2,:)=xi(2,:)/xi(1,:)*(2*pi)/box ! k_phy
    xi(3,:)=xi(3,:)/xi(1,:) ! Delta_LL
    xi(4,:)=xi(4,:)/xi(1,:) ! Delta_RR
    xi(5,:)=xi(5,:)/xi(1,:) ! Delta_LR ! cross power
    xi(6,:)=xi(6,:)/xi(1,:) ! kernel
    xi(7,:)=xi(7,:)/xi(1,:) ! kernel
    xi(8,:)=xi(5,:)/sqrt(xi(3,:)*xi(4,:)) ! r
    xi(9,:)=sqrt(xi(4,:)/xi(3,:)) ! b
    xi(10,:)=xi(8,:)**4/xi(9,:)**2 * xi(4,:) ! P_RR*r^4/b^2 reco power
  endif

  open(15,file='power_fields.dat',status='replace',access='stream')

  !! Wiener
  print*, 'Wiener filter delta_L'
  cxyz=0
  do k=1,npen
  do j=1,ng
  do i=1,ng*nn/2+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ng+j
    ig=i
    kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
#ifdef linear_kbin
      ibin=nint(kr)
#else
      rbin=4.0/log(2.)*log(kr/0.95)
      ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
#endif
    ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
    cxyz(i,j,k)=cx1(i,j,k)*xi(8,ibin)**2
  enddo
  enddo
  enddo
  call pencil_fft_backward
  write(15) r3 ! Wiener filtered delta_L

  if (xi_input(1,1)/=0) then
    print*, 'xi_input(1,1) =',xi_input(1,1)
    print*,'Detected nonzero input filter, use delta_R filter for delta_L.'
    xi_input(9,:)=1
  else
    print*,'Input filter is zero, use computed filter for delta_R.'
    xi_input=xi !
  endif

  cxyz=0
  do k=1,npen
  do j=1,ng
  do i=1,ng*nn/2+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ng+j
    ig=i
    kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
#ifdef linear_kbin
      ibin=nint(kr)
#else
      rbin=4.0/log(2.)*log(kr/0.95)
      ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
#endif
    ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
    cxyz(i,j,k)=cx2(i,j,k)*xi_input(8,ibin)**2
  enddo
  enddo
  enddo
  call pencil_fft_backward
  write(15) r3 ! Wiener filtered delta

  close(15)

endsubroutine

endmodule
