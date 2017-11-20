! read in checkpoint and compute the gravitational potential
! CIC potential in Lagrangian space
program potential
  use pencil_fft
  implicit none
  save

  logical,parameter :: correct_kernel=.true.
  logical,parameter :: write_potential=.true.

  integer(8),parameter :: npnode=nf**3 ! only true for this project
  real,parameter :: density_buffer=1.2
  integer(8),parameter :: npmax=npnode*density_buffer

  integer(8) i,j,k,l,i_dim,iq(3),nplocal,itx,ity,itz,nlast,ip,np,idx1(3),idx2(3)
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)
  real(4) rho_grid(0:ng+1,0:ng+1,0:ng+1)[*]
  real phi(0:nf+1,0:nf+1,0:nf+1)[*],phiq(nf,nf,nf),phi_temp
  real(4) mass_p,pos0(3),pos1(3),dx1(3),dx2(3)
  integer(izipx) xp(3,npmax)
  integer(8) pid(npmax)
  integer(8) ind,dx,dxy,kg,mg,jg,ig,ii,jj,kk
  integer(8) ileft,iright,nlen,g(3)
  real kr,kx,ky,kz,phi8,temp8[*]
  real dsp(3,ng,ng,ng)
  complex cx_temp(nyquest+1,nf,npen)

  call geometry
  if (head) then
    print*, 'potential on resolution:'
    print*, 'ng=',ng
    print*, 'ng*nn=',ng*nn
  endif

  if (head) then
    print*, 'checkpoint at:'
    open(16,file='../main/redshifts.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      print*, z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    print*,''
  endif
  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]

  call create_penfft_plan

  do cur_checkpoint=1,n_checkpoint
    ! read in fine grid cic density field
    open(10,file=output_name('den'),status='old',access='stream')
    !open(10,file=output_dir()//'delta_L'//output_suffix(),status='old',access='stream')
    read(10) r3
    close(10)

    ! solve Poisson equation
    call pencil_fft_forward
    cx_temp=cxyz

    if (head) print*, ''
    if (head) print*, 'Potential field'
    do k=1,npen
    do j=1,nf
    do i=1,nyquest+1
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*nf+j
      ig=i
      kz=mod(kg+nyquest-1,nf_global)-nyquest
      ky=mod(jg+nyquest-1,nf_global)-nyquest
      kx=ig-1
      kz=2*sin(pi*kz/nf_global)
      ky=2*sin(pi*ky/nf_global)
      kx=2*sin(pi*kx/nf_global)
      kr=kx**2+ky**2+kz**2
      kr=max(kr,1.0/nf_global**2) ! avoid kr being 0
      cxyz(i,j,k)=-4*pi/kr
    enddo
    enddo
    enddo
    if (head) cxyz(1,1,1)=0 ! DC frequency
    sync all

    if (correct_kernel) then
      if (head) print*, 'correct kernel'
      call pencil_fft_backward

      temp8=0
      if (image==1) temp8=temp8+r3(9,1,1)+r3(1,9,1)+r3(1,1,9)
      !if (image==1) print*, r3(9,1,1),r3(1,9,1),r3(1,1,9)
      sync all
      if (icx==nn .and. icy==1 .and. icz==1) temp8=temp8+r3(nf-7,1,1)
      !if (icx==nn .and. icy==1 .and. icz==1) print*, r3(nf-7,1,1)
      sync all
      if (icx==1 .and. icy==nn .and. icz==1) temp8=temp8+r3(1,nf-7,1)
      !if (icx==1 .and. icy==nn .and. icz==1) print*, r3(1,nf-7,1)
      sync all
      if (icx==1 .and. icy==1 .and. icz==nn) temp8=temp8+r3(1,1,nf-7)
      !if (icx==1 .and. icy==1 .and. icz==nn) print*, r3(1,1,nf-7)
      sync all
      phi8=0
      do i=1,nn**3
        phi8=phi8+temp8[i]
      enddo
      sync all
      phi8=phi8/6
      if (head) print*,'phi8 =',phi8
      sync all
      ! Construct Ewald potential kernel in real space
      if (head) print*, 'Construct Ewald potential kernel in real space'
      do k=1,nf
      do j=1,nf
      do i=1,nf
        kg=k+nf*(icz-1)
        jg=j+nf*(icy-1)
        ig=i+nf*(icx-1)
        kx=mod(kg+nyquest-1,nf_global)-nyquest
        ky=mod(jg+nyquest-1,nf_global)-nyquest
        kz=mod(ig+nyquest-1,nf_global)-nyquest
        kr=sqrt(kx**2+ky**2+kz**2)
        if (kr>8) then
          r3(i,j,k)=r3(i,j,k)-(phi8+1/8.)
        elseif (kr>0) then
          r3(i,j,k)=-1/kr
        elseif (kr==0) then
          r3(i,j,k)=-2.5
        endif
      enddo
      enddo
      enddo
      sync all
      !print*,'ewarld kernel',image,sum(r3*1d0)
      call pencil_fft_forward
    endif
    ! Complex multiply delta_L with potential kernel
    cxyz=real(cxyz)*cx_temp
    call pencil_fft_backward
    phi=0
    phi(1:nf,1:nf,1:nf)=r3 ! phi1
    phi(0,:,:)=phi(nf,:,:)[image1d(inx,icy,icz)]
    phi(nf+1,:,:)=phi(1,:,:)[image1d(ipx,icy,icz)]
    sync all
    phi(:,0,:)=phi(:,nf,:)[image1d(icx,iny,icz)]
    phi(:,nf+1,:)=phi(:,1,:)[image1d(icx,ipy,icz)]
    sync all
    phi(:,:,0)=phi(:,:,nf)[image1d(icx,icy,inz)]
    phi(:,:,nf+1)=phi(:,:,1)[image1d(icx,icy,ipz)]
    sync all

    ! write potential
    if (write_potential) then
      open(11,file=output_name('phi1'),status='replace',access='stream')
      write(11) phi(1:nf,1:nf,1:nf)
      close(11)
    endif
    sync all

    ! transform potential into Lagrangian space
    !if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    !open(12,file=output_name('zip2'),status='old',action='read',access='stream')
    !read(12) sim
    !if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
    !  print*, 'zip format incompatable'
    !  close(12)
    !  stop
    !endif
    !read(12) rhoc(1:nt,1:nt,1:nt,:,:,:) ! coarse grid density
    !close(12)
    !mass_p=sim%mass_p
    !if (head) print*, 'mass_p =',mass_p
    !nplocal=sim%nplocal
    !if (head) print*, 'nplocal =',nplocal
    !open(10,file=output_name('zip0'),status='old',action='read',access='stream')
    !read(10) xp(:,:nplocal) ! particle Eulerian positions
    !close(10)
    !open(10,file=output_name('zipid'),status='old',action='read',access='stream')
    !read(10) pid(:nplocal) ! particle Eulerian positions
    !close(10)
    !
    open(10,file=output_name('dsp'),status='old',action='read',access='stream')
    read(10) dsp(1,:,:,:)
    read(10) dsp(2,:,:,:)
    read(10) dsp(3,:,:,:)
    close(10)

    do k=1,nf
    do j=1,nf
    do i=1,nf
      pos1 = real((/i,j,k/)-1)+0.5 + dsp(:,i,j,k)
      pos1=modulo(pos1,real(nf))
      pos1=pos1-0.5
      idx1=floor(pos1)+1
      idx2=idx1+1
      dx1=idx1-pos1
      dx2=1-dx1

      phi_temp=0
      phi_temp=phi_temp+dx1(1)*dx1(2)*dx1(3)*phi(idx1(1),idx1(2),idx1(3))
      phi_temp=phi_temp+dx2(1)*dx1(2)*dx1(3)*phi(idx2(1),idx1(2),idx1(3))
      phi_temp=phi_temp+dx1(1)*dx2(2)*dx1(3)*phi(idx1(1),idx2(2),idx1(3))
      phi_temp=phi_temp+dx1(1)*dx1(2)*dx2(3)*phi(idx1(1),idx1(2),idx2(3))
      phi_temp=phi_temp+dx1(1)*dx2(2)*dx2(3)*phi(idx1(1),idx2(2),idx2(3))
      phi_temp=phi_temp+dx2(1)*dx1(2)*dx2(3)*phi(idx2(1),idx1(2),idx2(3))
      phi_temp=phi_temp+dx2(1)*dx2(2)*dx1(3)*phi(idx2(1),idx2(2),idx1(3))
      phi_temp=phi_temp+dx2(1)*dx2(2)*dx2(3)*phi(idx2(1),idx2(2),idx2(3))
      phiq(i,j,k)=phi_temp
    enddo
    enddo
    enddo

    open(10,file=output_name('phiq'),status='replace',access='stream')
      write(10) phiq
    close(10)

  enddo
  call destroy_penfft_plan
end
