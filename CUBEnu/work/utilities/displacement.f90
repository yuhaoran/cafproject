!! add -DRECONSTRUCTION to compute delta_E from dsp
program displacement
  use parameters
#ifdef Emode
  use powerspectrum
#endif
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  integer(8),parameter :: npnode=nf**3 ! only true for this project
  real,parameter :: density_buffer=1.2
  integer(8),parameter :: npmax=npnode*density_buffer
  integer(8) i,j,k,l,i_dim,iq(3),pid8,itx,ity,itz,nlast,ip,np
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt),rhoc0(nt,nt,nt,nnt,nnt,nnt)
  integer(4) rho0(ng,ng,ng) !!!! for checking there is one and only one particle per fine grid
  real dsp(3,ng,ng,ng),pos0(3),pos1(3),dpos(3)

  integer(izipx),allocatable :: xp(:,:)
  integer(8),allocatable :: pid(:)

#ifdef Emode
  integer dim_1,dim_2,dim_3,ig,jg,kg
  real kr,kx(3),cube1(ng,ng,ng),cube0(ng,ng,ng),xi(10,nbin)[*]
  complex ekx(3),cdiv(ng*nn/2+1,ng,npen),pdim
#endif

#ifdef dspE
  integer(8) idx1(3),idx2(3)
  real(4) dx1(3),dx2(3)
  complex kpsi(ng*nn/2+1,ng,npen)
  real(4) rho_grid(0:ng+1,0:ng+1,0:ng+1)
#endif

  call geometry

  if (head) then
    print*, 'Displacement field analysis on resolution:'
    print*, 'ng=',ng
  endif
  sync all

  if (head) then
    print*, 'checkpoint at:'
    open(16,file='../main/z_checkpoint.txt',status='old')
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
  sync all

#ifdef Emode
  call create_penfft_plan
#endif

  do cur_checkpoint= n_checkpoint,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    !call print_header(sim); stop
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, 'zip format incompatable'
      close(11)
      stop
    endif
    if (head) then
      print*, 'nplocal =',sim%nplocal
    endif
    !cdm
    allocate(xp(3,sim%nplocal))
    allocate(pid(sim%nplocal))
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)

    open(14,file=output_name('id'),status='old',action='read',access='stream')
    print*, output_name('id')
    read(14) pid ! particle Lagrangian positions
    close(14)
    print*,'check PID range:',minval(pid),maxval(pid)

    rho0=0
    dsp=0
    nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      if (head) print*, 'CIC interpolation on tile',int(itx,1),int(ity,1),int(itz,1)
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pid8=pid(ip)-1
          iq(3)=pid8/int(nf_global,4)**2
          iq(2)=(pid8-iq(3)*int(nf_global,4)**2)/int(nf_global,4)
          iq(1)=modulo(pid8,int(nf_global,4))

          pos0=iq+0.5
          pos1=nt*((/itx,ity,itz/)-1) + (/i,j,k/)-1 + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          pos1=real(ng)*((/icx,icy,icz/)-1) + pos1*real(ng)/real(nc)

          dpos=pos1-pos0
          dpos=modulo(dpos+ng*nn/2,real(ng*nn))-ng*nn/2
          dpos=dpos*real(nf)/real(ng)

          dsp(:,iq(1)+1,iq(2)+1,iq(3)+1)=dpos
          rho0(iq(1)+1,iq(2)+1,iq(3)+1)=rho0(iq(1)+1,iq(2)+1,iq(3)+1)+1
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    deallocate(xp,pid)

    print*, 'check: min,max,sum of rho0 = '
    print*, minval(rho0),maxval(rho0),sum(rho0)
    if (minval(rho0)<1 .or. maxval(rho0)>1 .or. sum(rho0)/=sim%nplocal) then
      print*, 'Warning'
      !stop
    endif

    do i_dim=1,3
      print*, 'dsp: dimension',int(i_dim,1),'min,max values ='
      print*, minval(dsp(i_dim,:,:,:)), maxval(dsp(i_dim,:,:,:))
    enddo

    if (head) print*,'Write dsp into file'
    open(15,file=output_name('dsp'),status='replace',access='stream')
      write(15) dsp(1,1:ng,1:ng,1:ng)
      write(15) dsp(2,1:ng,1:ng,1:ng)
      write(15) dsp(3,1:ng,1:ng,1:ng)
    close(15)

#ifdef Emode
      if (head) print*,'Start computing delta_E'
      !cphi=0
      cdiv=0
      do i_dim=1,3
        if (head) print*,'Start working on dim',int(i_dim,1)
        r3=dsp(i_dim,1:ng,1:ng,1:ng)
        if (head) print*,'start forward tran'
        call pencil_fft_forward
        if (head) print*,'loop over k'
        ! cxyz is the fourier of dsp(i_dim,1:ng,1:ng,1:ng)
        do k=1,npen
        do j=1,ng
        do i=1,ng*nn/2+1
          kg=(nn*(icz-1)+icy-1)*npen+k
          jg=(icx-1)*ng+j
          ig=i
          kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
          kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
          ekx=exp(2*pi*(0,1)*kx/ng)
          dim_1=i_dim
          dim_2=mod(dim_1,3)+1
          dim_3=mod(dim_2,3)+1
          pdim=(ekx(dim_1)-1)*(ekx(dim_2)+1)*(ekx(dim_3)+1)/4
          !cphi(i,j,k)=cphi(i,j,k)+cxyz(i,j,k)*pdim/(-4*sum(sin(pi*kx/ng)**2)+0.000001)
          cdiv(i,j,k)=cdiv(i,j,k)+cxyz(i,j,k)*pdim
        enddo
        enddo
        enddo
      enddo ! i_dim

      if (head) then
        !cphi(1,1,1)=0
        cdiv(1,1,1)=0
      endif
      sync all

      !! reconstructed delta
      cxyz=cdiv
      if (head) print*,'start backward tran'
      call pencil_fft_backward
      cube1=-r3
      if (head) print*,'Write delta_E into file'
      open(15,file=output_name('delta_E'),status='replace',access='stream')
      write(15) cube1
      close(15)
      sync all

      if (head) print*,'Read delta_L from file'
      !! linear delta
      open(15,file=output_dir()//'delta_L'//output_suffix(),status='old',access='stream')
      read(15) cube0
      close(15)
      call cross_power(xi,cube0,cube1)
      open(15,file=output_name('power_LE'),status='replace',access='stream')
      write(15) xi
      close(15)
#endif

#ifdef dspE
      kpsi=0
      do i_dim=1,3
        r3=dsp(i_dim,1:ng,1:ng,1:ng)
        call pencil_fft_forward
        do k=1,npen
        do j=1,ng
        do i=1,ng*nn/2+1
          kg=(nn*(icz-1)+icy-1)*npen+k
          jg=(icx-1)*ng+j
          ig=i
          kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
          kpsi(i,j,k)=kpsi(i,j,k)+kx(i_dim)*cxyz(i,j,k)
        enddo
        enddo
        enddo
      enddo
      print*, 'sum of kpsi',sum(kpsi)

      open(16,file=output_name('dsp_E'),status='replace',access='stream')
      print*,'compute and write into',output_name('dsp_E')
      do i_dim=1,3
        do k=1,npen
        do j=1,ng
        do i=1,ng*nn/2+1
          kg=(nn*(icz-1)+icy-1)*npen+k
          jg=(icx-1)*ng+j
          ig=i
          kx=mod((/ig,jg,kg/)+ng/2-1,ng)-ng/2
          cxyz(i,j,k)=kpsi(i,j,k)*kx(i_dim)/sum(kx**2)
        enddo
        enddo
        enddo
        if (head) then
          !cphi(1,1,1)=0
          cxyz(1,1,1)=0
        endif
        print*, 'sum of cxyz',minval(abs(cxyz)),maxval(abs(cxyz))
        call pencil_fft_backward
        print*, 'r3',minval(r3),maxval(r3)
        write(16) r3
      enddo
      close(16)

      ! CIC to get the Eulerian density field
      open(16,file=output_name('dsp_E'),access='stream')
      read(16) dsp(1,:,:,:)
      read(16) dsp(2,:,:,:)
      read(16) dsp(3,:,:,:)
      close(16)
      if (head) print*, 'CIC interpolation by dsp_E'
      rho_grid=0
      do k=1,ng
      do j=1,ng
      do i=1,ng

          pos1=[i,j,k]-0.5+dsp(:,i,j,k)
          pos1=modulo(pos1-1,real(ng))+1
          pos1=pos1-0.5

          !pos1=250;


          idx1=floor(pos1)+1
          idx2=idx1+1
          dx1=idx1-pos1
          dx2=1-dx1

          !print*, idx1;stop

          rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
          rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
          rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
          rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
          rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
          rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
          rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
          rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)

      enddo
      enddo
      enddo
      rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)
      rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)
      rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)
      rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)
      rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)
      rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)
      rho_grid=rho_grid-1

      print*, minval(rho_grid(1:ng,1:ng,1:ng)),maxval(rho_grid(1:ng,1:ng,1:ng)),sum(rho_grid(1:ng,1:ng,1:ng)*1d0)
      open(16,file=output_name('delta_cE'),status='replace',access='stream')
      print*,'compute and write into',output_name('delta_cE')
      write(16) rho_grid(1:ng,1:ng,1:ng)
      close(16)

      ! check
      open(16,file=output_name('dsp'),access='stream')
      read(16) dsp(1,:,:,:)
      read(16) dsp(2,:,:,:)
      read(16) dsp(3,:,:,:)
      close(16)
      if (head) print*, 'CIC interpolation by dsp_E'
      rho_grid=0
      do k=1,ng
      do j=1,ng
      do i=1,ng

        pos1=[i,j,k]-0.5+dsp(:,i,j,k)
        pos1=modulo(pos1-1,real(ng))+1
        pos1=pos1-0.5

        idx1=floor(pos1)+1
        idx2=idx1+1
        dx1=idx1-pos1
        dx2=1-dx1

        rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
        rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
        rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
        rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
        rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
        rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
        rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
        rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)

      enddo
      enddo
      enddo
      rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)
      rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)
      rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)
      rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)
      rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)
      rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)
      rho_grid=rho_grid-1
      print*, minval(rho_grid(1:ng,1:ng,1:ng)),maxval(rho_grid(1:ng,1:ng,1:ng)),sum(rho_grid(1:ng,1:ng,1:ng)*1d0)
      open(16,file=output_name('delta_check'),status='replace',access='stream')
      print*,'compute and write into',output_name('delta_check')
      write(16) rho_grid(1:ng,1:ng,1:ng)
      close(16)

#endif

    !open(15,file='../output/universe1/image1/0.000dsp_1.bin',access='stream')
    !read(15) dsp0(1,1:ng,1:ng,1:ng)
    !read(15) dsp0(2,1:ng,1:ng,1:ng)
    !read(15) dsp0(3,1:ng,1:ng,1:ng)
    !close(15)
    !print*, 'max dsp-dsp0', maxval(abs(dsp-dsp0))

  enddo
#ifdef Emode
  call destroy_penfft_plan
#endif

  print*,'displacement done'
end
