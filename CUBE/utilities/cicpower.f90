program cicpower
  use parameters
  use pencil_fft
  use powerspectrum
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  integer(8),parameter :: npnode=nf**3 ! only true for this project
  real,parameter :: density_buffer=1.2
  integer(8),parameter :: npmax=npnode*density_buffer
  integer(8) i,j,k,l,i_dim,iq(3),nplocal,itx,ity,itz,nlast,ip,np,idx1(3),idx2(3)
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)
  real(4) rho_grid(0:ng+1,0:ng+1,0:ng+1)[*],rho0(ng,ng,ng),rho1(ng,ng,ng)
  real(4) mass_p,pos1(3),dx1(3),dx2(3)
  integer(izipx) x(3,npmax)

  real xi(10,nbin)[*]

  call geometry

  if (head) then
    print*, 'cicpower on resolution:'
    print*, 'ng=',ng
  endif
  sync all

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
  sync all

  call create_penfft_plan
  
  do cur_checkpoint= 1,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    open(12,file=output_name('zip2'),status='old',action='read',access='stream')
    read(12) sim
    ! check zip format and read rhoc
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, 'zip format incompatable'
      close(12)
      stop
    endif
    read(12) rhoc ! coarse grid density
    close(12)
    mass_p=sim%mass_p
    if (head) print*, 'mass_p =',mass_p

    nplocal=sim%nplocal
    if (head) print*, 'nplocal =',nplocal

    open(10,file=output_name('zip0'),status='old',action='read',access='stream')
    read(10) x(:,:nplocal) ! particle Eulerian positions
    close(10)

    rho0=0
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
          pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(x(:,ip)+ishift,izipx)+rshift)*x_resolution
          pos1=pos1*real(ng)/real(nc) - 0.5

          idx1=floor(pos1)+1
          idx2=idx1+1
          dx1=idx1-pos1
          dx2=1-dx1

          rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
          rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
          rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
          rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
          rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
          rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
          rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
          rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all
    print*, 'sum of rho_grid (with buffer) = '
    print*, sum(rho_grid*1d0)

    if (head) print*, 'Start sync from buffer regions'
    sync all
    rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)[image1d(inx,icy,icz)]
    rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)[image1d(ipx,icy,icz)]
    sync all
    rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)[image1d(icx,iny,icz)]
    rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)[image1d(icx,ipy,icz)]
    sync all
    rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)[image1d(icx,icy,inz)]
    rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)[image1d(icx,icy,ipz)]
    sync all

    rho1=rho_grid(1:ng,1:ng,1:ng)
    print*, 'check: min,max,sum of rho_grid = '
    print*, minval(rho1),maxval(rho1),sum(rho1*1d0)

    ! convert to density contrast
    rho1=rho1/(sum(rho1*1d0)/ng/ng/ng)-1

    if (head) print*,'Write rho_grid into file'
    open(15,file=output_name('delta_N'),status='replace',access='stream')
    write(15) rho1
    close(15)

    ! cross correlate correct density field
    open(15,file='../output/universe1/image1/0.000delta_N_1.bin',access='stream')
    read(15) rho0
    close(15)
    print*, 'mean error', sum(abs(rho1-rho0)*1d0)/ng/ng/ng

    call cross_power(xi,rho1,rho0)
    open(15,file=output_name('cicpower'),status='replace',access='stream')
    write(15) xi
    close(15)

  enddo
  call destroy_penfft_plan
  print*,'cicpower done'
end
