program displacement
  use parameters
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  integer(8),parameter :: npnode=nf**3 ! only true for this project
  real,parameter :: density_buffer=1.2
  integer(8),parameter :: npmax=npnode*density_buffer
  integer(8) i,j,k,l,i_dim,iq(3),pid8,nplocal,itx,ity,itz,nlast,ip,np
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt),rhoc0(nt,nt,nt,nnt,nnt,nnt)
  integer(4) rho0(ng,ng,ng) !!!! for checking there is one and only one particle per fine grid
  real dsp(3,ng,ng,ng),dsp0(3,ng,ng,ng),mass_p,pos0(3),pos1(3),dpos(3)

  integer(izipx) xp(3,npmax)
  integer(8)   pid(npmax)

  call geometry

  if (head) then
    print*, 'Displacement field analysis on resolution:'
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
    read(10) xp(:,:nplocal) ! particle Eulerian positions
    close(10)

    open(14,file=output_name('zipid'),status='old',action='read',access='stream')
    read(14) pid(:nplocal) ! particle Lagrangian positions
    close(14)
    print*,'check PID range:'
    print*,minval(pid(:nplocal)),maxval(pid(:nplocal)); !stop
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

    print*, 'check: min,max,sum of rho0 = '
    print*, minval(rho0),maxval(rho0),sum(rho0)
    if (minval(rho0)<1 .or. maxval(rho0)>1 .or. sum(rho0)/=nplocal) then
      print*, 'Warning'
      stop
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

    open(15,file='../output/universe1/image1/0.000dsp_1.bin',access='stream')
    read(15) dsp0(1,1:ng,1:ng,1:ng)
    read(15) dsp0(2,1:ng,1:ng,1:ng)
    read(15) dsp0(3,1:ng,1:ng,1:ng)
    close(15)
    print*, 'max dsp-dsp0', maxval(abs(dsp-dsp0))
    
  enddo
  print*,'displacement done'
end
