program ang_mom_corr
  use parameters
  use halo_output, only: type_halo_catalog, type_halo_info
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim

  integer,parameter :: nc_halo_max=128
  integer,parameter :: newbox_max=300
  integer,parameter :: nlist=5*(nc_halo_max+1)**3
  integer,parameter :: n_peak_max=5*nc_halo_max**3
  integer,parameter :: max_halo_np=5*(nc_halo_max+1)**3/3
  integer(8) i,j,k,l
  integer(4) ihalo,nhalo,ip,np,ncheck
  character(20) str_z,str_i
  integer pid(max_halo_np),ipos(3)
  real qpos(3,max_halo_np),qpos_mean(3),vq(3),dx(3),vf,scale_factor,dx_mean(3)
  real spin_q(3),spin_x(3),spin_t(3),spin_u(3)
  real inertia(3,3),tide(3,3),tidu(3,3),torque(3,3)
  real force_cdm(3),force_neu(3),torque_c,torque_u
  real theta_qx(3000),theta_qt(3000),theta_tx(3000),theta_cu(3000),theta_ff(3000)
  real phi(0:nf+1,0:nf+1,0:nf+1)[*]
  real phu(0:nf+1,0:nf+1,0:nf+1)[*]

  type(type_halo_info) halo_info
  type(type_halo_catalog) halo_catalog
  real,allocatable :: corr_info(:,:)

  call geometry
  if (head) then
    print*, 'angular momemtum correlation on resolution:'
    print*, 'ng=',ng
    print*, 'ng*nn=',ng*nn
  endif
  sync all

  if (head) then
    print*, 'use halo catalog at:'
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

  cur_checkpoint=1
  print*,output_name('phi1')
  open(21,file=output_name('phi1'),status='old',access='stream')
  !read(21) phi(1:nf,1:nf,1:nf)
  !open(21,file='../../output/universe2/image1/5.000_phi_z5_1.bin',status='old',access='stream')
  !open(21,file='../../output/universe2/image1/0.000_phiE_1.bin',status='old',access='stream')
  read(21) phi(1:nf,1:nf,1:nf)  !! attention 1 layer buffer also loaded
  close(21)
  print*, output_name('phi1_nu')
  open(21,file=output_name('phi1_nu'),status='old',access='stream')
  !open(21,file='../../output/universe2/image1/5.000_tf5.000_phi1_nu_1.bin',status='old',access='stream')
  read(21) phu(1:nf,1:nf,1:nf)
  close(21)
  vf=vfactor(1/(1+z_checkpoint(cur_checkpoint)))

  if (head) print*, '  buffer phi'
  phi(0,:,:)=phi(nf,:,:)[image1d(inx,icy,icz)]
  phi(nf+1,:,:)=phi(1,:,:)[image1d(ipx,icy,icz)]
  sync all
  phi(:,0,:)=phi(:,nf,:)[image1d(icx,iny,icz)]
  phi(:,nf+1,:)=phi(:,1,:)[image1d(icx,ipy,icz)]
  sync all
  phi(:,:,0)=phi(:,:,nf)[image1d(icx,icy,inz)]
  phi(:,:,nf+1)=phi(:,:,1)[image1d(icx,icy,ipz)]
  sync all
  phu(0,:,:)=phu(nf,:,:)[image1d(inx,icy,icz)]
  phu(nf+1,:,:)=phu(1,:,:)[image1d(ipx,icy,icz)]
  sync all
  phu(:,0,:)=phu(:,nf,:)[image1d(icx,iny,icz)]
  phu(:,nf+1,:)=phu(:,1,:)[image1d(icx,ipy,icz)]
  sync all
  phu(:,:,0)=phu(:,:,nf)[image1d(icx,icy,inz)]
  phu(:,:,nf+1)=phu(:,:,1)[image1d(icx,icy,ipz)]
  sync all

  do cur_checkpoint= n_checkpoint,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    scale_factor=1/(1+z_checkpoint(cur_checkpoint))
    print*,output_name('halo')
    open(11,file=output_name('halo'),status='old',access='stream')
    open(12,file=output_name('halo_pid'),status='old',access='stream')
    open(13,file=output_name('halo_init_spin'),status='replace',access='stream')
    
    read(11) halo_catalog

    print*,' N_halos_global =',halo_catalog%nhalo_tot
    print*,' N_halos_local  =',halo_catalog%nhalo
    print*,' Overdensity    =',halo_catalog%den_odc
    nhalo=halo_catalog%nhalo
    allocate(corr_info(2,nhalo))

    do ihalo=1,nhalo
      read(11) halo_info
      read(12) np
      read(12) pid(:np)
      read(12) ncheck
      if (ncheck/=0) stop "halo_pid file error"
      if (minval(pid(:np))<1) stop "pids are not all positive numbers"
      pid(:np)=pid(:np)-1
      dx_mean=0
      do ip=1,np
        qpos(3,ip)=pid(ip)/nf_global**2
        qpos(2,ip)=(pid(ip)-(pid(ip)/nf_global**2)*nf_global**2)/nf_global
        qpos(1,ip)=modulo(pid(ip),nf_global)
        qpos(:,ip)=qpos(:,ip)+0.5
        dx=qpos(:,ip)-halo_info%x_mean
        dx=modulo(dx+nf_global/2,real(nf_global))-nf_global/2
        dx_mean=dx_mean+dx
      enddo
      dx_mean=dx_mean/np
      qpos_mean=halo_info%x_mean+dx_mean
      qpos_mean=modulo(qpos_mean,real(nf_global))
      if(ihalo==1) then
        print*,'check halo position'
        print*,halo_info%x_mean
        print*,dx_mean
        print*,halo_info%x_mean+dx_mean
        print*,qpos_mean
      endif
      spin_q=0; spin_u=0
      force_cdm=0; force_neu=0
      torque_c=0; torque_u=0
      dx_mean=0
      inertia=0
      tide=0
      tidu=0
      do ip=1,np
        ipos(3)=pid(ip)/nf_global**2
        ipos(2)=(pid(ip)-ipos(3)*nf_global**2)/nf_global
        ipos(1)=modulo(pid(ip),nf_global)
        ipos=ipos+1

        dx=ipos-0.5-qpos_mean
        dx=modulo(dx+nf_global/2,real(nf_global))-nf_global/2
        dx_mean=dx_mean+dx
        !CDM ===================================================================
        vq(1)=phi(ipos(1)+1,ipos(2),ipos(3))-phi(ipos(1)-1,ipos(2),ipos(3))
        vq(2)=phi(ipos(1),ipos(2)+1,ipos(3))-phi(ipos(1),ipos(2)-1,ipos(3))
        vq(3)=phi(ipos(1),ipos(2),ipos(3)+1)-phi(ipos(1),ipos(2),ipos(3)-1)
        vq=-vq/(8*pi)/Dgrow(1/(1+z_checkpoint(1)))*vf
        spin_q(1)=spin_q(1)+dx(2)*vq(3)-dx(3)*vq(2)
        spin_q(2)=spin_q(2)+dx(3)*vq(1)-dx(1)*vq(3)
        spin_q(3)=spin_q(3)+dx(1)*vq(2)-dx(2)*vq(1)
        force_cdm=force_cdm+vq
        !neutrinos
        vq(1)=phu(ipos(1)+1,ipos(2),ipos(3))-phu(ipos(1)-1,ipos(2),ipos(3))
        vq(2)=phu(ipos(1),ipos(2)+1,ipos(3))-phu(ipos(1),ipos(2)-1,ipos(3))
        vq(3)=phu(ipos(1),ipos(2),ipos(3)+1)-phu(ipos(1),ipos(2),ipos(3)-1)
        vq=-vq/(8*pi)/Dgrow(1/(1+z_checkpoint(1)))*vf
        spin_u(1)=spin_u(1)+dx(2)*vq(3)-dx(3)*vq(2)
        spin_u(2)=spin_u(2)+dx(3)*vq(1)-dx(1)*vq(3)
        spin_u(3)=spin_u(3)+dx(1)*vq(2)-dx(2)*vq(1)
        force_neu=force_neu+vq

        ! initial inertia
        inertia(1,1)=inertia(1,1)+dx(1)**2
        inertia(2,2)=inertia(2,2)+dx(2)**2
        inertia(3,3)=inertia(3,3)+dx(3)**2
        inertia(1,2)=inertia(1,2)+dx(1)*dx(2)
        inertia(2,3)=inertia(2,3)+dx(2)*dx(3)
        inertia(3,1)=inertia(3,1)+dx(3)*dx(1)

        ! initial tidal tensor
        tide(1,1)=tide(1,1)+phi(ipos(1)+1,ipos(2),ipos(3))-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1)-1,ipos(2),ipos(3))
        tide(2,2)=tide(2,2)+phi(ipos(1),ipos(2)+1,ipos(3))-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1),ipos(2)-1,ipos(3))
        tide(3,3)=tide(3,3)+phi(ipos(1),ipos(2),ipos(3)+1)-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1),ipos(2),ipos(3)-1)
        tide(1,2)=tide(1,2)+(phi(ipos(1)+1,ipos(2)+1,ipos(3))+phi(ipos(1)-1,ipos(2)-1,ipos(3))&
                            -phi(ipos(1)+1,ipos(2)-1,ipos(3))-phi(ipos(1)-1,ipos(2)+1,ipos(3)))/4
        tide(2,3)=tide(2,3)+(phi(ipos(1),ipos(2)+1,ipos(3)+1)+phi(ipos(1),ipos(2)-1,ipos(3)-1)&
                            -phi(ipos(1),ipos(2)+1,ipos(3)-1)-phi(ipos(1),ipos(2)-1,ipos(3)+1))/4
        tide(3,1)=tide(3,1)+(phi(ipos(1)+1,ipos(2),ipos(3)+1)+phi(ipos(1)-1,ipos(2),ipos(3)-1)&
                            -phi(ipos(1)+1,ipos(2),ipos(3)-1)-phi(ipos(1)-1,ipos(2),ipos(3)+1))/4

        tidu(1,1)=tidu(1,1)+phu(ipos(1)+1,ipos(2),ipos(3))-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1)-1,ipos(2),ipos(3))
        tidu(2,2)=tidu(2,2)+phu(ipos(1),ipos(2)+1,ipos(3))-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1),ipos(2)-1,ipos(3))
        tidu(3,3)=tidu(3,3)+phu(ipos(1),ipos(2),ipos(3)+1)-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1),ipos(2),ipos(3)-1)
        tidu(1,2)=tidu(1,2)+(phu(ipos(1)+1,ipos(2)+1,ipos(3))+phu(ipos(1)-1,ipos(2)-1,ipos(3))&
                            -phu(ipos(1)+1,ipos(2)-1,ipos(3))-phu(ipos(1)-1,ipos(2)+1,ipos(3)))/4
        tidu(2,3)=tidu(2,3)+(phu(ipos(1),ipos(2)+1,ipos(3)+1)+phu(ipos(1),ipos(2)-1,ipos(3)-1)&
                            -phu(ipos(1),ipos(2)+1,ipos(3)-1)-phu(ipos(1),ipos(2)-1,ipos(3)+1))/4
        tidu(3,1)=tidu(3,1)+(phu(ipos(1)+1,ipos(2),ipos(3)+1)+phu(ipos(1)-1,ipos(2),ipos(3)-1)&
                            -phu(ipos(1)+1,ipos(2),ipos(3)-1)-phu(ipos(1)-1,ipos(2),ipos(3)+1))/4


        spin_x=halo_info%ang_mom
      enddo
      dx_mean=dx_mean/np ! check mean
      if(maxval(abs(dx_mean))>5) then
        stop "dx incorrect"
      endif

      inertia(2,1)=inertia(1,2)
      inertia(3,2)=inertia(2,3)
      inertia(1,3)=inertia(3,1)
      tide(2,1)=tide(1,2)
      tide(3,2)=tide(2,3)
      tide(1,3)=tide(3,1)
      tidu(2,1)=tidu(1,2)
      tidu(3,2)=tidu(2,3)
      tidu(1,3)=tidu(3,1)
      tide=tide/np
      tidu=tidu/np
      torque=matmul(inertia,tide)
      !torque=matmul(tidu,tide) ! reconstructed

      spin_t(1)=-torque(2,3)+torque(3,2)
      spin_t(2)=-torque(3,1)+torque(1,3)
      spin_t(3)=-torque(1,2)+torque(2,1)

      write(13) spin_q,spin_t,spin_x
      theta_qx(ihalo)=sum(spin_q*spin_x)/sqrt(sum(spin_q**2))/sqrt(sum(spin_x**2))
      theta_qt(ihalo)=sum(spin_q*spin_t)/sqrt(sum(spin_q**2))/sqrt(sum(spin_t**2))
      theta_tx(ihalo)=sum(spin_t*spin_x)/sqrt(sum(spin_t**2))/sqrt(sum(spin_x**2))
      theta_cu(ihalo)=sum(spin_q*spin_u)/sqrt(sum(spin_q**2))/sqrt(sum(spin_u**2))
      torque_c=torque_c+sqrt(sum(spin_q**2))
      torque_u=torque_u+sqrt(sum(spin_u**2))
      theta_ff(ihalo)=sum(force_cdm*force_neu)/sqrt(sum(force_cdm**2))/sqrt(sum(force_neu**2))

      if (ihalo==1) then
        print*,'halo',ihalo
        print*,'inertia'
        print*,inertia
        print*,'tide'
        print*,tide
        print*,'spin'
        print*,spin_t
      endif

    enddo

    print*,'mean qx correlation =',sum(theta_qx(:nhalo))/nhalo
    print*,'mean qt correlation =',sum(theta_qt(:nhalo))/nhalo
    print*,'mean tx correlation =',sum(theta_tx(:nhalo))/nhalo
    print*,'mean cu correlation =',sum(theta_cu(:nhalo))/nhalo
    print*,'mean ff correlation =',sum(theta_ff(:nhalo))/nhalo
    print*,'torque_c/torque_u=',torque_u/torque_c

    close(11);close(12);close(13)
    deallocate(corr_info)
  enddo
  sync all

  contains
  function vfactor(a)
    implicit none
    real, parameter :: np = -(1./4.)+(5./4.)*sqrt(1-24.*omega_nu/omega_m/25.) !~1-3f/5
    real :: a
    real :: H,km,lm
    real :: vfactor
    lm=omega_l/omega_m
    km=(1-omega_m-omega_l)/omega_m
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H
    vfactor=vfactor*np!(1.-3.*(omega_nu/omega_m)/5.)
  endfunction vfactor
  function Dgrow(a)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real a
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=a*ga/g
  end function Dgrow
end
