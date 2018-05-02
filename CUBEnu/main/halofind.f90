module halo_output
  implicit none

  type halo_type
    real hpos(3)
    real mass_vir,mass_odc,radius_vir,radius_odc,v_disp
    real x_mean(3),v_mean(3),ang_mom(3),var_x(3)
  endtype
endmodule

subroutine halofind
  use omp_lib
  use variables
  use neutrinos
  use halo_output
  implicit none
  save

  type(halo_type) halo_info
  integer,parameter :: nc_halo_max=128
  integer,parameter :: ngrid_max=300
  integer,parameter :: nlist=5*(nc_halo_max+1)**3
  integer,parameter :: max_maxima=5*nc_halo_max**3
  integer,parameter :: max_halo_np=5*(nc_halo_max+1)**3/3
  integer, parameter :: search_ratio = 4
  integer, parameter :: refine_ratio = 5
  logical,parameter :: NGP = .true. ! NGP mass assignment by default
  ! set .true. to limit initial halo mass calculation to exiting after only
  ! completing an entire radial shell
  logical, parameter      :: complete_shell = .true.
  logical physical_halo
  real,parameter :: den_peak_cutoff=100 ! lower density peak threshold for determining which peaks should be inspected as halos
  integer,parameter :: min_halo_particles=100 ! lower mass cutoff for cataloging a halo, in particles
  real,parameter :: halo_odc=200 ! spherical overdensity cutoff for determining halo extent


  integer i0,j0,k0,l0,idx1(3),idx2(3),irtot
  integer(1) hpart_odc(np_image_max),hpart_vir(np_image_max)
  integer(8) ilist_odc(max_halo_np),ilist_vir(max_halo_np),i_vir,i_odc
  integer idist(3,nlist),isortdist(nlist),idist_tmp(3,nlist),isortpeak(max_maxima)
  integer isortpos_vir(max_halo_np),isortpos_odc(max_halo_np)
  integer iloc,np_vir,np_odc,np_search
  integer crbox(3,2),csbox(3,2),frbox(3,2),itile(3),ngrid(3)
  real d1,d2,r1,r2,w1,w2,odci,odcj,dx(3),dv(3)
  real finegrid(ngrid_max,ngrid_max,ngrid_max)
  real xv_vir(6,max_halo_np),xv_odc(6,max_halo_np),r_vir(max_halo_np),r_odc(max_halo_np)
  real dx1(3),dx2(3),pos1(3),rr,rdist(nlist),halo_vir,xflat,den_peak(max_maxima),ipeak(3,max_maxima),amtot,denmax
  real hpos(3),mass_proxy,dgrid,rrefine,rsearch,dr(3)
  real rhof(1-nfb:nft+nfb,1-nfb:nft+nfb,1-nfb:nft+nfb),halo_mesh_mass(max_maxima)

  if (head) print*, ''
  if (head) print*, 'halofind'

  ! initialize_halofinder
  ! Loop through a box of length 2*nc_halo_max
  ! if cell is within sphere of radius = box length / 2
  ! include distince in rdist at entry ii
  ! ordered bottom left to top right
  print*, 'initialize_halofinder'
  ii=0
  do i0=-nc_halo_max,nc_halo_max
  do j0=-nc_halo_max,nc_halo_max
  do k0=-nc_halo_max,nc_halo_max
    !rr=norm2((/real(i0),real(j0),real(k0)/))
    rr=norm2([real :: i0,j0,k0])
    if (rr>nc_halo_max) cycle
    ii=ii+1
    if (ii>nlist) then
      print*, 'ii exceeded ',nlist
      stop
    endif
    idist(:,ii)=[i0,j0,k0]
    rdist(ii)=rr
  enddo
  enddo
  enddo
  irtot=ii
  ! sorts the rdist array from lowest to highest radial position
  ! from center of sphere saves rdist array position values in idist
  isortdist(:irtot)=(/ (i0,i0=1,irtot) /)
  call indexedsort(irtot,rdist,isortdist)
  idist(:,:irtot)=idist(:,isortdist(:irtot))

  ! Find halo candidates based on local overdensities for each tile
  ! Determine which candidates are to be considered halos and write their properties to file.
  ! Determine Delta_vir using equation (6) of Bryan et al. (1997) for a flat universe
  print*, 'scale_factor =',1./(1+z_halofind(cur_halofind))
  xflat=(omega_m/a**3)/(omega_m/a**3+omega_l)-1.
  halo_vir=18.*pi**2+82.*xflat-39.*xflat**2 ! <178
  if (head) print*, '  scale_factor,xflat,halo_vir =',a,xflat,halo_vir
  sync all
  if (halo_vir>halo_odc) stop 'Warning: halo_vir > halo_odc'

  if (head) then
    print*, '  mass_p       =',sim%mass_p_cdm
    print*, '  nplocal      =',sim%nplocal
    print*, '  nplocal_nu   =',sim%nplocal_nu
    print*, '  np_image_max =',np_image_max
  endif

  ! open halo file
  open(11,file=output_name_halo('halo'),status='replace',access='stream')
  open(12,file=output_name_halo('fail'),status='replace',access='stream')
  ! Determine which candidates are to be considered halos and write their properties to file.
  nhalo=0
  n_search_fail=0 ! Ticker recording how many particle searches end up empty handed
  hpart_odc=0 ! Initialize so that no particles are yet part of a halo
  hpart_vir=0
  write(11) nhalo_tot,nhalo,halo_vir,halo_odc

  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    print*, '  on tile',int(itx,1),int(ity,1),int(itz,1)
    print*, '  fine mass assignment'
    !call find_halo_candidates(tile, n_candidate)
    !subroutine find_halo_candidates(tile, ic)
    rhof=0; isortpeak=0; den_peak=0; n_candidate=0; n_candidate_real=0;
    !$omp paralleldo default(shared) &
    !$omp& private(k0,j0,i0,np,nzero,l0,ip,pos1,idx1)
    do k0=2-ncb,nt+ncb-1
    do j0=2-ncb,nt+ncb-1
    do i0=2-ncb,nt+ncb-1
      np=rhoc(i0,j0,k0,itx,ity,itz)
      nzero=idx_b_r(j0,k0,itx,ity,itz)-sum(rhoc(i0:,j0,k0,itx,ity,itz))
      do l0=1,np
        ip=nzero+l0
        if (NGP) then
          pos1=ncell*([i0,j0,k0]-1) + ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          idx1=floor(pos1)+1
          rhof(idx1(1),idx1(2),idx1(3))=rhof(idx1(1),idx1(2),idx1(3))+sim%mass_p_cdm
        else
          pos1=ncell*([i0,j0,k0]-1) + ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution - 0.5
          idx1=floor(pos1)+1
          idx2=idx1+1
          dx1=idx1-pos1
          dx2=1-dx1
          rhof(idx1(1),idx1(2),idx1(3))=rhof(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx2(1),idx1(2),idx1(3))=rhof(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx1(1),idx2(2),idx1(3))=rhof(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx1(1),idx1(2),idx2(3))=rhof(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          rhof(idx1(1),idx2(2),idx2(3))=rhof(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
          rhof(idx2(1),idx1(2),idx2(3))=rhof(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          rhof(idx2(1),idx2(2),idx1(3))=rhof(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx2(1),idx2(2),idx2(3))=rhof(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
        endif
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo

    ! find maxima
    print*, '  find_halo_candidates'
    do k0=1-ncell,nft+ncell  ! 1 more coarse cell layer
    do j0=1-ncell,nft+ncell
    do i0=1-ncell,nft+ncell
      denmax=maxval(rhof(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1))
      if (denmax==rhof(i0,j0,k0).and.denmax>den_peak_cutoff) then
        !print*,maxloc(rhof(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1))
        if (n_candidate>max_maxima-2) stop 'too many maxima'
        ! Find the fine mesh mass of this peak
        amtot=0
        do ii=1,irtot ! keep adding mass until mean density is less than halo_odc
          ix=i0+idist(1,ii); if(ix<=1-nfb .or. ix>=nft+nfb-1) cycle
          iy=j0+idist(2,ii); if(iy<=1-nfb .or. iy>=nft+nfb-1) cycle
          iz=k0+idist(3,ii); if(iz<=1-nfb .or. iz>=nft+nfb-1) cycle
          ! skip the outerest layer, due to CIC mass assignment
          amtot=amtot+rhof(ix,iy,iz)
          if (complete_shell .and. rdist(ii+1)==rdist(ii)) cycle
          if (ii==19 .and. amtot/ii<halo_odc) exit ! upto sqrt(2)
        enddo
        ! Consider this a halo candidate if the fine mesh mass is large enough
        if (amtot>sim%mass_p_cdm*min_halo_particles/2.) then
          n_candidate=n_candidate+1
          n_candidate_real=n_candidate_real+merge(1,0,minval([i0,j0,k0])>=1 .and. maxval([i0,j0,k0])<=nft)
          ! NGP maximum
          !ipeak(:,n_candidate)=[i0,j0,k0]-0.5
          ! parabolic interpolation maximum
          ipeak(1,n_candidate)=(i0-0.5)+0.5*(rhof(i0+1,j0,k0)-rhof(i0-1,j0,k0))&
          /(-rhof(i0-1,j0,k0)+2*rhof(i0,j0,k0)-rhof(i0+1,j0,k0))
          ipeak(2,n_candidate)=(j0-0.5)+0.5*(rhof(i0,j0+1,k0)-rhof(i0,j0-1,k0))&
          /(-rhof(i0,j0-1,k0)+2*rhof(i0,j0,k0)-rhof(i0,j0+1,k0))
          ipeak(3,n_candidate)=(k0-0.5)+0.5*(rhof(i0,j0,k0+1)-rhof(i0,j0,k0-1))&
          /(-rhof(i0,j0,k0-1)+2*rhof(i0,j0,k0)-rhof(i0,j0,k0+1))

          den_peak(n_candidate)=denmax
          halo_mesh_mass(n_candidate)=amtot
        endif
      endif
    enddo
    enddo
    enddo
    print*, '  real/total candidates =',n_candidate_real,n_candidate

    ! sort density maxima
    isortpeak(:n_candidate) = [(i0,i0=1,n_candidate)]
    call indexedsort(n_candidate,den_peak(:),isortpeak(:))
    ipeak(:,:n_candidate)=ipeak(:,isortpeak(:n_candidate))
    halo_mesh_mass(:n_candidate)=halo_mesh_mass(isortpeak(:n_candidate))

    do iloc=n_candidate,1,-1
      ! determine searching regions
      mass_proxy=halo_mesh_mass(iloc)
      hpos=ipeak(:,iloc)

      ! find_halo_particles
      ! Search for particles by looking at local particle distribution
      !call find_halo_particles(halo_vir, mass_proxy, hpos(:), r_vir, i_vir, 1)
      !call find_halo_particles(halo_odc, mass_proxy, hpos(:), r_odc, i_odc, 2)

      ! Use the mass proxy to guess at the size of the halo
      ! This will be the radius within which the refined density peak will be found
      rrefine=(0.75/pi/halo_vir*mass_proxy)**(1./3.)
      ! This (larger) radius will be used to store particle positions determine which are part of the halo
      rsearch = search_ratio * rrefine
      !print*,rrefine,rsearch
      !print*, nft
      itile=[itx,ity,itz]
      crbox(:,1)=floor((hpos-rrefine)/ncell)+1
      crbox(:,2)=floor((hpos+rrefine)/ncell)+1
      csbox(:,1)=floor((hpos-rsearch)/ncell)+1
      csbox(:,2)=floor((hpos+rsearch)/ncell)+1
      frbox(:,1)=ncell*(crbox(:,1) - 1) ! Boundary of the refined region in fine mesh cell units, integer coordinates
      frbox(:,2)=ncell*crbox(:,2)
      ngrid=refine_ratio*(frbox(:,2)-frbox(:,1)) ! Number of extra refined cells in this region
      dgrid    = 1./refine_ratio ! and their spacing
      if (maxval(ngrid)>ngrid_max) stop 'ngrid>ngrid_max'

      ! find halo particles
      np_search=0
      np_vir=0; np_odc=0
      finegrid=0
      xv_vir=0; xv_odc=0
      !open(13,file='ilist_vir.dat')
      !print*,'search particle'
      do k0=csbox(3,1),csbox(3,2)
      do j0=csbox(2,1),csbox(2,2)
      do i0=csbox(1,1),csbox(1,2)
        !nlast=cum(i0-1,j0,k0,itile(1),itile(2),itile(3))
        np=rhoc(i0,j0,k0,itile(1),itile(2),itile(3))
        nzero=idx_b_r(j0,k0,itx,ity,itz)-sum(rhoc(i0:,j0,k0,itx,ity,itz))
        do l0=1,np
          !ip=nlast+l0
          !if(ip/=nzero+l0) then
          !   print*,'1',nzero,nlast
          !   stop
          !endif
          ip=nzero+l0
          pos1=ncell*([i0,j0,k0]-1)+ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          dr=pos1-hpos
          rr=norm2(dr)
          if (rr<rsearch) then
            if (hpart_vir(ip)==0) then
              np_vir=np_vir+1
              xv_vir(1:3,np_vir)=pos1
              xv_vir(4:6,np_vir)=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost)) &
                                 +vfield(:,i0,j0,k0,itx,ity,itz)
              ilist_vir(np_vir)=ip
              if (rr < rrefine) then
                idx1=int((pos1-frbox(:,1))/dgrid)+1
                finegrid(idx1(1),idx1(2),idx1(3))=finegrid(idx1(1),idx1(2),idx1(3))+1
              endif
            endif
            if (hpart_odc(ip)==0) then
              np_odc=np_odc+1
              xv_odc(1:3,np_odc)=pos1
              xv_odc(4:6,np_odc)=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost)) &
                                +vfield(:,i0,j0,k0,itx,ity,itz)
              ilist_odc(np_odc)=ip
            endif
          endif
        enddo
      enddo
      enddo
      enddo

      if(np_odc>max_halo_np) stop 'np_search>max_halo_np'
      hpos=frbox(:,1)+dgrid*(maxloc(finegrid)-0.5) ! Find refined mesh density maximum, wrt tile
      physical_halo=(minval(hpos)>0 .and. maxval(hpos)<nft)
      ! sort by radius -------------------------------------------------------
      r_vir(1:np_vir)=norm2(xv_vir(1:3,1:np_vir)-spread(hpos,2,np_vir),1)
      isortpos_vir(:np_vir)=[(i0,i0=1,np_vir)]
      call indexedsort(np_vir,r_vir,isortpos_vir)
      r_odc(1:np_odc)=norm2(xv_odc(1:3,1:np_odc)-spread(hpos,2,np_odc),1)
      isortpos_odc(:np_odc)=[(i0,i0=1,np_odc)]
      call indexedsort(np_odc,r_odc,isortpos_odc)
      ! calculate halo radius ------------------------------------------------
      i_vir=0; i_odc=0
      do ii=2,np_vir
        odcj=0.75/pi*ii*sim%mass_p_cdm/r_vir(ii)**3
        if (odcj<halo_vir) then
          odci=0.75/pi*(ii-1)*sim%mass_p_cdm/r_vir(ii-1)**3
          r2 = log10(r_vir(ii))
          r1 = log10(r_vir(ii-1))
          d2 = log10(odcj)
          d1 = log10(odci)
          w1 = log10(halo_vir) - d2
          w2 = d1 - log10(halo_vir)
          halo_info%radius_vir = 10**((w1 * r1 + w2 * r2) / (d1 - d2))
          i_vir=ii-1
          exit
        endif
      enddo
      do ii=2,np_odc
        odcj=0.75/pi*ii*sim%mass_p_cdm/r_odc(ii)**3
        if (odcj<halo_odc) then
          odci=0.75/pi*(ii-1)*sim%mass_p_cdm/r_odc(ii-1)**3
          r2 = log10(r_odc(ii))
          r1 = log10(r_odc(ii-1))
          d2 = log10(odcj)
          d1 = log10(odci)
          w1 = log10(halo_odc) - d2
          w2 = d1 - log10(halo_odc)
          halo_info%radius_odc = 10**((w1 * r1 + w2 * r2) / (d1 - d2))
          i_odc=ii-1
          exit
        endif
      enddo

      if (i_vir==0 .or. i_odc==0) then ! check number of particles
        n_search_fail=n_search_fail+1
        print*,'  search fail:'
        print*, ipeak(:,iloc)+(itile-1)*nft
        print*,'  r',rrefine,rsearch
        write(12) ipeak(:,iloc)+(itile-1)*nft
        !stop
      endif

      if (i_vir >= min_halo_particles) then
        hpart_vir(ilist_vir(1:i_vir))=1
        hpart_odc(ilist_odc(1:i_odc))=1
        if (physical_halo) then  ! if the maximum is in physical regions
          nhalo=nhalo+1
          halo_info%hpos=hpos+(itile-1)*nft
          halo_info%mass_vir=sim%mass_p_cdm*i_vir
          halo_info%mass_odc=sim%mass_p_cdm*i_odc
          halo_info%x_mean=sum(xv_vir(1:3,isortpos_vir(1:i_vir)),2)/i_vir
          halo_info%var_x=sum((xv_vir(1:3,isortpos_vir(1:i_vir))-spread(halo_info%x_mean,2,i_vir))**2,2)/(i_vir-1)
          halo_info%v_mean=sum(xv_vir(4:6,isortpos_vir(1:i_vir)),2)/i_vir
          halo_info%v_disp=sqrt(sum((xv_vir(4:6,isortpos_vir(1:i_vir))-spread(halo_info%v_mean,2,i_vir))**2)/(i_vir-1))
          halo_info%ang_mom=0
          do ii=1,i_vir
            dx=xv_vir(1:3,ii)-halo_info%hpos
            dv=xv_vir(4:6,ii)-halo_info%v_mean
            halo_info%ang_mom(1)=halo_info%ang_mom(1)+dx(2)*dv(3)-dx(3)*dv(2)
            halo_info%ang_mom(2)=halo_info%ang_mom(2)+dx(3)*dv(1)-dx(1)*dv(3)
            halo_info%ang_mom(3)=halo_info%ang_mom(3)+dx(1)*dv(2)-dx(2)*dv(1)
          enddo
          write(11) halo_info ! if the maximum is in physical regions
        endif ! physical_halo
      endif
    enddo ! iloc
    print*,'  found nhalo =',nhalo
  enddo
  enddo
  enddo ! end looping over tiles
  sync all
  if (head) then
    nhalo_tot=0
    do ii=1,nn**3
      nhalo_tot=nhalo_tot+nhalo[ii]
    enddo
  endif
  sync all
  nhalo_tot=nhalo_tot[1]
  sync all
  if (head) then
    print*, '  found',nhalo_tot,'halos'
  endif
  rewind(11)
  write(11) nhalo_tot,nhalo
  close(11)
  close(12)

  cur_halofind=cur_halofind+1
  halofind_step=.false.

  sync all

endsubroutine
