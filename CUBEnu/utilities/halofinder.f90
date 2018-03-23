module halo_output
  implicit none

  type halo_type
    real hpos(3)
    real mass_vir,mass_odc,radius_vir,radius_odc,v_disp
    real x_mean(3),v_mean(3),ang_mom(3),var_x(3)
  endtype
endmodule

program halofinder
  use omp_lib
  use parameters
  use halo_output
  !use buffer_particle_subroutines, only :: buffer_vp
  implicit none

  type(halo_type) halo_info
  integer,parameter :: nc_halo_max=128
  integer,parameter :: ngrid_max=300
  integer,parameter :: nlist=5*(nc_halo_max+1)**3
  integer,parameter :: max_maxima=5*nc_halo_max**3
  integer,parameter :: max_halo_np=5*(nc_halo_max+1)**3/3
  integer, parameter :: search_ratio = 4
  integer, parameter :: refine_ratio = 5
  real finegrid(ngrid_max,ngrid_max,ngrid_max)
  real xv_vir(6,max_halo_np),xv_odc(6,max_halo_np),r_vir(max_halo_np),r_odc(max_halo_np)
  logical(1) select_vir(max_halo_np),select_odc(max_halo_np)
  integer(8) ilist_odc(max_halo_np),ilist_vir(max_halo_np)

  integer(8),parameter :: np_image=(nc*np_nc)**3*merge(2,1,body_centered_cubic) ! average number of particles per image
  integer(8),parameter :: np_image_max=np_image*(nte*1./nt)**3*image_buffer
  integer(8),parameter :: np_tile_max=np_image/nnt**3*(nte*1./nt)**3*tile_buffer
  integer(izipx) xp(3,np_image_max)[*]
  integer(izipv) vp(3,np_image_max)[*]
  integer(1) hpart_odc(np_image_max),hpart_vir(np_image_max)
  integer(4) rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  real(4) vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt) ! cannot have >7 dims
  integer(8) cum(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]

  real rhof(1-nfb:nft+nfb,1-nfb:nft+nfb,1-nfb:nft+nfb)

  integer idist(3,nlist),isortdist(nlist),idist_tmp(3,nlist),isortpeak(max_maxima)
  integer isortpos_vir(max_halo_np),isortpos_odc(max_halo_np)
  real ipeak(3,max_maxima),den_peak(max_maxima),halo_mesh_mass(max_maxima)
  real rr,rdist(nlist),mass_p,pos1(3),dx1(3),dx2(3),amtot,denmax,scale_factor,xflat,halo_vir,radius_vir,radius_odc
  real hpos(3),mass_proxy,mass_odc,mass_vir,rrefine,rsearch,dgrid,dr(3)
  real :: odci, odcj, r1, r2, d1, d2, w1, w2, I_ij(6), E0, E_tmp
  integer i0,j0,k0,l0,ii,irtot,itx,ity,itz,idx1(3),idx2(3),ix,iy,iz,nhalo[*],nhalo_tot[*],iloc,itile(3)
  integer crbox(3,2),csbox(3,2),frbox(3,2),ngrid(3),np_search,imass_vir,imass_odc,np_vir,np_odc
  integer n_candidate[*],n_candidate_real[*],n_candidate_global,n_candidate_global_real
  integer(8) nlast,ip,np,nplocal,nplocal_nu,n_search_fail[*],i_vir,i_odc
  character(20) str_z,str_i

  real(4) :: v_disp,sigma_vi
  real(4), dimension(3) :: x_mean, x2_mean, var_x, v_mean, v2_mean, offset, dx,dv, l_CM, ang_mom
  real(4), dimension(3) :: r_wrt_halo, v_wrt_halo, v2_wrt_halo

  logical,parameter :: NGP = .true. ! NGP mass assignment by default
  ! set .true. to limit initial halo mass calculation to exiting after only
  ! completing an entire radial shell
  logical, parameter      :: complete_shell = .true.
  logical physical_halo
  real,parameter :: den_peak_cutoff=100 ! lower density peak threshold for determining which peaks should be inspected as halos
  integer,parameter :: min_halo_particles=100 ! lower mass cutoff for cataloging a halo, in particles
  real,parameter :: halo_odc=200 ! spherical overdensity cutoff for determining halo extent
  call geometry

  if (head) then
    print*, 'halofinder'
    print*, 'checkpoint at:'
    open(16,file='../main/redshifts.txt',status='old')
    do i0=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i0)
      print*, z_checkpoint(i0)
    enddo
    71 n_checkpoint=i0-1
    close(16)
    print*,''
  endif

  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  sync all

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
  !idist_tmp(:,:irtot)=idist(:,isortdist(:irtot))
  !idist(:,:irtot)=idist_tmp(:,:irtot)

  ! Find halo candidates based on local overdensities for each tile
  n_checkpoint=1
  do cur_checkpoint= 1,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    ! Determine which candidates are to be considered halos and write their properties to file.
    ! Determine Delta_vir using equation (6) of Bryan et al. (1997) for a flat universe
    scale_factor=1./(1+z_checkpoint(cur_checkpoint))
    xflat=(omega_m/scale_factor**3)/(omega_m/scale_factor**3+omega_l)-1.
    halo_vir=18.*pi**2+82.*xflat-39.*xflat**2 ! <178
    if (head) print*, '  scale_factor,xflat,halo_vir =',scale_factor,xflat,halo_vir
    sync all
    if (halo_vir>halo_odc) stop 'Warning: halo_vir > halo_odc'

    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    !call print_header(sim); stop
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, '  zip format incompatable'
      close(11)
      stop
    endif
    nplocal=sim%nplocal
    nplocal_nu=sim%nplocal_nu
    mass_p=sim%mass_p_cdm
    !mass_p=1.
    sigma_vi=sim%sigma_vi
    if (head) then
      print*, '  mass_p       =',mass_p
      print*, '  nplocal      =',nplocal
      print*, '  nplocal_nu   =',nplocal_nu
      print*, '  np_image_max =',np_image_max
    endif

    ! read zip checkpoints
    xp=0;vp=0;rhoc=0;vfield=0;
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp(:,:nplocal)
    close(11)
    open(11,file=output_name('vp'),status='old',action='read',access='stream')
    read(11) vp(:,:nplocal)
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc(1:nt,1:nt,1:nt,:,:,:)
    close(11)
    open(11,file=output_name('vc'),status='old',action='read',access='stream')
    read(11) vfield(:,1:nt,1:nt,1:nt,:,:,:)
    close(11)
    ! read PID & EID
    sync all
    print*, '  sum rhoc =', sum(rhoc)
    call buffer_np(rhoc)
    call buffer_vc(vfield)
    call redistribute_cdm
    call buffer_xp
    call buffer_vp
    sync all
    print*, '  sum rhoc =', sum(rhoc)
    print*, ''

    ! Open halo file
    open(11,file=output_name('halo'),status='replace',access='stream')
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
      rhof=0; isortpeak=0; den_peak=0;
      n_candidate=0; n_candidate_real=0;
      do k0=-1,nt+2 ! 2 more coarse cell layers
      do j0=-1,nt+2
      do i0=-1,nt+2
        nlast=cum(i0-1,j0,k0,itx,ity,itz)
        np=rhoc(i0,j0,k0,itx,ity,itz)
        do l0=1,np
          ip=nlast+l0
          if (NGP) then
            pos1=ncell*([i0,j0,k0]-1) + ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
            idx1=floor(pos1)+1
            rhof(idx1(1),idx1(2),idx1(3))=rhof(idx1(1),idx1(2),idx1(3))+mass_p
          else
            pos1=ncell*([i0,j0,k0]-1) + ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution - 0.5
            idx1=floor(pos1)+1
            idx2=idx1+1
            dx1=idx1-pos1
            dx2=1-dx1
            rhof(idx1(1),idx1(2),idx1(3))=rhof(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
            rhof(idx2(1),idx1(2),idx1(3))=rhof(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
            rhof(idx1(1),idx2(2),idx1(3))=rhof(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
            rhof(idx1(1),idx1(2),idx2(3))=rhof(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
            rhof(idx1(1),idx2(2),idx2(3))=rhof(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
            rhof(idx2(1),idx1(2),idx2(3))=rhof(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
            rhof(idx2(1),idx2(2),idx1(3))=rhof(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
            rhof(idx2(1),idx2(2),idx2(3))=rhof(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p
          endif
        enddo
      enddo
      enddo
      enddo

      ! find maxima
      print*, '  find_halo_candidates'
      do k0=1-ncell,nft+ncell  ! 1 more coarse cell layer
      do j0=1-ncell,nft+ncell
      do i0=1-ncell,nft+ncell
        denmax=maxval(rhof(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1))
        if (denmax==rhof(i0,j0,k0).and.denmax>den_peak_cutoff) then
          if (n_candidate>max_maxima-2) stop 'too many maxima'
          ! Find the fine mesh mass of this peak
          amtot=0
          do ii=1,irtot ! keep adding mass until mean density is less than halo_odc
            ix=i0+idist(1,ii); if(ix<5-nfb .or. ix>nft+nfb-4) cycle
            iy=j0+idist(2,ii); if(iy<5-nfb .or. iy>nft+nfb-4) cycle
            iz=k0+idist(3,ii); if(iz<5-nfb .or. iz>nft+nfb-4) cycle
            amtot=amtot+rhof(ix,iy,iz)
            if (complete_shell .and. rdist(ii+1)==rdist(ii)) cycle
            if (ii>18 .and. amtot/ii<halo_odc) exit
          enddo
          ! Consider this a halo candidate if the fine mesh mass is large enough
          if (amtot>mass_p*min_halo_particles/2.) then
            n_candidate=n_candidate+1
            n_candidate_real=n_candidate_real+merge(1,0,minval([i0,j0,k0])>=1 .and. maxval([i0,j0,k0])<=nft)
            ipeak(:,n_candidate)=[i0,j0,k0]-0.5
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
        !physical_halo=(minval(hpos)>0 .and. maxval(hpos)<nft)

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
          nlast=cum(i0-1,j0,k0,itile(1),itile(2),itile(3))
          np=rhoc(i0,j0,k0,itile(1),itile(2),itile(3))
          do l0=1,np
            ip=nlast+l0
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
          odcj=0.75/pi*ii*mass_p/r_vir(ii)**3
          if (odcj<halo_vir) then
            odci=0.75/pi*(ii-1)*mass_p/r_vir(ii-1)**3
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
          odcj=0.75/pi*ii*mass_p/r_odc(ii)**3
          if (odcj<halo_odc) then
            odci=0.75/pi*(ii-1)*mass_p/r_odc(ii-1)**3
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

print*,'  r',rrefine,rsearch;stop


        if (i_vir==0 .or. i_odc==0) then ! check number of particles
          n_search_fail=n_search_fail+1
          print*,'  search fail:'
          print*, ii,i_vir,i_odc
          print*,iloc,np_vir,np_odc
          print*,'  r',rrefine,rsearch
          !stop
        endif
        !print*,'  new hpos =',halo_info%hpos
        !open(12,file='finegrid.dat',access='stream')
        !write(12) finegrid
        !close(12)
        !close(13)
        !print*,'stat'
        if (i_vir >= min_halo_particles) then
          hpart_vir(ilist_vir(1:i_vir))=1
          hpart_odc(ilist_odc(1:i_odc))=1
          if (physical_halo) then  ! if the maximum is in physical regions
            nhalo=nhalo+1
            halo_info%hpos=hpos+(itile-1)*nft
            halo_info%mass_vir=mass_p*i_vir
            halo_info%mass_odc=mass_p*i_odc
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
      enddo ! do iloc=n_candidate,1,-1
      print*, '  found nhalo =',nhalo
    enddo
    enddo
    enddo !! itz
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
  enddo ! do cur_checkpoint= 1,n_checkpoint
  sync all
  if (head) print*, 'halofinder done'


contains
  subroutine buffer_np(rhoc)
    use parameters
    implicit none
    save

    integer(4) rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
    integer(8),parameter :: unit8=1
    integer(8) itest
    if (head) then
      print*, 'buffer_np'
    endif
    !x
    rhoc(:0,:,:,1,:,:)=rhoc(nt-ncb+1:nt,:,:,nnt,:,:)[image1d(inx,icy,icz)]
    rhoc(:0,:,:,2:,:,:)=rhoc(nt-ncb+1:nt,:,:,:nnt-1,:,:)
    rhoc(nt+1:,:,:,nnt,:,:)=rhoc(1:ncb,:,:,1,:,:)[image1d(ipx,icy,icz)]
    rhoc(nt+1:,:,:,:nnt-1,:,:)=rhoc(1:ncb,:,:,2:,:,:)
    sync all
    !y
    rhoc(:,:0,:,:,1,:)=rhoc(:,nt-ncb+1:nt,:,:,nnt,:)[image1d(icx,iny,icz)]
    rhoc(:,:0,:,:,2:,:)=rhoc(:,nt-ncb+1:nt,:,:,1:nnt-1,:)
    rhoc(:,nt+1:,:,:,nnt,:)=rhoc(:,1:ncb,:,:,1,:)[image1d(icx,ipy,icz)]
    rhoc(:,nt+1:,:,:,:nnt-1,:)=rhoc(:,1:ncb,:,1:,2:,:)
    sync all
    !z
    rhoc(:,:,:0,:,:,1)=rhoc(:,:,nt-ncb+1:nt,:,:,nnt)[image1d(icx,icy,inz)]
    rhoc(:,:,:0,:,:,2:)=rhoc(:,:,nt-ncb+1:nt,:,:,:nnt-1)
    rhoc(:,:,nt+1:,:,:,nnt)=rhoc(:,:,1:ncb,:,:,1)[image1d(icx,icy,ipz)]
    rhoc(:,:,nt+1:,:,:,:nnt-1)=rhoc(:,:,1:ncb,:,:,2:)
    sync all
  endsubroutine

  subroutine buffer_vc(vfield)
    use parameters
    implicit none
    save

    real(4) vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    ! the following variables are introduced because
    ! gcc only allows <= 7 ranks in arrays
    real(4) vtransx(3,ncb,nt+2*ncb,nt+2*ncb,nnt,nnt)[*]
    real(4) vtransy(3,nt+2*ncb,ncb,nt+2*ncb,nnt,nnt)[*]
    real(4) vtransz(3,nt+2*ncb,nt+2*ncb,ncb,nnt,nnt)[*]
    if (head) print*, 'buffer_vc'
    !x
    vtransx=vfield(:,nt-ncb+1:nt,:,:,nnt,:,:)
    sync all
    vfield(:,:0,:,:,1,:,:)=vtransx(:,:,:,:,:,:)[image1d(inx,icy,icz)]
    sync all
    vfield(:,:0,:,:,2:,:,:)=vfield(:,nt-ncb+1:nt,:,:,:nnt-1,:,:)

    vtransx=vfield(:,1:ncb,:,:,1,:,:)
    sync all
    vfield(:,nt+1:,:,:,nnt,:,:)=vtransx(:,:,:,:,:,:)[image1d(ipx,icy,icz)]
    sync all
    vfield(:,nt+1:,:,:,:nnt-1,:,:)=vfield(:,1:ncb,:,:,2:,:,:)
    sync all

    !y
    vtransy=vfield(:,:,nt-ncb+1:nt,:,:,nnt,:)
    sync all
    vfield(:,:,:0,:,:,1,:)=vtransy(:,:,:,:,:,:)[image1d(icx,iny,icz)]
    sync all
    vfield(:,:,:0,:,:,2:,:)=vfield(:,:,nt-ncb+1:nt,:,:,1:nnt-1,:)

    vtransy=vfield(:,:,1:ncb,:,:,1,:)
    sync all
    vfield(:,:,nt+1:,:,:,nnt,:)=vtransy(:,:,:,:,:,:)[image1d(icx,ipy,icz)]
    sync all
    vfield(:,:,nt+1:,:,:,:nnt-1,:)=vfield(:,:,1:ncb,:,1:,2:,:)
    sync all

    !z
    vtransz=vfield(:,:,:,nt-ncb+1:nt,:,:,nnt)
    sync all
    vfield(:,:,:,:0,:,:,1)=vtransz(:,:,:,:,:,:)[image1d(icx,icy,inz)]
    sync all
    vfield(:,:,:,:0,:,:,2:)=vfield(:,:,:,nt-ncb+1:nt,:,:,:nnt-1)

    vtransz=vfield(:,:,:,1:ncb,:,:,1)
    sync all
    vfield(:,:,:,nt+1:,:,:,nnt)=vtransz(:,:,:,:,:,:)[image1d(icx,icy,ipz)]
    sync all
    vfield(:,:,:,nt+1:,:,:,:nnt-1)=vfield(:,:,:,1:ncb,:,:,2:)
    sync all
  endsubroutine

  subroutine redistribute_cdm()
    implicit none
    save
    integer(8),parameter :: unit8=1
    integer(8) nshift,nlen,ifrom,checkxp0,checkxp1,iy,iz
    real overhead_image[*]

    if (head) then
      print*,''
      print*, 'redistribute_cdm'
    endif
    ! check
    overhead_image=sum(rhoc*unit8)/real(np_image_max,8)
    sync all
    do i0=1,nn**3
      overhead_image=max(overhead_image,overhead_image[i0])
    enddo
    if (head) then
      print*, '  image overhead',overhead_image*100,'% full'
      print*, '  comsumed image_buffer =',overhead_image*image_buffer,'/',image_buffer
    endif
    sync all

    if (overhead_image>1d0) then
      print*, '  error: too many particles in this image+buffer'
      print*, '  ',sum(rhoc*unit8),np_image_max
      print*, '  on',this_image()
      print*, '  please set image_buffer larger'
      stop
    endif

    ! shift to right
    checkxp0=sum(xp*int(1,kind=8))
    nshift=np_image_max-nplocal
    xp(:,nshift+1:np_image_max)=xp(:,1:nplocal)
    xp(:,1:nshift)=0
    vp(:,nshift+1:np_image_max)=vp(:,1:nplocal)
    vp(:,1:nshift)=0
# ifdef PID
      pid(nshift+1:np_image_max)=pid(1:nplocal)
      pid(1:nshift)=0
# endif
    checkxp1=sum(xp*int(1,kind=8))
    print*, '  ',np_image_max,nplocal,nshift
    if (checkxp0/=checkxp1) then
      print*, '  error in shifting right',image,checkxp0,checkxp1
      stop
    endif

    ! shift back
    cum=cumsum6(rhoc)
    ifrom=nshift

    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt ! loop over tiles
      do iz=1,nt ! loop over z slab
      do iy=1,nt ! loog over y slot
        ! nlast is the last particle's index
        nlast=cum(nt,iy,iz,itx,ity,itz)
        ! nlen is the number of particle in this x-slot
        nlen=nlast-cum(0,iy,iz,itx,ity,itz)
        xp(:,nlast-nlen+1:nlast)=xp(:,ifrom+1:ifrom+nlen)
        vp(:,nlast-nlen+1:nlast)=vp(:,ifrom+1:ifrom+nlen)
        xp(:,ifrom+1:ifrom+nlen)=0
        vp(:,ifrom+1:ifrom+nlen)=0
#     ifdef PID
          pid(nlast-nlen+1:nlast)=pid(ifrom+1:ifrom+nlen)
          pid(ifrom+1:ifrom+nlen)=0
#     endif
        ifrom=ifrom+nlen
      enddo
      enddo
    enddo
    enddo
    enddo

    checkxp1=sum(xp*int(1,kind=8))
    if (checkxp0/=checkxp1) then
      print*, '  error in shifting back',image,checkxp0,checkxp1
      stop
    endif
    sync all
  endsubroutine redistribute_cdm

  subroutine buffer_xp
    implicit none
    save
    integer(8) nshift,nlen,ifrom,mlast,iy,iz

    if (head) print*, 'buffer_xp'
    ! buffer x direction
    ! sync x- buffer with node on the left
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,1 ! do only tile_x=1
      do iz=1,nt
      do iy=1,nt
        nlast=cum(0,iy,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,iy,iz,itx,ity,itz)+rhoc(1-ncb,iy,iz,itx,ity,itz)
        mlast=cum(nt,iy,iz,nnt,ity,itz)[image1d(inx,icy,icz)]
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute local x- buffer
    do itz=1,nnt
    do ity=1,nnt
    do itx=2,nnt ! skip tile_x=1
      do iz=1,nt
      do iy=1,nt
        nlast=cum(0,iy,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,iy,iz,itx,ity,itz)+rhoc(1-ncb,iy,iz,itx,ity,itz)
        mlast=cum(nt,iy,iz,itx-1,ity,itz)
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! sync x+ buffer with node on the right
    do itz=1,nnt
    do ity=1,nnt
    do itx=nnt,nnt ! do only tile_x=nnt
      do iz=1,nt
      do iy=1,nt
        nlast=cum(nt+ncb,iy,iz,itx,ity,itz)
        nlen=nlast-cum(nt,iy,iz,itx,ity,itz)
        mlast=cum(ncb,iy,iz,1,ity,itz)[image1d(ipx,icy,icz)]
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute local x+ buffer
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt-1 ! skip tile_x=nnt
      do iz=1,nt
      do iy=1,nt
        nlast=cum(nt+ncb,iy,iz,itx,ity,itz)
        nlen=nlast-cum(nt,iy,iz,itx,ity,itz)
        mlast=cum(ncb,iy,iz,itx+1,ity,itz)
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! buffer y direction
    ! sync y-
    do itz=1,nnt
    do ity=1,1 ! do only ity=1
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,0,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc(1-ncb,1-ncb,iz,itx,ity,itz)
        mlast=cum(nt+ncb,nt,iz,itx,nnt,itz)[image1d(icx,iny,icz)]
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute y-
    do itz=1,nnt
    do ity=2,nnt ! skip ity=1
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,0,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc(1-ncb,1-ncb,iz,itx,ity,itz)
        mlast=cum(nt+ncb,nt,iz,itx,ity-1,itz)
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
      enddo
    enddo
    enddo
    enddo
    sync all

    ! sync y+
    do itz=1,nnt
    do ity=nnt,nnt ! do only ity=nnt
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,nt+ncb,iz,itx,ity,itz)
        nlen=nlast-cum(nt+ncb,nt,iz,itx,ity,itz)
        mlast=cum(nt+ncb,ncb,iz,itx,1,itz)[image1d(icx,ipy,icz)]
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute y+
    do itz=1,nnt
    do ity=1,nnt-1 ! skip ity=nnt
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,nt+ncb,iz,itx,ity,itz)
        nlen=nlast-cum(nt+ncb,nt,iz,itx,ity,itz)
        mlast=cum(nt+ncb,ncb,iz,itx,ity+1,itz)
        xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
      enddo
    enddo
    enddo
    enddo
    sync all

    ! buffer z direction
    ! sync z-
    do itz=1,1 ! do only itz=1
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,0,itx,ity,itz)
      nlen=nlast-cum(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc(1-ncb,1-ncb,1-ncb,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,nt,itx,ity,nnt)[image1d(icx,icy,inz)]
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
    enddo
    enddo
    enddo
    sync all

    ! redistribute z-
    do itz=2,nnt ! skip itz=1
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,0,itx,ity,itz)
      nlen=nlast-cum(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc(1-ncb,1-ncb,1-ncb,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,nt,itx,ity,itz-1)
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
    enddo
    enddo
    enddo
    sync all

    ! sync z+
    do itz=nnt,nnt ! do only itz=nnt
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
      nlen=nlast-cum(nt+ncb,nt+ncb,nt,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,ncb,itx,ity,1)[image1d(icx,icy,ipz)]
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
    enddo
    enddo
    enddo
    sync all

    ! redistribute z+
    do itz=1,nnt-1 ! skip itz=nnt
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
      nlen=nlast-cum(nt+ncb,nt+ncb,nt,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,ncb,itx,ity,itz+1)
      xp(:,nlast-nlen+1:nlast)=xp(:,mlast-nlen+1:mlast)
    enddo
    enddo
    enddo
    sync all

  endsubroutine buffer_xp

  subroutine buffer_vp
    !use variables
    implicit none
    save
    integer(8) nshift,nlen,ifrom,mlast

# ifdef PID
      if (head) print*, 'buffer_vp (vp & pid)'
# else
      if (head) print*, 'buffer_vp'
# endif
    ! buffer x direction
    ! sync x- buffer with node on the left
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,1 ! do only tile_x=1
      do iz=1,nt
      do iy=1,nt
        nlast=cum(0,iy,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,iy,iz,itx,ity,itz)+rhoc(1-ncb,iy,iz,itx,ity,itz)
        mlast=cum(nt,iy,iz,nnt,ity,itz)[image1d(inx,icy,icz)]
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(inx,icy,icz)]
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute local x- buffer
    do itz=1,nnt
    do ity=1,nnt
    do itx=2,nnt ! skip tile_x=1
      do iz=1,nt
      do iy=1,nt
        nlast=cum(0,iy,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,iy,iz,itx,ity,itz)+rhoc(1-ncb,iy,iz,itx,ity,itz)
        mlast=cum(nt,iy,iz,itx-1,ity,itz)
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! sync x+ buffer with node on the right
    do itz=1,nnt
    do ity=1,nnt
    do itx=nnt,nnt ! do only tile_x=nnt
      do iz=1,nt
      do iy=1,nt
        nlast=cum(nt+ncb,iy,iz,itx,ity,itz)
        nlen=nlast-cum(nt,iy,iz,itx,ity,itz)
        mlast=cum(ncb,iy,iz,1,ity,itz)[image1d(ipx,icy,icz)]
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(ipx,icy,icz)]
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute local x+ buffer
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt-1 ! skip tile_x=nnt
      do iz=1,nt
      do iy=1,nt
        nlast=cum(nt+ncb,iy,iz,itx,ity,itz)
        nlen=nlast-cum(nt,iy,iz,itx,ity,itz)
        mlast=cum(ncb,iy,iz,itx+1,ity,itz)
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all

    ! buffer y direction
    ! sync y-
    do itz=1,nnt
    do ity=1,1 ! do only ity=1
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,0,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc(1-ncb,1-ncb,iz,itx,ity,itz)
        mlast=cum(nt+ncb,nt,iz,itx,nnt,itz)[image1d(icx,iny,icz)]
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,iny,icz)]
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute y-
    do itz=1,nnt
    do ity=2,nnt ! skip ity=1
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,0,iz,itx,ity,itz)
        nlen=nlast-cum(1-ncb,1-ncb,iz,itx,ity,itz)+rhoc(1-ncb,1-ncb,iz,itx,ity,itz)
        mlast=cum(nt+ncb,nt,iz,itx,ity-1,itz)
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
      enddo
    enddo
    enddo
    enddo
    sync all

    ! sync y+
    do itz=1,nnt
    do ity=nnt,nnt ! do only ity=nnt
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,nt+ncb,iz,itx,ity,itz)
        nlen=nlast-cum(nt+ncb,nt,iz,itx,ity,itz)
        mlast=cum(nt+ncb,ncb,iz,itx,1,itz)[image1d(icx,ipy,icz)]
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,ipy,icz)]
      enddo
    enddo
    enddo
    enddo
    sync all

    ! redistribute y+
    do itz=1,nnt
    do ity=1,nnt-1 ! skip ity=nnt
    do itx=1,nnt
      do iz=1,nt
        nlast=cum(nt+ncb,nt+ncb,iz,itx,ity,itz)
        nlen=nlast-cum(nt+ncb,nt,iz,itx,ity,itz)
        mlast=cum(nt+ncb,ncb,iz,itx,ity+1,itz)
#   ifdef PID
        pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
        vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
      enddo
    enddo
    enddo
    enddo
    sync all

    ! buffer z direction
    ! sync z-
    do itz=1,1 ! do only itz=1
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,0,itx,ity,itz)
      nlen=nlast-cum(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc(1-ncb,1-ncb,1-ncb,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,nt,itx,ity,nnt)[image1d(icx,icy,inz)]
# ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
# endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,icy,inz)]
    enddo
    enddo
    enddo
    sync all

    ! redistribute z-
    do itz=2,nnt ! skip itz=1
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,0,itx,ity,itz)
      nlen=nlast-cum(1-ncb,1-ncb,1-ncb,itx,ity,itz)+rhoc(1-ncb,1-ncb,1-ncb,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,nt,itx,ity,itz-1)
# ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
# endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
    enddo
    enddo
    enddo
    sync all

    ! sync z+
    do itz=nnt,nnt ! do only itz=nnt
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
      nlen=nlast-cum(nt+ncb,nt+ncb,nt,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,ncb,itx,ity,1)[image1d(icx,icy,ipz)]
# ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
# endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)[image1d(icx,icy,ipz)]
    enddo
    enddo
    enddo
    sync all

    ! redistribute z+
    do itz=1,nnt-1 ! skip itz=nnt
    do ity=1,nnt
    do itx=1,nnt
      nlast=cum(nt+ncb,nt+ncb,nt+ncb,itx,ity,itz)
      nlen=nlast-cum(nt+ncb,nt+ncb,nt,itx,ity,itz)
      mlast=cum(nt+ncb,nt+ncb,ncb,itx,ity,itz+1)
#   ifdef PID
      pid(nlast-nlen+1:nlast)=pid(mlast-nlen+1:mlast)
#   endif
      vp(:,nlast-nlen+1:nlast)=vp(:,mlast-nlen+1:mlast)
    enddo
    enddo
    enddo
    sync all
  endsubroutine buffer_vp


  function cumsum6(rho_input)
    implicit none
    integer(4) rho_input(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8) cumsum6(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8) nsum,ihx,ihy,ihz,igx,igy,igz
    nsum=0
    do ihz=1,nnt
    do ihy=1,nnt
    do ihx=1,nnt
    do igz=1-ncb,nt+ncb
    do igy=1-ncb,nt+ncb
    do igx=1-ncb,nt+ncb
      nsum=nsum+rho_input(igx,igy,igz,ihx,ihy,ihz)
      cumsum6(igx,igy,igz,ihx,ihy,ihz)=nsum
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  endfunction

end
!! test indexedsort
!integer n
!integer iarr(4)
!real rarr(4)
!n=4
!iarr=(/1,2,3,4/)
!rarr=(/40,20,30,10/)
!call indexedsort(n,rarr(1:n),iarr(1:n))
!print*, iarr
!print*, rarr
!!!!!!!!!!!!!!!!!!!
