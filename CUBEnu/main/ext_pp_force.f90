module pp_force
  use variables
  !use neutrinos
  implicit none
  save


contains

subroutine ext_pp_force
  use omp_lib
  use variables
  !use neutrinos
  implicit none
  save

  integer(8) ntest,ip_offset,ipll1,ipll2

  if (head) then
    print*, ''
    print*, 'ext_pp_force'
    print*, '  pp_range=',int(pp_range,1)
    print*, '  ext_pp_force over',int(nnt**3,2),'tiles'
    print*, '  np_pp_max=',np_pp_max
  endif
  nth=omp_get_max_threads()
  print*,'  max num_threads =',nth

  itest1=0
  npairs=0
  f2_max_pp=0
  print*, '  sim%nplocal =',sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    if (head) print*, '  tile:',int(itx,1),int(ity,1),int(itz,1)
    call system_clock(t1,t_rate)
    if (np_pp_max<sum(rhoc(0:nt+1,0:nt+1,0:nt+1,itx,ity,itz))) then
      print*, 'np_pp_max too small'
      print*, np_pp_max, sum(rhoc(0:nt+1,0:nt+1,0:nt+1,itx,ity,itz))
      stop
    endif
    ! the offset to use tile-based linked-list
    ip_offset=cum(-1,0,0,itx,ity,itz)
    !print*,'  ip_offset',ip_offset
    hoc=0; ll=0
    do igz=0,nt+1
    do igy=0,nt+1
    do igx=0,nt+1
      nlast=cum(igx-1,igy,igz,itx,ity,itz)
      np=rhoc(igx,igy,igz,itx,ity,itz)
      do lp=1,np ! loop over cdm particles
        ip1=nlast+lp
        xvec1=ncell*((/igx,igy,igz/)-1)+ncell*(int(xp(:,ip1)+ishift,izipx)+rshift)*x_resolution
        ivec1=floor(xvec1)+1
        ipll1=ip1-ip_offset
        ll(ipll1)=hoc(ivec1(1),ivec1(2),ivec1(3))
        hoc(ivec1(1),ivec1(2),ivec1(3))=ipll1
      enddo
    enddo
    enddo
    enddo

    !$omp paralleldo &
    !$omp& default(shared) &
    !$omp& private(igz,igy,igx,ip1,xvec1,vreal,f_tot,kk,jj,ii,ip2) &
    !$omp& private(xvec2,xvec21,rmag,rcut,pcut,force_pp,ntest) &
    !$omp& reduction(+:itest1,npairs) &
    !$omp& reduction(max:f2_max_pp)
    do igz=1,nft
    do igy=1,nft
    do igx=1,nft
      ipll1=hoc(igx,igy,igz)
      do while (ipll1/=0)
        itest1=itest1+1
        ip1=ipll1+ip_offset
        xvec1=ncell*floor(((/igx,igy,igz/)-1)/4.)+ncell*(int(xp(:,ip1)+ishift,izipx)+rshift)*x_resolution
        vreal=tan(pi*real(vp(:,ip1))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        f_tot=0; ntest=0
        do kk=igz-pp_range,igz+pp_range
        do jj=igy-pp_range,igy+pp_range
        do ii=igx-pp_range,igx+pp_range
          ipll2=hoc(ii,jj,kk)
          do while (ipll2/=0) ! loop over particle 2
            ntest=ntest+1
            npairs=npairs+1
            ip2=ipll2+ip_offset
            xvec2=ncell*floor(((/ii,jj,kk/)-1)/4.)+ncell*(int(xp(:,ip2)+ishift,izipx)+rshift)*x_resolution
            xvec21=xvec2-xvec1
            rmag=sqrt(sum(xvec21**2))
            rmag=merge(1d0,rmag,rmag==0)
            rcut=rmag/nf_cutoff
            pcut=1-(7./4*rcut**3)+(3./4*rcut**5)
            force_pp=sim%mass_p_cdm*(xvec21/rmag**3)*pcut
            force_pp=merge(force_pp,force_pp*0,rmag>rsoft)
            f_tot=f_tot+force_pp
            ipll2=ll(ipll2)
          enddo !! do while (ip2/=0)
        enddo !! kk
        enddo !! jj
        enddo !! ii
        vreal=vreal+f_tot*a_mid*dt/6/pi
#       ifdef FORCETEST
          print*,'v_pp',ip,vreal
          print*,'vfield',vfield(:,i,j,k,itx,ity,itz)
#       endif
        vp(:,ip1)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
        f2_max_pp=max(f2_max_pp,sum(f_tot**2))
        ipll1=ll(ipll1)
      enddo !! do while (ip1/=0)
    enddo
    enddo
    enddo
    !$omp endparalleldo
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, '  itest1 =',itest1
  enddo
  enddo
  enddo !! itz
  ! for reference, in pm.f90
  !dt_fine=sqrt( 1.0 / (sqrt(maxval(f2_max_fine))*a_mid*GG) )
  !dt_coarse=sqrt( real(ncell) / (sqrt(f2_max_coarse)*a_mid*GG) )
  !dt_pp=sqrt(0.1*rsoft) / max(sqrt(maxval(f2_max_pp))*a_mid*GG,1e-3)
  sim%dt_pp=sqrt(1.0) / max(sqrt(f2_max_pp)*a_mid*GG,1e-3)
  sync all
  do i=1,nn**3
    sim%dt_pp=min(sim%dt_pp,sim[i]%dt_pp)
  enddo
  sync all

  if (head) then
    print*,'  max of f_pp', sqrt(f2_max_pp)
    print*,'  dt_pp',sim%dt_pp
    print*,'  updated',itest1,'particles'
    print*,'  average pairs',npairs/real(itest1)
  endif
  sync all

  if (itest1/=sim%nplocal) then
    print*, 'itest1/=sim%nplocal'
    print*, itest1,sim%nplocal
    stop
  endif
  sync all
endsubroutine


endmodule
