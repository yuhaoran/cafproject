program ps

use iso_fortran_env , only : int64
Implicit None
character(len=250)  filename
integer nst, Nd, Nd21, boxsize, Nd2, Klen
real(8) PI

Nd=500
Nd21=Nd/2+1
Nd2=Nd/2
boxsize=500 ! Mpc/h
PI = 3.141592654
Klen = Nd*Nd*Nd21*2;

filename = "DeltaKi_1500" ! input the filename of the deltaK
call read_deltaki()
call Initial_GFKXYZ()

nst=1500
call cal_numps()

contains

Subroutine cal_numps()

  integer i,j,k,id
  real(8) kk,volume,hre,him
  real(8) pi2L,Nor,tt1,tt2, ipp
  real(8) pssim(3, 10000), DeltaKi(Klen), K3d(Nd)
  integer  nps(10000), pid, stat
  character(len=6) nsts

  print*, "Now output the power spectrum and Hamilton forces for %d initial condition", nst
   
   open (unit = 2, file = filename, action='read', iostat=stat,  status = 'old', access = 'stream', form='unformatted')
   if (stat .ne. 0) then
       print*, "can't open file `%s'\n", filename
       stop 1
   end if
   read(2) nst
   read(2) ipp
   do i = 1, Nd*Nd*Nd21
    read(2) DeltaKi(2*i-1), DeltaKi(2*i)
   end do
   close(2)

  pssim = 0.0
  nps = 0
  pi2L = 2*PI/boxsize  
  volume = boxsize**3
  do i = 1, Nd
   do j = 1, Nd
    do k = 1, Nd21
      id=((i-1)*Nd+(j-1))*Nd21+k
      kk=sqrt(K3d(i)*K3d(i)+K3d(j)*K3d(j)+K3d(k)*K3d(k))
      pid=(kk/pi2L) + 1
      if(pid < Nd2)then
       tt1=DeltaKi(id*2-1)*DeltaKi(id*2-1)+DeltaKi(id*2)*DeltaKi(id*2)
       if(k == 1)then
        pssim(1, pid)=pssim(1, pid)+kk
        pssim(2, pid)=pssim(2, pid)+tt1
        nps(pid)=nps(pid)+1
       else
        pssim(1, pid)=pssim(1, pid)+2*kk
        pssim(2, pid)=pssim(2, pid)+2*tt1
        nps(pid)=nps(pid)+2
       end if
      end if  
    end do
   end do
  end do

  write( nsts,'(i6)' ) nst
  filename = 'powspe/PSHF_'//trim(adjustl(nsts))//'.dat'
  
  open (unit = 1, file = filename, action='write', iostat=stat,  status = 'replace')
  if (stat .ne. 0) then
      print*, "can't open file `%s'\n", filename
      stop 1
  end if

  do i = 2, Nd21
      pssim(1, i)=pssim(1, i)/nps(i)
      pssim(2, i)=pssim(2, i)/nps(i)
      !tt1=PowerSpe_Z0(pssim[i][0])	
      write(1,*) pssim(1, i), pssim(2, i)*volume, tt1, nps(i)
  end do
  close(1) 

  print*, "output power spectrum and Hamilton forces into file %s done", filename

endsubroutine cal_numps

Subroutine read_deltaki()

   integer nst,i,j,ccid, stat
   real(8) ipp, DeltaKi(Klen)
   print*,  "Now read DeltaKi data from previously dumped file:", filename
   
   open (unit = 2, file = filename, action='read', iostat=stat,  status = 'old', access = 'stream', form='unformatted')
   if (stat .ne. 0) then
       print*, "can't open file `%s'\n", filename
       stop 1
   end if

   read(2) nst
   read(2) ipp
   do i = 1, Nd*Nd*Nd21
    read(2) DeltaKi(2*i-1), DeltaKi(2*i)
   end do
   close(2)
   
   print*, "Ip is ", ipp, "nstep is ", nst

endsubroutine read_deltaki

Subroutine Initial_GFKXYZ()
   integer i;
   real(8) pi2L;
   real(8) K3d(Nd);
   pi2L=2*PI/boxsize;
   do i = 1, Nd
      K3d(i)=(i-1)*pi2L;
      if(i-1 > Nd2)then
     	K3d(i)=(i-1-Nd)*pi2L;
      end if
   end do
endsubroutine Initial_GFKXYZ

endprogram