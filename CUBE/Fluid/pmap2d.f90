  integer, parameter :: nx=512,nz=nx,nt=426
  real, dimension(nx,nz) :: u0
  character*3 fn3

  open(12,file='u0.dat',access='stream')
  rmax=0
  rmin=1000
  do it=1,nt
  read(12) u0
  rmax=max(rmax,maxval(u0))
  rmin=min(rmin,minval(u0))
  enddo
  rewind(12)
  write(*,*) rmax,rmin
  do it=1,nt
  write(fn3,'(i3.3)')it-1
  read(12) u0
  call pmap('Movie/u0.'//fn3//'.pgm',u0,nx,nz,1)
  enddo

contains
  subroutine pmap(fn,rmap1,nx,ny,iscale)
  real rmap(nx,ny),rmap1(nx,ny)
  integer*2, dimension(nx,ny) :: imap
  integer*1, dimension(nx,ny) :: imap1
  character(len=*):: fn
  integer npix,mypos

  npix=min(ny/2-1,nx/2-1,300)
  
  
  rmap=rmap1
  iscale1=iscale
  do while (iscale1 > 1)      
     rmap=sign((sqrt(abs(rmap))),rmap)
     iscale1=iscale1-1
  end do
!  rmax=maxval(rmap)
!  rmin=minval(rmap)
  write(*,*) trim(fn),rmax,rmin
  imap=255*(rmap-rmin)/(rmax-rmin)
  imap1=127*(rmap-rmin)/(rmax-rmin)
  open(10,file=fn)
  write(10,'(2hP5)')
  write(10,*)nx,ny
  write(10,*) 255
!  write(10,*) 127
  INQUIRE(UNIT=10, POS=mypos)
  close(10)
  open(10,file=fn, access='stream',position='append')
!  write(10,pos=mypos) int(imap,1)
  write(10) int(imap,1)
  close(10)
end subroutine pmap

end program
