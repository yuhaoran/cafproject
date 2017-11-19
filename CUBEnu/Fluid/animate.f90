!mkdir Movie
!convert -delay 20 -loop 10 Movie/u1*.pgm u1.gif
  integer, parameter :: n=100,nt=1000,ny=100
  real, dimension(n,nt) :: u0,u1
  real, dimension(n,ny) :: rmap
  character*3 fn3
  character*1 fn1

  do ifl=0,1
     write(fn1,'(i1)'), ifl
     open(10,file='u'//fn1//'.dat',access='stream')
     read(10) u0
     rmax=maxval(u0)
     rmin=minval(u0)
     write(*,*) rmin,rmax
     do it=1,nt
        write(fn3,'(i3.3)') it-1
        rmap=0
        do i=1,n
           iy=ny-nint((u0(i,it)-rmin)*ny/(rmax-rmin))
           if (iy<1) write(*,*) i,it,iy,u0(i,it)
           iy=max(1,iy)
           iy=min(ny,iy)
           rmap(i,iy)=1
        end do
        call pmap('Movie/u'//fn1//'.'//fn3//'.pgm',rmap,n,ny,1)        
     end do
  end do
contains
  subroutine pmap(fn,rmap1,nx,ny,iscale)
  real rmap(nx,ny),rmap1(nx,ny)
  integer*2, dimension(nx,ny) :: imap
  integer*1, dimension(nx,ny) :: imap1
  character(len=*):: fn
  integer npix,mypos
  real rmax,rmin

  npix=min(ny/2-1,nx/2-1,300)
  
  
  rmap=rmap1
  iscale1=iscale
  do while (iscale1 > 1)      
     rmap=sign((sqrt(abs(rmap))),rmap)
     iscale1=iscale1-1
  end do
  rmax=maxval(rmap)
  rmin=minval(rmap)
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
