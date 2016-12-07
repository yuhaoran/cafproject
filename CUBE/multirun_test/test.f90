program test

Implicit None
integer seed(12), i, n
integer(4) :: t
real a(5)
call random_seed(size=n) ! get the size of random seed
call system_clock(t)
  do i = 1, n
      seed(i) = t+i
      !print*,t,seed(i)
  end do

call random_seed(put=seed) ! generate seed using system time
call random_number(a)
print*,"Test1"
print*,a

endprogram