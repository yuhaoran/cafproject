rm -f *.o *.out

#mpif90 -O3 -cpp -fcoarray=single -mcmodel=medium test.f90 -lfftw3f -lm -ldl

mpif90 -O3 -cpp -fcoarray=single test.f90 -lfftw3f -lm -ldl

# work at least for nc=768

# on macbook pro, adding -mcmodel=medium causes problems:
#/var/folders/0g/jjmy2nk55731vdwm6mjtqz_m0000gn/T//ccmIq8rz.s:457:2: error: unsupported symbol modifier in relocation
#        movabsq $535822336+_equiv.1.3440@GOTOFF, %rax

