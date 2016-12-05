module penfft_config
  use parameters
  implicit none
  save
  
  !! parameters !!
  ! integer,parameter  nn ! in parameters.f90

! ng is the fft grid number /node/dim

#ifdef penfft_fine
  integer,parameter :: ng=nf
#else

# ifdef penfft_1x
    integer,parameter :: ng=nc
# elif penfft_2x
    integer,parameter :: ng=nc*2
# elif penfft_4x
    integer,parameter :: ng=nc*4
# elif penfft_8x
    integer,parameter :: ng=nc*8
# elif penfft_32
    integer,parameter :: ng=32
# endif

#endif

  integer,parameter :: ngpen=ng/nn
  
contains



endmodule
