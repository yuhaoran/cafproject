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

  if (head) print*, 'halofind'

  open(10,file=output_name_halo('halo'),status='replace',access='stream')
  write(10) halo_info
  close(10)



  cur_halofind=cur_halofind+1
  halofind_step=.false.

  sync all


endsubroutine
