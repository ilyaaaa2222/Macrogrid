program macrogrid_solver
   implicit none
   integer, parameter :: io = 10

   integer, parameter :: SUBGRID_SIZES_COUNT = 8
   integer, dimension(SUBGRID_SIZES_COUNT), parameter :: subgrid_sizes = &
      [10, 18, 34, 66, 130, 258, 514, 1026]

   real*8, dimension(SUBGRID_SIZES_COUNT), parameter :: factors_const = &
      [1.52729d0, 1.69504d0, 1.82964d0, 1.90932d0, 1.95305d0, 1.97616d0, 1.98798d0, 1.99394d0]

   real*8, dimension(SUBGRID_SIZES_COUNT), parameter :: factors_log = &
      [1.52129d0, 1.69524d0, 1.82919d0, 1.90899d0, 1.95307d0, 1.97618d0, 1.98793d0, 1.99403d0]

   integer, parameter :: SUBGRID_CFGS_COUNT = 5
   integer, dimension(SUBGRID_CFGS_COUNT), parameter :: subgrid_configs = &
      [1, 2, 3, 4, 5]

   integer, parameter :: MACRIGRID_CFGS_COUNT = 2
   integer, dimension(2, MACRIGRID_CFGS_COUNT), parameter :: macrogrid_configs = reshape([ &
      2, 2, &
      4, 4 &
      ], shape=[2, MACRIGRID_CFGS_COUNT])

   integer, parameter :: max_iter_subgrid   = 100000
   integer, parameter :: max_iter_interface = 100000
   real*8,  parameter :: eps_subgrid        = 1.0d-8
   real*8,  parameter :: eps_interface      = 1.0d-8

   real*8 :: macrogrid(10, 10, 1026, 1026)

   integer :: i_cfg, i_size, cfg_x, cfg_y, sub_sz, iter
   real*8  :: omega, time, error

   external :: simple_iteration, conjugate_residuals,simple_iteration_withBoost
   external :: original_sor, tiling_sor, subtiling_sor, subtiling_sor_test_version
   external ::  subtiling_sor_8
   external :: initialize_constant_boundary, initialize_logarithmic_boundary
   external :: compute_constant_boundary_error, compute_logarithmic_boundary_error

   open(io, file='test.txt', status='unknown', action='write')

   write(io, *) "simple_iteration_one_iter_withBoost with original_sor"
   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)
      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_sizes(subgrid_configs(i_size))
         omega = factors_log(subgrid_configs(i_size))
         call run_test(io, simple_iteration_withBoost, original_sor, initialize_logarithmic_boundary, &
            compute_logarithmic_boundary_error, &
            macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, 1, max_iter_interface)
      end do
   end do
   write(io, *) " "
   write(*, *) "hui"
   write(io, *) "simple_iteration_one_iter with original_sor"
   do i_cfg = 1, MACRIGRID_CFGS_COUNT
      cfg_x = macrogrid_configs(1, i_cfg)
      cfg_y = macrogrid_configs(2, i_cfg)
      do i_size = 1, SUBGRID_CFGS_COUNT
         sub_sz = subgrid_sizes(subgrid_configs(i_size))
         omega = factors_log(subgrid_configs(i_size))
         call run_test(io, simple_iteration, original_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface,1, max_iter_interface)
      end do
   end do
   !write(io, *) " "
   !write(*, *) "hui"
   !write(io, *) "simple_iteration_one_iter with tiling_sor"
   !do i_cfg = 1, MACRIGRID_CFGS_COUNT
      !cfg_x = macrogrid_configs(1, i_cfg)
      !cfg_y = macrogrid_configs(2, i_cfg)
      !do i_size = 1, SUBGRID_CFGS_COUNT
         !sub_sz = subgrid_sizes(subgrid_configs(i_size))
         !omega = factors_log(subgrid_configs(i_size))
         !call run_test(io, simple_iteration, tiling_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            !macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, 1, max_iter_interface)
      !end do
   !end do
   !write(io, *) " "
   !write(*, *) "hui"
   !write(io, *) "conjugate_residuals with original_sor"
   !do i_cfg = 1, MACRIGRID_CFGS_COUNT
      !cfg_x = macrogrid_configs(1, i_cfg)
      !cfg_y = macrogrid_configs(2, i_cfg)
      !do i_size = 1, SUBGRID_CFGS_COUNT
         !sub_sz = subgrid_sizes(subgrid_configs(i_size))
         !omega = factors_log(subgrid_configs(i_size))
         !call run_test(io, conjugate_residuals, original_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            !macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, max_iter_subgrid, max_iter_interface)
      !end do
   !end do
   !write(io, *) " "
   !write(*, *) "hui"
   !write(io, *) "conjugate_residuals with tiling_sor"
   !do i_cfg = 1, MACRIGRID_CFGS_COUNT
      !cfg_x = macrogrid_configs(1, i_cfg)
      !cfg_y = macrogrid_configs(2, i_cfg)
      !do i_size = 1, SUBGRID_CFGS_COUNT
         !sub_sz = subgrid_sizes(subgrid_configs(i_size))
         !omega = factors_log(subgrid_configs(i_size))
         !call run_test(io, conjugate_residuals, tiling_sor, initialize_logarithmic_boundary, compute_logarithmic_boundary_error, &
            !macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, max_iter_subgrid, max_iter_interface)
      !end do
   !end do
   !write(io, *) " "
   !write(*, *) "hui"
   !write(io, *) "conjugate_residuals with subtiling_sor_8"
   !do i_cfg = 1, MACRIGRID_CFGS_COUNT
      !cfg_x = macrogrid_configs(1, i_cfg)
      !cfg_y = macrogrid_configs(2, i_cfg)
      !do i_size = 1, SUBGRID_CFGS_COUNT
         !sub_sz = subgrid_sizes(subgrid_configs(i_size))
         !if (sub_sz == 10) then
            !cycle
         !end if
         !omega = factors_log(subgrid_configs(i_size))
         !call run_test(io, conjugate_residuals, subtiling_sor_8, initialize_logarithmic_boundary, &
            !compute_logarithmic_boundary_error, &
            !macrogrid, cfg_x, cfg_y, sub_sz, omega, eps_subgrid, eps_interface, max_iter_subgrid, max_iter_interface)
      !end do
   !end do
   !write(io, *) " "
   !write(*, *) "hui"

   close(io)

end program macrogrid_solver

subroutine run_test(io, macrigrid_solver, subgrid_solver, &
   boundary_condition, compute_error,                 &
   macrogrid, macrogrid_size_x, macrogrid_size_y,     &
   subgrid_size, omega,                               &
   eps_subgrid, eps_interface,                        &
   max_iter_subgrid, max_iter_interface)
   implicit none
   external :: macrigrid_solver, subgrid_solver, boundary_condition, compute_error
   integer, intent(in) :: io, macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter_subgrid, max_iter_interface
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(in)   :: omega, eps_subgrid, eps_interface
   real*8   :: error, time
   integer  :: iter

   call boundary_condition(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   call macrigrid_solver(subgrid_solver, macrogrid, macrogrid_size_x, macrogrid_size_y, &
      subgrid_size, omega, eps_subgrid, eps_interface, &
      max_iter_subgrid, max_iter_interface, error, time, iter)
   call compute_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, error)

   write(io,'(A,I5,A,I5,A,I5,A,F14.8,A,F14.8,A,I6)') &
      "Macrogrid_size_x&Macrogrid_size_y&Subgrid_size&Error&Run_time&Iterations#", &
      macrogrid_size_x, "&", macrogrid_size_y, "&", subgrid_size, "&", error, "&", time, "&", iter

end subroutine run_test

subroutine print_macrogrid(macrogrid, macrogrid_size, subgrid_size)
   implicit none
   integer, intent(in) :: macrogrid_size, subgrid_size
   real*8, intent(in) :: macrogrid(macrogrid_size, macrogrid_size, subgrid_size, subgrid_size)
   integer :: i, j, x, y

   do i = 1, macrogrid_size
      do y = 1, subgrid_size
         do j = 1, macrogrid_size
            do x = 1, subgrid_size
               write(*, '(10F8.3)', advance='no') macrogrid(j, i, x, y)
            end do
         end do
         print *, ' '
      end do
   end do
   print *, ' '

end subroutine print_macrogrid

subroutine compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(out) :: dx, dy
   integer :: global_size_x
   integer :: global_size_y
   global_size_x = macrogrid_size_x*subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y*subgrid_size - (macrogrid_size_y - 1)
   dx = 1.0d0 / dble(global_size_x - 1)
   dy = 1.0d0 / dble(global_size_y - 1)
end subroutine compute_dxdy

subroutine get_interface(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, u)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size
   real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(out) :: u(interface_size)

   integer :: iX, iY, i1, j1, k, i

   i = 1
   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y - 1
         do k = 2, subgrid_size - 1

            u(i) = macrogrid(iX, iY, k, subgrid_size)
            i = i + 1

         end do
      end do
   end do

   do iX = 1, macrogrid_size_x - 1
      do iY = 1, macrogrid_size_y
         do k = 2, subgrid_size - 1

            u(i) = macrogrid(iX, iY, subgrid_size, k)
            i = i + 1

         end do
      end do
   end do

end subroutine get_interface

subroutine set_interface(macrogrid, u, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size
   real*8,  intent(in) :: u(interface_size)
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   integer :: iX, iY, i1, j1, k, i

   i = 1
   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y - 1
         do k = 2, subgrid_size - 1

            macrogrid(iX, iY, k, subgrid_size) = u(i)
            macrogrid(iX, iY+1, k, 1) = u(i)
            i = i + 1

         end do
      end do
   end do

   do iX = 1, macrogrid_size_x - 1
      do iY = 1, macrogrid_size_y
         do k = 2, subgrid_size - 1

            macrogrid(iX, iY, subgrid_size, k) = u(i)
            macrogrid(iX+1, iY, 1, k) = u(i)
            i = i + 1

         end do
      end do
   end do

end subroutine set_interface

subroutine s(macrogrid, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, s_vec)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size
   real*8,  intent(in) :: dx, dy
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(out) :: s_vec(interface_size)

   integer :: iX, iY, i1, j1, k, i

   i = 1
   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y - 1
         do k = 2, subgrid_size - 1

            s_vec(i) = (3.0d0*macrogrid(iX, iY, k, subgrid_size) - 4.0d0*macrogrid(iX, iY, k, subgrid_size-1 ) + &
               macrogrid(iX, iY, k, subgrid_size-2) - &
               (- 3.0d0*macrogrid(iX, iY+1, k, 1) + 4.0d0*macrogrid(iX, iY+1, k, 2 ) - macrogrid(iX, iY+1, k, 3 ))) / (2.0d0 * dy)
            i = i + 1

         end do
      end do
   end do

   do iX = 1, macrogrid_size_x - 1
      do iY = 1, macrogrid_size_y
         do k = 2, subgrid_size - 1

            s_vec(i) = (3.0d0*macrogrid(iX, iY, subgrid_size, k) - 4.0d0*macrogrid(iX,iY, subgrid_size-1, k) + &
               macrogrid(iX,iY, subgrid_size-2, k) - &
               (- 3.0d0*macrogrid(iX+1, iY, 1, k) + 4.0d0*macrogrid(iX+1, iY, 2, k) - macrogrid(iX+1, iY, 3, k))) / (2.0d0 * dx)
            i = i + 1

         end do
      end do
   end do

end subroutine s

subroutine Cv(subgrid_solver, macrogrid, v, b, dx, dy, macrogrid_size_x, macrogrid_size_y, &
   subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Cv_vec)
   implicit none
   external :: subgrid_solver
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, max_iter_subgrid
   real*8,  intent(in) :: v(interface_size), b(interface_size), dx, dy, omega, eps_subgrid
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(out) :: Cv_vec(interface_size)

   real*8 :: subgrid_error
   integer :: iX, iY, i1, j1, k, i

   call set_interface(macrogrid, v, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size)

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
            omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
      end do
   end do

   call s(macrogrid, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, Cv_vec)

   i = 1
   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y - 1
         do k = 2, subgrid_size - 1

            Cv_vec(i) = Cv_vec(i) - b(i)
            i = i + 1

         end do
      end do
   end do

   do iX = 1, macrogrid_size_x - 1
      do iY = 1, macrogrid_size_y
         do k = 2, subgrid_size - 1

            Cv_vec(i) = Cv_vec(i) - b(i)
            i = i + 1

         end do
      end do
   end do

end subroutine Cv

subroutine dot_product(a, b, size, result)
   implicit none
   integer, intent(in) :: size
   real*8, intent(in) :: a(size), b(size)
   real*8, intent(out) :: result
   integer :: i

   result = 0.0d0
   do i = 1, size
      result = result + a(i) * b(i)
   end do

end subroutine dot_product

subroutine conjugate_residuals(subgrid_solver, macrogrid,&
   macrogrid_size_x, macrogrid_size_y,                   &
   subgrid_size, omega,                                  &
   eps_subgrid, eps_interface,                           &
   max_iter_subgrid, max_iter_interface,                 &
   interface_error, time, iter)
   implicit none
   external :: subgrid_solver
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter_subgrid, max_iter_interface
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(in)   :: omega, eps_subgrid, eps_interface
   real*8,  intent(out)   :: interface_error, time
   integer,  intent(out)   :: iter

   integer :: interface_size, i, iX, iY, k
   real(8), dimension((macrogrid_size_x*(macrogrid_size_y-1) + macrogrid_size_y*(macrogrid_size_x-1))*(subgrid_size-2)) :: &
      b, u_new, u_old, r_old, r_new, p_old, p_new, Ar_old, Ar_new, Ap_old, Ap_new

   real(8) :: alpha, beta, dx, dy, temp0, temp1, subgrid_error, old_value, new_val, start_time, end_time

   interface_size = (macrogrid_size_x*(macrogrid_size_y-1) + macrogrid_size_y*(macrogrid_size_x-1))*(subgrid_size-2)
   call compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)

   call cpu_time(start_time)

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
            omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
      end do
   end do
   call s(macrogrid, dx, dy, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size, b)

   u_old = 0.0d0

   call Cv(subgrid_solver, macrogrid, u_old, b, dx, dy, macrogrid_size_x, macrogrid_size_y, &
      subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_old)
   r_new = - b - Ar_old

   p_new = r_new

   call Cv(subgrid_solver, macrogrid, r_new, b, dx, dy, macrogrid_size_x, macrogrid_size_y, &
      subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_new)
   call dot_product(Ar_new, r_new, interface_size, temp0)
   call dot_product(Ar_new, Ar_new, interface_size, temp1)
   alpha = temp0 / temp1

   r_old = r_new - alpha * Ar_new

   call Cv(subgrid_solver, macrogrid, r_old, b, dx, dy, macrogrid_size_x, macrogrid_size_y, &
      subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_old)
   call dot_product(Ar_old, r_old, interface_size, temp0)
   call dot_product(Ar_new, r_new, interface_size, temp1)
   beta = temp0 / temp1

   p_old = r_old + beta * p_new
   Ap_old = Ar_old + beta * Ar_new
   u_old = alpha * p_new

   do iter = 1, max_iter_interface

      call dot_product(Ar_old, r_old, interface_size, temp0)
      call dot_product(Ap_old, Ap_old, interface_size, temp1)
      alpha = temp0 / temp1

      u_new = u_old + alpha * p_old

      interface_error = 0.0d0
      do i = 1, interface_size
         interface_error = interface_error + dabs(u_new(i) - u_old(i))
      end do

      if (interface_error < eps_interface) then
         exit
      end if

      r_new = r_old - alpha * Ap_old

      call Cv(subgrid_solver, macrogrid, r_new, b, dx, dy, macrogrid_size_x, macrogrid_size_y, &
         subgrid_size, interface_size, omega, eps_subgrid, max_iter_subgrid, Ar_new)
      call dot_product(Ar_new, r_new, interface_size, temp0)
      call dot_product(Ar_old, r_old, interface_size, temp1)
      beta = temp0 / temp1

      p_old = r_new + beta * p_old
      Ap_old = Ar_new + beta * Ap_old

      r_old = r_new
      Ar_old = Ar_new
      u_old = u_new

   end do

   call cpu_time(end_time)
   time = end_time - start_time

   call set_interface(macrogrid, u_new, macrogrid_size_x, macrogrid_size_y, subgrid_size, interface_size)
   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
            omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
      end do
   end do

   do iX = 1, macrogrid_size_x - 1
      do iY = 1, macrogrid_size_y - 1

         old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
         new_val = ( &
            macrogrid(iX, iY, subgrid_size-1, subgrid_size) +  &
            macrogrid(iX+1, iY, 2, subgrid_size)  + &
            macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
            macrogrid(iX, iY+1, subgrid_size, 2 )  &
            ) / 4.d0

         macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = new_val
         macrogrid(iX+1, iY,   1,              subgrid_size)   = new_val
         macrogrid(iX,   iY+1, subgrid_size,   1)              = new_val
         macrogrid(iX+1, iY+1, 1,              1)              = new_val
      end do
   end do

end subroutine conjugate_residuals

subroutine simple_iteration(subgrid_solver, macrogrid,   &
   macrogrid_size_x, macrogrid_size_y,                   &
   subgrid_size, omega,                                  &
   eps_subgrid, eps_interface,                           &
   max_iter_subgrid, max_iter_interface,                 &
   interface_error, time, iter)
   implicit none
   external :: subgrid_solver
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter_subgrid, max_iter_interface
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(in)   :: omega, eps_subgrid, eps_interface
   real*8,  intent(out)   :: interface_error, time
   integer,  intent(out)   :: iter

   integer :: iX, iY, k
   real*8  :: dx, dy, old_value, new_val, subgrid_error, start_time, end_time

   call compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)
   call cpu_time(start_time)

   do iter = 1, max_iter_interface

      interface_error = 0.0d0

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
               omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
         end do
      end do

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               old_value = macrogrid(iX, iY, k, subgrid_size)

               new_val = (-macrogrid(iX, iY,   k, subgrid_size-2 ) - macrogrid(iX, iY+1,   k, 3 )  + &
                  4.0d0*macrogrid(iX, iY,   k,   subgrid_size-1 ) +4.0d0*macrogrid(iX, iY+1, k, 2 ))/6.0d0

               macrogrid(iX, iY,   k, subgrid_size) = new_val
               macrogrid(iX, iY+1, k, 1)            = new_val

               interface_error = interface_error + dabs(new_val - old_value)
            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               old_value = macrogrid(iX, iY, subgrid_size, k)

               new_val = (-macrogrid(iX,iY, subgrid_size-2, k) + 4.0d0*macrogrid(iX+1, iY, 2, k ) + &
                  4.0d0*macrogrid(iX,   iY, subgrid_size-1, k) - macrogrid(iX+1, iY, 3, k) &
                  ) / 6.d0

               macrogrid(iX,   iY, subgrid_size, k) = new_val
               macrogrid(iX+1, iY, 1,            k) = new_val

               interface_error = interface_error + dabs(new_val - old_value)
            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y - 1

            old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
            new_val = ( &
               macrogrid(iX, iY, subgrid_size-1, subgrid_size) +  &
               macrogrid(iX+1, iY, 2, subgrid_size)  + &
               macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
               macrogrid(iX, iY+1, subgrid_size, 2 )  &
               ) / 4.d0

            macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = new_val
            macrogrid(iX+1, iY,   1,              subgrid_size)   = new_val
            macrogrid(iX,   iY+1, subgrid_size,   1)              = new_val
            macrogrid(iX+1, iY+1, 1,              1)              = new_val

            interface_error = interface_error + dabs(new_val - old_value)
         end do
      end do

      if (interface_error < eps_interface) then
         exit
      end if

   end do

   call cpu_time(end_time)
   time = end_time - start_time

end subroutine simple_iteration

subroutine simple_iteration_withBoost(subgrid_solver, macrogrid,   &
   macrogrid_size_x, macrogrid_size_y,                   &
   subgrid_size, omega,                                  &
   eps_subgrid, eps_interface,                           &
   max_iter_subgrid, max_iter_interface,                 &
   interface_error, time, iter)
   implicit none
   external :: subgrid_solver
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size, max_iter_subgrid, max_iter_interface
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(in)   :: omega, eps_subgrid, eps_interface
   real*8,  intent(out)   :: interface_error, time
   integer,  intent(out)   :: iter

   integer :: iX, iY, k
   real*8  :: dx, dy, old_value, new_val, subgrid_error, start_time, end_time, eps_Luyst
   real*8  :: p1, p2, norm1, norm2
   logical :: is_start
 
   eps_Luyst = 1.0d-4
   is_start = .true.

   call compute_dxdy(macrogrid_size_x, macrogrid_size_y, subgrid_size, dx, dy)
   call cpu_time(start_time)

   do iter = 1, max_iter_interface

      interface_error = 0.0d0

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y
            call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
               omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
         end do
      end do

      p2 = 0.0d0
      norm2 = 0.0d0

      do iX = 1, macrogrid_size_x
         do iY = 1, macrogrid_size_y - 1
            do k = 2, subgrid_size - 1

               old_value = macrogrid(iX, iY, k, subgrid_size)

               new_val = (-macrogrid(iX, iY,   k, subgrid_size-2 ) - macrogrid(iX, iY+1,   k, 3 )  + &
                  4.0d0*macrogrid(iX, iY,   k,   subgrid_size-1 ) +4.0d0*macrogrid(iX, iY+1, k, 2 ))/6.0d0

               macrogrid(iX, iY,   k, subgrid_size) = new_val
               macrogrid(iX, iY+1, k, 1)            = new_val

               norm2 =norm2 + (new_val - old_value)*(new_val - old_value)
               interface_error = interface_error + dabs(new_val - old_value)
            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y
            do k = 2, subgrid_size - 1

               old_value = macrogrid(iX, iY, subgrid_size, k)

               new_val = (-macrogrid(iX,iY, subgrid_size-2, k) + 4.0d0*macrogrid(iX+1, iY, 2, k ) + &
                  4.0d0*macrogrid(iX,   iY, subgrid_size-1, k) - macrogrid(iX+1, iY, 3, k) &
                  ) / 6.d0

               macrogrid(iX,   iY, subgrid_size, k) = new_val
               macrogrid(iX+1, iY, 1,            k) = new_val

               norm2 =norm2 + (new_val - old_value)*(new_val - old_value)
               interface_error = interface_error + dabs(new_val - old_value)
            end do
         end do
      end do

      do iX = 1, macrogrid_size_x - 1
         do iY = 1, macrogrid_size_y - 1

            old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
            new_val = ( &
               macrogrid(iX, iY, subgrid_size-1, subgrid_size) +  &
               macrogrid(iX+1, iY, 2, subgrid_size)  + &
               macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
               macrogrid(iX, iY+1, subgrid_size, 2 )  &
               ) / 4.d0

            macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = new_val
            macrogrid(iX+1, iY,   1,              subgrid_size)   = new_val
            macrogrid(iX,   iY+1, subgrid_size,   1)              = new_val
            macrogrid(iX+1, iY+1, 1,              1)              = new_val

            norm2 =norm2 + (new_val - old_value)*(new_val - old_value)
            interface_error = interface_error + dabs(new_val - old_value)
         end do
      end do 

      if (is_start .eqv. .true. .and. iter > 1) then
         p2 = sqrt(norm2/norm1)
      end if

      if (is_start .eqv. .true. .and. abs(p2-p1) < eps_Luyst .and. iter > 2 ) then
         interface_error = 0.0d0
         write(*,*) p2
         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y
               call subgrid_solver(macrogrid(iX,iY,:,:), subgrid_size, &
                  omega, eps_subgrid, max_iter_subgrid, subgrid_error, dx, dy)
            end do
         end do

         p2 = 0.0d0
         norm2 = 0.0d0

         do iX = 1, macrogrid_size_x
            do iY = 1, macrogrid_size_y - 1
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, k, subgrid_size)

                  new_val = (-macrogrid(iX, iY,   k, subgrid_size-2 ) - macrogrid(iX, iY+1,   k, 3 )  + &
                     4.0d0*macrogrid(iX, iY,   k,   subgrid_size-1 ) +4.0d0*macrogrid(iX, iY+1, k, 2 ))/6.0d0

                  new_val = (new_val - p2*old_value)/(1-p2)

                  macrogrid(iX, iY,   k, subgrid_size) = new_val
                  macrogrid(iX, iY+1, k, 1)            = new_val

                  norm2 =norm2 + (new_val - old_value)*(new_val - old_value)
                  interface_error = interface_error + dabs(new_val - old_value)
               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y
               do k = 2, subgrid_size - 1

                  old_value = macrogrid(iX, iY, subgrid_size, k)

                  new_val = (-macrogrid(iX,iY, subgrid_size-2, k) + 4.0d0*macrogrid(iX+1, iY, 2, k ) + &
                     4.0d0*macrogrid(iX,   iY, subgrid_size-1, k) - macrogrid(iX+1, iY, 3, k) &
                     ) / 6.d0
                  new_val = (new_val - p2*old_value)/(1-p2)

                  macrogrid(iX,   iY, subgrid_size, k) = new_val
                  macrogrid(iX+1, iY, 1,            k) = new_val

                  norm2 =norm2 + (new_val - old_value)*(new_val - old_value)
                  interface_error = interface_error + dabs(new_val - old_value)
               end do
            end do
         end do

         do iX = 1, macrogrid_size_x - 1
            do iY = 1, macrogrid_size_y - 1

               old_value = macrogrid(iX, iY, subgrid_size, subgrid_size)
               new_val = ( &
                  macrogrid(iX, iY, subgrid_size-1, subgrid_size) +  &
                  macrogrid(iX+1, iY, 2, subgrid_size)  + &
                  macrogrid(iX, iY, subgrid_size, subgrid_size-1) + &
                  macrogrid(iX, iY+1, subgrid_size, 2 )  &
                  ) / 4.d0
               new_val = (new_val - p2*old_value)/(1-p2)

               macrogrid(iX,   iY,   subgrid_size,   subgrid_size)   = new_val
               macrogrid(iX+1, iY,   1,              subgrid_size)   = new_val
               macrogrid(iX,   iY+1, subgrid_size,   1)              = new_val
               macrogrid(iX+1, iY+1, 1,              1)              = new_val

               norm2 =norm2 + (new_val - old_value)*(new_val - old_value)
               interface_error = interface_error + dabs(new_val - old_value)
            end do
         end do 
         is_start = .false.
         
      end if
      norm1 = norm2
      p1 = p2
      

      if (interface_error < eps_interface) then
         exit
      end if

   end do

   call cpu_time(end_time)
   time = end_time - start_time

end subroutine simple_iteration_withBoost

subroutine initialize_constant_boundary(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   integer :: global_size_x, global_size_y
   integer :: iX, iY, i1, j1, lX, lY

   global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               lX = i1 + (iX-1)*subgrid_size - (iX-1)
               lY = j1 + (iY-1)*subgrid_size - (iY-1)

               if (lX == 1 .or. lX == global_size_x .or. &
                  lY == 1 .or. lY == global_size_y) then
                  macrogrid(iX, iY, i1, j1) = 1.0d0
               else
                  macrogrid(iX, iY, i1, j1) = 0.0d0
               end if
            end do
         end do
      end do
   end do
end subroutine initialize_constant_boundary

subroutine initialize_logarithmic_boundary(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size)
   implicit none
   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(inout) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)

   integer :: global_size_x, global_size_y
   integer :: iX, iY, i1, j1, lX, lY
   real*8  :: boundary_val

   global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               lX = i1 + (iX-1)*subgrid_size - (iX-1)
               lY = j1 + (iY-1)*subgrid_size - (iY-1)

               if (lX == 1 .or. lX == global_size_x .or. &
                  lY == 1 .or. lY == global_size_y) then
                  boundary_val = log( sqrt( (x_min + len*(lX-1)/(global_size_x-1))**2 + &
                     (y_min + len*(lY-1)/(global_size_y-1))**2 )*R2/(R1*R1) ) / &
                     log(R2/R1)

                  macrogrid(iX, iY, i1, j1) = boundary_val
               else
                  macrogrid(iX, iY, i1, j1) = 0.0d0
               end if
            end do
         end do
      end do
   end do
end subroutine initialize_logarithmic_boundary

subroutine compute_constant_boundary_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, error)
   implicit none
   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(out) :: error

   integer :: iX, iY, i1, j1
   real*8  :: diff

   error = 0.0d0

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               diff = abs(macrogrid(iX,iY,i1,j1) - 1.0d0)
               if (diff > error) error = diff
            end do
         end do
      end do
   end do

end subroutine compute_constant_boundary_error

subroutine compute_logarithmic_boundary_error(macrogrid, macrogrid_size_x, macrogrid_size_y, subgrid_size, error)
   implicit none
   real*8, parameter :: R1 = 0.1d0, R2 = 1.0d0
   real*8, parameter :: x_min = 0.3d0, y_min = 0.0d0
   real*8, parameter :: len = 0.4d0

   integer, intent(in) :: macrogrid_size_x, macrogrid_size_y, subgrid_size
   real*8,  intent(in) :: macrogrid(macrogrid_size_x,macrogrid_size_y,subgrid_size,subgrid_size)
   real*8,  intent(out) :: error

   integer :: global_size_x, global_size_y
   integer :: iX, iY, i1, j1, lX, lY
   real*8  :: val_exact, diff

   global_size_x = macrogrid_size_x * subgrid_size - (macrogrid_size_x - 1)
   global_size_y = macrogrid_size_y * subgrid_size - (macrogrid_size_y - 1)

   error = 0.0d0

   do iX = 1, macrogrid_size_x
      do iY = 1, macrogrid_size_y
         do i1 = 1, subgrid_size
            do j1 = 1, subgrid_size
               lX = i1 + (iX-1)*subgrid_size - (iX-1)
               lY = j1 + (iY-1)*subgrid_size - (iY-1)

               val_exact = log( sqrt( (x_min + len*(lX-1)/(global_size_x-1))**2 + &
                  (y_min + len*(lY-1)/(global_size_y-1))**2 )*R2/(R1*R1) ) / &
                  log(R2/R1)

               diff = abs(macrogrid(iX,iY,i1,j1) - val_exact)
               if (diff > error) error = diff
            end do
         end do
      end do
   end do

end subroutine compute_logarithmic_boundary_error

subroutine original_sor(u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer :: iter, i, l0, l1, l2, l3
   real*8 :: invdx2, invdy2, f0, f1, u_old

   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

   do iter = 1, max_iter

      error = 0.0d0
      i = u_size + 2

      do l1 = 3, u_size
         do l0 = 3, u_size

            u_old = u(i)
            u(i) = f0*u_old + &
               f1*((u(i-1) + u(i+1))*invdx2 + &
               (u(i-u_size) + u(i+u_size))*invdy2)

            error = error + abs(u(i) - u_old)
            i = i + 1

         end do
         i = i + 2
      end do

      if (error < eps) then
         return
      end if

   end do
end subroutine original_sor

subroutine tiling_sor(u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, parameter :: tile_size = 4

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer :: tile_count, iter, i, l0, l1, l2, l3
   real*8 :: invdx2, invdy2, f0, f1, u_old

   tile_count = (u_size - 2)/tile_size

   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

   do iter = 1, max_iter

      error = 0.0d0
      i = u_size + 2

      do l3 = 1, tile_count
         do l2 = 1, tile_count
            do l1 = 1, tile_size
               do l0 = 1, tile_size

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  error = error + abs(u(i) - u_old)
                  i = i + 1

               end do
               i = i + u_size - tile_size
            end do
            i = i - tile_size*(u_size - 1)
         end do
         i = i + u_size*(tile_size - 1) + 2
      end do

      if (error < eps) then
         return
      end if

   end do
end subroutine tiling_sor

subroutine subtiling_sor(u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, parameter :: tile_size = 4
   integer, parameter :: tile_level = tile_size

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer :: tile_count, iter, i, l0, l1, l2, l3, l4
   real*8 :: invdx2, invdy2, f0, f1, u_old

   tile_count = (u_size - 2)/tile_size

   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

   do iter = 1, max_iter

      error = 0.0d0
      i = u_size + 2

      do l2 = 0, tile_level
         do l1 = 1, tile_size - l2
            do l0 = 1, tile_size - l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
                  error = error + abs(u(i) - u_old)
               end if
               i = i + 1

            end do
            i = i + u_size - l0 + 1
         end do
         i = i - u_size*(l1 - 1)
      end do
      i = i + tile_size

      do l3 = 2, tile_count - 1
         do l2 = 0, tile_level
            do l1 = 1, tile_size - l2
               do l0 = 1, tile_size

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*(l1 - 1) - 1
         end do
         i = i + tile_size + tile_level + 1
      end do

      do l2 = 0, tile_level
         do l1 = 1, tile_size - l2
            do l0 = 1, tile_size + l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
                  error = error + abs(u(i) - u_old)
               end if
               i = i + 1

            end do
            i = i + u_size - l0 + 1
         end do
         i = i - u_size*(l1 - 1) - 1
      end do
      i = i + tile_size + tile_level + u_size*(tile_size - 1) + 3

      do l4 = 2, tile_count - 1
         do l2 = 0, tile_level
            do l1 = 1, tile_size
               do l0 = 1, tile_size - l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1
         end do
         i = i + u_size*(tile_level + 1) + tile_size

         do l3 = 2, tile_count - 1
            do l2 = 0, tile_level
               do l1 = 1, tile_size
                  do l0 = 1, tile_size

                     u_old = u(i)
                     u(i) = f0*u_old + &
                        f1*((u(i-1) + u(i+1))*invdx2 + &
                        (u(i-u_size) + u(i+u_size))*invdy2)

                     if (l2.eq.tile_level) then
                        error = error + abs(u(i) - u_old)
                     end if
                     i = i + 1

                  end do
                  i = i + u_size - l0 + 1
               end do
               i = i - u_size*l1 - 1
            end do
            i = i + u_size*(tile_level + 1) + tile_size + tile_level + 1
         end do

         do l2 = 0, tile_level
            do l1 = 1, tile_size
               do l0 = 1, tile_size + l2

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1 - 1
         end do
         i = i + u_size*(tile_level + 1) + tile_size + tile_level + u_size*(tile_size - 1) + 3
      end do

      do l2 = 0, tile_level
         do l1 = 1, tile_size + l2
            do l0 = 1, tile_size - l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
                  error = error + abs(u(i) - u_old)
               end if
               i = i + 1

            end do
            i = i + u_size - l0 + 1
         end do
         i = i - u_size*l1
      end do
      i = i + u_size*(tile_level + 1) + tile_size

      do l3 = 2, tile_count - 1
         do l2 = 0, tile_level
            do l1 = 1, tile_size + l2
               do l0 = 1, tile_size

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)

                  if (l2.eq.tile_level) then
                     error = error + abs(u(i) - u_old)
                  end if
                  i = i + 1

               end do
               i = i + u_size - l0 + 1
            end do
            i = i - u_size*l1 - 1
         end do
         i = i + u_size*(tile_level + 1) + tile_size + tile_level + 1
      end do

      do l2 = 0, tile_level
         do l1 = 1, tile_size + l2
            do l0 = 1, tile_size + l2

               u_old = u(i)
               u(i) = f0*u_old + &
                  f1*((u(i-1) + u(i+1))*invdx2 + &
                  (u(i-u_size) + u(i+u_size))*invdy2)

               if (l2.eq.tile_level) then
                  error = error + abs(u(i) - u_old)
               end if
               i = i + 1

            end do
            i = i + u_size - l0 + 1
         end do
         i = i - u_size*l1 - 1
      end do

      if (error < eps) then
         return
      end if

   end do
end subroutine subtiling_sor

subroutine subtiling_sor_test_version(u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none

   integer, parameter :: tile_size = 4
   integer, parameter :: tile_level = tile_size

   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error

   integer, dimension(:,:), allocatable :: subtile_indices
   integer :: tile_count, subtile_count, subtile_idx, iter, i, l0, l1, l2, l3, l4, s0, s1, s2, s3
   real*8 :: invdx2, invdy2, f0, f1, u_old

   tile_count = (u_size - 2)/tile_size
   subtile_count = tile_count*tile_count*(tile_level+1)
   allocate(subtile_indices(subtile_count, 4))
   subtile_idx = 1

   do l4 = 1, tile_count

      if (l4.eq.1) then
         s1 = -1
      else if (l4.eq.tile_count) then
         s1 = 1
      else
         s1 = 0
      end if

      if (l4.eq.1) then
         s3 = 0
      else
         s3 = 1
      end if

      do l3 = 1, tile_count

         if (l3.eq.1) then
            s0 = -1
         else if (l3.eq.tile_count) then
            s0 = 1
         else
            s0 = 0
         end if

         if (l3.eq.1) then
            s2 = 0
         else
            s2 = 1
         end if

         do l2 = 0, tile_level
            subtile_indices(subtile_idx, 1) = (l4-1)*tile_size*u_size + (l3-1)*tile_size + u_size + 2 - l2*s3*u_size - l2*s2
            subtile_indices(subtile_idx, 2) = tile_size + l2*s1
            subtile_indices(subtile_idx, 3) = tile_size + l2*s0
            if (l2.eq.tile_level) then
               subtile_indices(subtile_idx, 4) = 1
            else
               subtile_indices(subtile_idx, 4) = 0
            end if
            subtile_idx = subtile_idx + 1

         end do
      end do
   end do

   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)

   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)

   do iter = 1, max_iter

      error = 0.0d0
      i = u_size + 2

      do subtile_idx = 1, subtile_count
         i = subtile_indices(subtile_idx, 1)
         l2 = subtile_indices(subtile_idx, 2)
         l3 = subtile_indices(subtile_idx, 3)
         l4 = subtile_indices(subtile_idx, 4)
         if (l4.eq.1) then
            do l1 = 1, l2
               do l0 = 1, l3

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + abs(u(i) - u_old)
                  i = i + 1

               end do
               i = i + u_size - l3
            end do
         else
            do l1 = 1, l2
               do l0 = 1, l3

                  u_old = u(i)
                  u(i) = f0*u_old + &
                     f1*((u(i-1) + u(i+1))*invdx2 + &
                     (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1

               end do
               i = i + u_size - l3
            end do
         end if
      end do

      if (error < eps) then
         exit
      end if

   end do

   deallocate(subtile_indices)

end subroutine subtiling_sor_test_version

subroutine subtiling_sor_8(u, u_size, factor, eps, max_iter, error, dx, dy)
   implicit none
   integer, parameter :: tile_level = 8
   integer :: tile_size = 8
   integer, intent(in) :: u_size, max_iter
   real*8,  intent(in) :: factor, eps, dx, dy
   real*8,  intent(inout) :: u(u_size*u_size)
   real*8,  intent(out)   :: error
   integer :: tile_count, iter, i, l0, l1, l2, l3, l4
   real*8 :: invdx2, invdy2, f0, f1, u_old
   tile_count = (u_size - 2)/tile_size
   invdx2 = 1.0d0/(dx*dx)
   invdy2 = 1.0d0/(dy*dy)
   f0 = 1.0d0 - factor
   f1 = factor / (2.0d0*invdx2 + 2.0d0*invdy2)
   do iter = 1, max_iter
      error = 0.0d0
      i = u_size + 2
      do l1 = 1, tile_size - 0
         do l0 = 1, tile_size - 0
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 1
         do l0 = 1, tile_size - 1
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 2
         do l0 = 1, tile_size - 2
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 3
         do l0 = 1, tile_size - 3
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 4
         do l0 = 1, tile_size - 4
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 5
         do l0 = 1, tile_size - 5
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 6
         do l0 = 1, tile_size - 6
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 7
         do l0 = 1, tile_size - 7
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1)
      do l1 = 1, tile_size - 8
         do l0 = 1, tile_size - 8
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            error = error + Abs(u_old - u(i))
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i + tile_size - u_size*(l1 - 1)
      do l3 = 2, tile_count - 1
         do l1 = 1, tile_size - 0
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 1
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 2
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 3
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 4
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 5
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 6
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 7
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - u_size*(l1 - 1) - 1
         do l1 = 1, tile_size - 8
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i + tile_level + tile_size - u_size*(l1 - 1)
      end do
      do l1 = 1, tile_size - 0
         do l0 = 1, tile_size + 0
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 1
         do l0 = 1, tile_size + 1
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 2
         do l0 = 1, tile_size + 2
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 3
         do l0 = 1, tile_size + 3
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 4
         do l0 = 1, tile_size + 4
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 5
         do l0 = 1, tile_size + 5
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 6
         do l0 = 1, tile_size + 6
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 7
         do l0 = 1, tile_size + 7
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - u_size*(l1 - 1) - 1
      do l1 = 1, tile_size - 8
         do l0 = 1, tile_size + 8
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            error = error + Abs(u_old - u(i))
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size + tile_level + tile_size*u_size + tile_size + 2
      do l4 = 2, tile_count - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size
         do l1 = 1, tile_size
            do l0 = 1, tile_size - 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + tile_size + u_size*(tile_level + 1)
         do l3 = 2, tile_count - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size - 1
            do l1 = 1, tile_size
               do l0 = 1, tile_size
                  u_old = u(i)
                  u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
                  error = error + Abs(u_old - u(i))
                  i = i + 1
               end do
               i = i - l0 + u_size + 1
            end do
            i = i - l1*u_size + tile_level + tile_size + u_size*(tile_level + 1)
         end do
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 0
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 1
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 2
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 3
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 4
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 5
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 6
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 7
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size
            do l0 = 1, tile_size + 8
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + tile_level*u_size + tile_level + tile_size*u_size + tile_size + 2
      end do
      do l1 = 1, tile_size + 0
         do l0 = 1, tile_size - 0
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 1
         do l0 = 1, tile_size - 1
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 2
         do l0 = 1, tile_size - 2
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 3
         do l0 = 1, tile_size - 3
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 4
         do l0 = 1, tile_size - 4
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 5
         do l0 = 1, tile_size - 5
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 6
         do l0 = 1, tile_size - 6
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 7
         do l0 = 1, tile_size - 7
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size
      do l1 = 1, tile_size + 8
         do l0 = 1, tile_size - 8
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            error = error + Abs(u_old - u(i))
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size + tile_size + u_size*(tile_level + 1)
      do l3 = 2, tile_count - 1
         do l1 = 1, tile_size + 0
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 1
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 2
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 3
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 4
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 5
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 6
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 7
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size - 1
         do l1 = 1, tile_size + 8
            do l0 = 1, tile_size
               u_old = u(i)
               u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
               error = error + Abs(u_old - u(i))
               i = i + 1
            end do
            i = i - l0 + u_size + 1
         end do
         i = i - l1*u_size + tile_level + tile_size + u_size*(tile_level + 1)
      end do
      do l1 = 1, tile_size + 0
         do l0 = 1, tile_size + 0
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 1
         do l0 = 1, tile_size + 1
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 2
         do l0 = 1, tile_size + 2
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 3
         do l0 = 1, tile_size + 3
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 4
         do l0 = 1, tile_size + 4
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 5
         do l0 = 1, tile_size + 5
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 6
         do l0 = 1, tile_size + 6
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 7
         do l0 = 1, tile_size + 7
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      do l1 = 1, tile_size + 8
         do l0 = 1, tile_size + 8
            u_old = u(i)
            u(i) = f0*u_old + f1*((u(i-1) + u(i+1))*invdx2 + (u(i-u_size) + u(i+u_size))*invdy2)
            error = error + Abs(u_old - u(i))
            i = i + 1
         end do
         i = i - l0 + u_size + 1
      end do
      i = i - l1*u_size - 1
      if (error < eps) then
         return
      end if
   end do
end subroutine subtiling_sor_8

