module changegrid

  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  integer, parameter :: dp = selected_real_kind(15,307)

  integer :: n_old, n_new
  type(C_PTR) :: planf, planb, new_planb, new_planf
  complex(C_DOUBLE_COMPLEX), dimension(:, :), allocatable :: rval, new_rval
  complex(C_DOUBLE_COMPLEX), dimension(:),    allocatable :: in, out, new_in, new_out

contains
  
  subroutine setup(n_old, n_new)
    implicit none
    integer, intent(in) :: n_old
    integer, intent(in) :: n_new
    allocate(in(n_old))
    allocate(out(n_old))
    allocate(new_in(n_new))
    allocate(new_out(n_new))
    allocate(rval(n_old, 6))
    allocate(new_rval(n_new,6))

    planf       = fftw_plan_dft_1d(n_old, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    planb       = fftw_plan_dft_1d(n_old, out, in, FFTW_BACKWARD, FFTW_ESTIMATE)
    new_planf   = fftw_plan_dft_1d(n_new, new_in, new_out, FFTW_FORWARD, FFTW_ESTIMATE)
    new_planb   = fftw_plan_dft_1d(n_new, new_out, new_in, FFTW_BACKWARD, FFTW_ESTIMATE)

  end subroutine setup


  subroutine readIn()

    integer :: m
    complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
    real(dp), dimension(n_old) :: in1, in2
    open(20, file = './POPI_in.dat')

    do m = 1, 6
       read(20, '(10000f30.15)') in1(:)
       read(20, '(10000f30.15)') in2(:)
       rval(:, m) = in1(:) + i * in2(:)
    end do

    close(20)
  end subroutine readIn

  subroutine change_grid()
    implicit none
    
    integer :: m
    complex(C_DOUBLE_COMPLEX), dimension(n_old, 6) :: fval
    complex(C_DOUBLE_COMPLEX), dimension(n_new, 6) :: new_fval
    
    new_fval = 0.0_dp
    do m = 1, 6
       in(:) = rval(:, m)
       call fftw_execute_dft(planf, in, out)
       fval(:, m) = out(:)/sqrt(real(n_old))
       new_fval(1:n_old/2+1, m) = fval(1:n_old/2+1, m)
       new_fval(n_new-n_old/2+2:n_new, m) = fval(n_old/2+2:n_old, m)
       new_out(:) = new_fval(:,m)
       call fftw_execute_dft(new_planb, new_out, new_in)
       new_rval(:,m) = new_in(:)/sqrt(real(n_old))
    end do
  end subroutine change_grid

  subroutine writeout
    integer :: m
    open(11, file = 'POPI_in_newgrid.dat', form = 'formatted')
    do m = 1, 6
       write(11, 1) real(new_rval(:,m))
       write(11, 1) DIMAG(new_rval(:,m))
    end do
    close(11)
1 format(100000f30.15)
  end subroutine writeout

end module

program change
  use changegrid
  read(*,*) n_old
  n_new = n_old + 2
  call setup(n_old, n_new)
  call readIn()
  call change_grid()
  call writeout()
  write(*,*) 'Done'
end program change
