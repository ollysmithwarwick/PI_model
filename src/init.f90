program initial 
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,307)
  integer :: N_x = 80
  real(dp), parameter :: pi = 4*atan(1.0_dp)
  real(dp) :: dx
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
  complex(dp), dimension(:), allocatable :: n1, phi1
  real(dp),    dimension(:), allocatable :: n0, phi0
  integer :: k


  allocate(n0(N_x))
  allocate(phi0(N_x))
  allocate(n1(N_x))
  allocate(phi1(N_x))

  dx = 1.0_dp / N_x
  n1(:) = 0.0_dp
  phi0(:) = 0.0_dp
  n0(:) = 0.0_dp
  phi1(:) = 0.0_dp

  ! Change this loop to change initial state. This is just a square initial condition
  do k = N_x/2-1, N_x/2+1
     n1(k) = 0.25_dp
     phi1(k) = 0.25_dp
  end do

  open(1, file = 'init.dat')

  do k = 0, N_x-1
     write(1, 2) n0(k+1), phi0(k+1), phi1(k+1), n1(k+1)
  end do
  close(1)
2 format (6F23.18)
end program initial
