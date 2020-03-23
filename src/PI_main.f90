include 'PI.f90'
program PI_program
  ! --------------------------------------------------------------------------!
! A program for solving the plasma interchange model using spectral methods !
! When compiling one must have FFTW3 loaded, include fftw3.f03 and link the fftw3 library !

! DECLARE VARIABLES !
  use PI_mod
  implicit none
  integer :: j
  write(*,*) 'pre_setup'
  call PI_setup(1)
  write(*,*) 'post_setup'
  call writeout(0)
  call fft(1)
!  do j = 1, N_x
!     do n = 1, 6
!        if (((real(fval(j, n)) <= 10E-15) .and. (dimag(fval(j,n)) <=10E-15))) then
!           fval(j,n) = 0.0_dp
!        end if
!     end do
!  end do
  call init_time
  j = 0
  call writeout(0)
  ! Iterate over timesteps
  do j = 1, N_t
     call time_step(dt)
     if (mod(j,steps_per_frame)== 0) then
        call write_frame(1)
     end if
  end do
  call write_frame(1)

  write(*,*) 'Done with ',njump,' jumps'
  close(11)
end program PI_program
