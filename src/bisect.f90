! Setup
! Read Parameters, e.g S, N_x etc
! Declare Initial Guess and guess step e.g b1 = 5, s0 = 1
! Create Matrix With Corresponding Initial Condtions
! Evolve initial state using PI.f90 step_forward until critical value reached (either laminar or turbulent)
!  if laminar then b2 = b1+s, lam1 = True(Laminar) 
!  else b2 = b1-s, lam1 = False(Turbulent)
! Evolve b2 state using PI.f90 step_forward until critical value reached
!  if laminar then lam2 = True
!  else lam2 = False
!  if lam2 != lam1 then continue
!  else repeat

! Repeats loop until one of b1 and b2 is turbulent, the other is not.
! Enter new loop - re-order b1 and b2 so b2 is the greater of the two
! b_new = b1+b2/2.0
! Evolve with init amplitude b_new
! If state turbulent then b2=b_new
! else b1 = b_new
include 'PI.f90'
module bisection_mod
  use, intrinsic :: iso_c_binding
  use PI_mod
  implicit none
  integer :: limit = 55
  real(dp) :: a0, a1, a2, a, step, lam, turb, t_min
  real(dp), dimension(:), allocatable :: amp
  complex(C_DOUBLE_COMPLEX), dimension(:,:), allocatable :: ival
  integer :: stat, bisect_count, j_t, j_f, mode, n
  namelist /bisinfo/ mode, limit, a0, step, lam, turb, t_min
contains 

  subroutine readBisectData
    open(21, file = 'Bisect/Input/bisect.in')
    read(21, nml = bisinfo)
    close(21)
  end subroutine readBisectData

  subroutine initial_state(a)
    implicit none
    real(dp) :: a
    rval(:,:) = 0.0_dp
    do n = 1, N_x
       rval(n,1:4) = 0.0_dp
    end do
    do n = 1, N_x
       rval(n,3:4) = a*exp(real(-(n-N_x/2+1)**2)/16)
    end do
  end subroutine initial_state

  subroutine new_state(a)
    implicit none
    real(dp) :: a, f
    f = a/a0
    fval(:,1:4) = ival(:,:)*f
  end subroutine new_state

  function amplitude(state)
    complex(dp), dimension(N_x,4) :: state
    real(dp) :: amplitude
    amplitude = 0.0_dp
    do n = 1, N_x
       amplitude = amplitude + state(n,3)*conjg(state(n,3))*N_x*dx
       amplitude = amplitude + state(n,4)*conjg(state(n,4))*N_x*dx
       amplitude = amplitude + state(n,2)*conjg(state(n,2))*N_x*dx
    end do
    amplitude = amplitude
  end function amplitude

  function ramplitude(state)
    complex(dp), dimension(N_x,4) :: state
    real(dp) :: ramplitude
    ramplitude = 0.0_dp
    do n = 1, N_x
       ramplitude = ramplitude + state(n,3)*conjg(state(n,3))*dx
       ramplitude = ramplitude + state(n,4)*conjg(state(n,4))*dx
       ramplitude = ramplitude + state(n,2)*conjg(state(n,2))*dx
    end do
    ramplitude = ramplitude
  end function ramplitude

end module bisection_mod

program bisection
  use PI_mod
  use bisection_mod
  implicit none
  call PI_setup(0)  ! Reads in parameters and initial conds for PI_mod
  call readBisectData
  allocate(ival(N_x, 4))          
  allocate(amp(nframes+1))     
!  mode = 1
!  a0 = 0.3_dp
  a1 = a0
!  t_min = 20
!  step = 0.1_dp
  stat = 0 ! No Bounds found
!  lam = 5d-1
!  turb = 100.0_dp
  write(*,*) dx, turb, lam
  call initial_state(a1)
  write(*,*) 'Initialised state'
  call init_time
  write(*,*) 'Init time'
  call fft(1)
  write(*,*) 'FFT'
  ival(:,:) = fval(:,1:4)
  open(31, file = 'Bisect/Output/bis.dat')

  ! if (mode == 2) then
  !    write(*,*) N_x/2, N_x/2+N_x/20
  !    write(*,*) amplitude(fval(:,1:4))
  !    write(*,*) ramplitude(rval(:, 1:4))
  !    stop
  ! end if

  if (mode == 2) then
     do j_t = 1, N_t
        call time_step(dt)
        write(31,'(f30.15)') amplitude(fval(:,1:4))
     end do
     write(*,*) njump
     stop
  end if
          
! ************************************ !
! TRY INITIAL GUESS                    !
! ************************************ !
  write(*,*) 'Main loop'
  do
     call time_step(dt)
     write(*,*) t, t_total, amplitude(fval(:,1:4))
     if (t_total > t_min) then
!           call backwardfft
!           call unshift
        if (amplitude(fval(:,1:4)) <= lam) then
           write(*,*) 'Lower bound found: ', a1, amplitude(fval(:,1:4))
           stat = 1 ! Lower bound found
           a = a1
           exit
        else if(amplitude(fval(:,1:4)) >= turb) then
           write(*,*) 'Upper bound found: ', a1, amplitude(fval(:,1:4))
           a2 = a1
           a = a1
           step = -step
           stat = 2 ! Upper bound found
           exit
        end if
     else
        cycle
     end if
  end do

! ************************************!
! FIND ANOTHER BOUND                  !
  do
     a = a + step
     write(*,*) 'Trying'
     write(*,*)  a
     call new_state(a)
     call init_time
     do
        call time_step(dt)
        if (t_total > t_min) then
           if (amplitude(fval(:,1:4)) <= lam) then
              write(*,*) 'Laminar',t_total, amplitude(fval(:,1:4))
              exit
           else if (amplitude(fval(:,1:4)) >= turb) then
              write(*,*) 'Turbulent',t_total, amplitude(fval(:,1:4))
              exit
           end if
        end if
     end do
        
     if (amplitude(fval(:,1:4)) <= lam) then
        a1 = a
        if (stat == 1) then
           cycle
        else
           stat = 3
           exit
        end if
     end if

     if (amplitude(fval(:,1:4)) >= turb) then
        a2 = a
        if (stat == 2) then
           cycle
        else
           stat = 3
           exit
        end if
     end if
  end do
  

  write(*,*) 'Lower limit: ', a1
  write(*,*) 'Upper limit: ', a2
  open(99, file = 'Bisect/Output/log.txt')

  ! Begin bisections !
  do bisect_count = 1, limit
     a = (a1 + a2)*0.5_dp                  ! Bisect the current upper and lower bounds
     write(*,*) 'Trying'
     write(*,'(f30.25)')  a
     write(99,*) 'Trying'
     write(99,'(f30.25)') a
     call new_state(a)                     ! Create a new initial state with amplitude a
     call init_time                        ! Initialise the time system in PI_mod
     write(*,*) bisect_count
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Loops until either a laminar or turbulent state is reached !
! waiting at least t_min time to check - to ensure low time  !
! transient growth is ignored.                               !

     amp(:) = 0.0_dp
     j_t = 0
     j_f = 1
     amp(1) = amplitude(fval(:,1:4))
     do
        j_t = j_t+1
        if ((mod(j_t, steps_per_frame) == 0) .and. (j_f < nframes + 1)) then
           j_f = j_f+1
           amp(j_f) = amplitude(fval(:,1:4))
        else if (j_f >= nframes + 1) then
!           write(*,*) 'j_f > nframes - consider using large N_t'   
        end if

        call time_step(dt)                     ! Takes a single step dt of the PI model  
        if (t_total > t_min) then
           if (amplitude(fval(:,1:4)) <= lam) then  ! If laminar, update lower bound
              write(*,*) 'Laminar',t_total, amplitude(fval(:,1:4)), amp(j_f)
              write(99,*) 'Laminar',t_total, amplitude(fval(:,1:4)), amp(j_f)
              a1 = a
              exit
           else if (amplitude(fval(:,1:4)) >= turb) then ! IF turbulent, update upper bound
              write(*,*) 'Turbulent',t_total, amplitude(fval(:,1:4)), amp(j_f)
              write(99,*) 'Turbulent',t_total, amplitude(fval(:,1:4)), amp(j_f)
              a2 = a
              exit
           end if
        end if
     end do
     write(31,'(10000f30.15)') amp
  end do

  if (mode == 1) then ! Calculate final trajectory
     j_f = t_total
     write(*,*) 'Storing final solution'
     call new_state(a)
     call init_time
     call fft(-1)
     call writeout(0)
     do j_t = 1, N_t
        call time_step(dt)
        if (mod(j_t, steps_per_frame) == 0) then
           call write_frame(0)
           !       write(31,'(f21.15)') amplitude(fval(:,1:4))
        end if
     end do

     open(11, file = 'Bisect/Output/out_turb.dat',      form = 'formatted')
     call new_state(a*1.01)
     call init_time
     call fft(-1)
     do m = 1, 4
        write(11,1) real(rval(:,m))
        write(11,1) DIMAG(rval(:,m))
     end do
     write(*,*) steps_per_frame

     do j_t = 1, N_t
        call time_step(dt)
        if(mod(j_t, steps_per_frame) == 0) then
           call fft(-1)
           do m = 1, 4
              write(11,1) real(rval(:,m))
              write(11,1) DIMAG(rval(:,m))
           end do
        end if
     end do
     close(11)
     stop
1    format(100000f30.15)

  end if

 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
 ! INTERPOLATION REGION ~~~~~~~~~~~~~~~ !

  a2 = a + 0.00000000000000003
  a1 = a - 0.00000000000000003
  write(*,*) a, a1, a2

 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
 ! UPPER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
  amp(:) = 0.0_dp
  j_t = 0
  j_f = 1
  call new_state(a2)
  call init_time
  amp(1) = amplitude(fval(:,1:4))
  do 
     j_t = j_t+1
     if (mod(j_t, steps_per_frame) == 0) then
        j_f = j_f+1
        amp(j_f) = amplitude(fval(:,1:4))
     end if
     call time_step(dt)                     ! Takes a single step dt of the PI model  
     if (t_total > t_min) then
        if (amplitude(fval(:,1:4)) <= lam) then  ! If laminar, update lower bound
           write(*,*) 'Laminar',t_total, amplitude(fval(:,1:4))
!           write(99,*) 'Laminar',t_total, amplitude(fval(:,1:4))
           exit
        else if (amplitude(fval(:,1:4)) >= turb) then ! IF turbulent, update upper bound
           write(*,*) 'Turbulent',t_total, amplitude(fval(:,1:4))
!           write(99,*) 'Turbulent',t_total, amplitude(fval(:,1:4))
           exit
        end if
     end if
  end do
  write(31, '(10000f25.15)') amp

  amp(:) = 0.0_dp
  j_t = 0
  j_f = 1
  call new_state(a1)
  call init_time
  amp(1) = amplitude(fval(:,1:4))
  do 
     j_t = j_t+1
     if (mod(j_t, steps_per_frame) == 0) then
        j_f = j_f+1
        amp(j_f) = amplitude(fval(:,1:4))
     end if
     call time_step(dt)                     ! Takes a single step dt of the PI model  
     if (t_total > t_min) then
        if (amplitude(fval(:,1:4)) <= lam) then  ! If laminar, update lower bound
           write(*,*) 'Laminar',t_total, amplitude(fval(:,1:4))
!           write(99,*) 'Laminar',t_total, amplitude(fval(:,1:4))
           exit
        else if (amplitude(fval(:,1:4)) >= turb) then ! IF turbulent, update upper bound
           write(*,*) 'Turbulent',t_total, amplitude(fval(:,1:4))
!          write(99,*) 'Turbulent',t_total, amplitude(fval(:,1:4))
           exit
        end if
     end if
  end do
  write(31, '(10000f21.15)') amp

  close(99)
contains
  
end program bisection




