!*************************************************************************
!  Example with one extra parameter T.
!  Can put parameter and constraint =0 when not required.
!- - - - - - - - - - - - - - -
!  Newton vector:
!    x(1:3)   = (/S, shift, phase/)
!    x(4: 6* N_x+ 3) = vector containing real and imag parts (for tilde quantities) of each quantity at each point
!
!  Extra constraints:
!    (F(x)-x). dx/dt = 0 .
!       no update along direction of trajectory.
!    
!    (F_(shift)(x) - x).dx/dX (X being the coordinate, x is the vector)
!       no update in the X direction (the solver is solving for a length translation anyway)
!
!    (F_(phase)(x) - x).dx/d(phi)
!       no update in a phase rotatin direction (solver solves for global phase transition)
!
!
!  Jacobian approximation:
!    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps Used to quickly approximate dF/dx, made more accurate by reducing eps
!    
!
!  File guess.in and guess.out:
!     S		initial guess for shear
!     shift     initial guess for translation length
!     phase     initial guess for phase shift
!     n_period  number of boundary condition periods (1 bc period = 2*pi/LS)  
!     ndts	number of timesteps taken in one period - if set to -1 this will split the total time into 0.01 unit chunks,
!               though as the total time changes with each guess (since S changes -> T = n_period*2*pi/L*S changes) the size
!               of each step changes - need to keep an eye to make sure it doesn't get too big
!     info_out  Prints out the info integer in guess.out and is ignored in guess.in (so you can copy the out file to the in file without changing anything)
!               0 If Converged
!               1 If the code is running - this should never be printed out, since it should only print once it is complete
!               2 If nits > nits_max
!               3 If trust region gets too small
!               4 If there was a memory error
!     nits      Number of newton iterations used to reach convergence (or whatever other outcome)
!               An approximate measure of how far the guess was from the solution, but depends on many factors
!*************************************************************************
 module orbit
! This module contains all the variables relevant to the NewtonHook method
   implicit none
   save
					!* Params.in:
   integer          :: mgmres, nits 	!m for gmres(m), max newton its
   double precision :: rel_err		!relative error
   double precision :: del, mndl, mxdl  !delta for hookstep-trust-region
   double precision :: gtol, epsJ	!gmres tolerance, eps in eqn above
					!* Typical params.in:
					!    100  100
					!    5d-10
					!   -1d0  1d-7  1d+7
					!    1d-3  1d-6

   integer :: n          !size of state + parameters
   double precision   :: tol
   integer :: info, ndts, n_int, file_stat = 0

   double precision :: dt__=0.001
 end module orbit

!*************************************************************************
include 'NewtonHook.f90' ! Newton-Krylov-Hookstep Algorithm
include 'GMRESm.f90'     ! GMRESm Algorithm
include 'PI.f90'         ! PI code
include 'init_mod.f90'   ! This is a module form of init.f90 - so it can be called repeatedly in the script rather than having to do it separately
 PROGRAM MAIN
!*************************************************************************
   use orbit
   use newton
   use PI_mod
   use initialise
   implicit none
   external :: getrhs, multJ, multJp, saveorbit
   double precision, external :: dotprod
   integer, external :: bestshift
   real, external :: bestphase
   double precision :: d
   integer :: i1, i2, i3
   double complex, dimension(:), allocatable :: final_f
!   character(4) :: cnum
   write(*,*) 'newton_PI: Setup'
   call PI_setup(0)
   write(*,*) 'newton_PI: Calling init'
   call init(N_x)
   write(*,*) 'newton_PI: Calling readinit'
   call readinit
   n = 6 * N_x + 3 ! Total number of variables, 3 global variables (S, shift, phase), 6*N_x system variables

   write(*,*) 'newton_PI: Reading params.in'

   open(99,status='old',file='POPI/Input/params.in')
   read(99,*) mgmres, nits
   read(99,*) rel_err
   read(99,*) del, mndl, mxdl
   read(99,*) gtol, epsJ
   close(99)

   write(*,*) 'newton_PI: Allocating memory'

   allocate(new_x(n))
   allocate(new_fx(n))
   allocate(final_f(N_x))
   allocate(pos(N_x, 4, 2))
   new_x = 0d0

   write(*,*) 'newton_PI: Memory allocated. Reading guess.in'

   open(99,status='old',file='POPI/Input/guess.in')
   read(99,*) new_x(1) ! S
   read(99,*) new_x(2) ! Spatial shift
   read(99,*) new_x(3) ! Phase shift
   read(99,*) n_int
!   new_x(3) = 0d0
   read(99,*) ndts
   close(99)

!   write(*,*) 'guess.in:', new_x(1), new_x(2), new_x(3), n_int ! new_x(1) = S
!   write(*,*) 'p       :', n_int, '*2*pi/(', L, '*', S, ') = ', n_int*2*pi*L_inv/S       
   
   if(ndts==-1) ndts = nint(real(n_int*2*pi*L_inv/new_x(1))/0.01d0) ! If ndts = -1 set it to be the number required to split into ~0.01 unit chunks
   ! Originally 0.001d0 but 0.01 seems to be accurate enough?

!   write(*,*) 'ndts    :', ndts

   steps_per_frame = int(ndts/nframes)  ! E.g if one has 1234 time steps but want 100 frames steps_per_frame becomes 12
   nframes = int(ndts/steps_per_frame)  ! With 12 steps per frame however we can actually squeeze more frames from it so nframes is updated to 102
   S = new_x(1)                         ! Set S to its initial guess new_x(1)

   call writeplotinfo                   ! Writes data to plotinfo.in for plotting later

!*****
!   LOAD STATE TO  new_x(3:) - this is your initial guess. The data structure can be determined from the allocation here, it's inverse is simple too
!*****
   write(*,*) 'newton_PI: Determining pos'
   call makepos(N_x)

   write(*,*) 'newton_PI: Writing new_x'
   do i1 = 1, N_x
      new_x(pos(i1,     NBAR, RL)) =  dble(rval(i1,     NBAR))
      new_x(pos(i1,        E, RL)) =  dble(rval(i1,        E))
      new_x(pos(i1, PHITILDE, RL)) =  dble(rval(i1, PHITILDE))
      new_x(pos(i1, PHITILDE, IM)) = dimag(rval(i1, PHITILDE))
      new_x(pos(i1,   NTILDE, RL)) =  dble(rval(i1,   NTILDE))
      new_x(pos(i1,   NTILDE, IM)) = dimag(rval(i1,   NTILDE))
   end do

   ! scale params by norm of state - i.e all errors are relative to the INITIAL GUESS not the final state

   d = dotprod(-1,new_x,new_x)
   tol  = rel_err * dsqrt(d)
!   tol  = rel_err
   del  = del     * dsqrt(d)
   mndl = mndl    * dsqrt(d)
   mxdl = mxdl    * dsqrt(d)

   info = 1

   call store_traj(n_int*2*pi*L_inv/S, ndts) ! Stores the initial guess' trajectory in out.dat - if an error occurs or the code is exited early this can be helpful

   if (new_x(2) == 0) then
      call getlshift
      call store_traj(n_int*2*pi*L_inv/S, ndts)
   end if

   if (new_x(3) == 0) then
      call getphase
      call store_traj(n_int*2*pi*L_inv/S, ndts)
   end if

   write(*,*) 'newton_PI: Calling newtonhook'
   call newtonhook(getrhs, multJ, multJp, saveorbit, dotprod, &
                   mgmres, n, gtol, tol, del, mndl, mxdl, nits, info) ! This is the main solving subroutine - treat it as a black box
                                                                      ! Solves f(x) = 0 by inverting df(x_n)/dx . s = f(x_n) to get s and taking a step in (roughly) this direction
                                                                      !    getrhs : calculates y = f(x_n) for a given x i.e given a state, Shear, shift etc will calculate the state at the later time
                                                                      !             y(1:3) = 0.0 for the RHS of the 3 extra conditions
                                                                      !    multJ  : Finds df(x_n)/dx . s. Uses Jacobian approximation. Also defines the LHS of each condition in components 1:3
                                                                      !    multJp : Preconditioner for multJ - not used.
                                                                      !    saveorbit : Does everything you want to do after each guess is tried
                                                                      !    dotprod: Finds the dot product of two vectors (n specifies whether to include the global variables or not)
                                                                      !    mgmres : dimension of mgmres/how many gmres iterations to do per newton iteration - I just set this to be the same as n
                                                                      !         n : Number of variables i.e 6*N_x + 3
                                                                      !    gtol   : Target tolerance for gmres (typically 1d-3)
                                                                      !    tol    : Target tolerance below which the amplitude of the residual will be considered to have converged.
                                                                      !    del    : Initial size of trust region (see numerical methods book by Dennis and Schabel, borrow from Olly Smith if you wish)
                                                                      !    mndl,mxdl : min and max allowed values of del - lower bound often triggers so need to work out if this is a good thing or not
                                                                      !    nits   : Max number of newton iterations to be performed
                                                                      !    info   : An output value that gives the status of the output. 0 is converged soln

   write(*,*) 'newton_PI: storing final solution in final.dat'

   open(98, status = 'unknown', file = 'POPI/Output/final.dat')

   do i1 = 1, n
      write(98, '(f26.18)') new_x(i1)
   end do
   close(98)

   write(*,*) 'newton_PI: storing final trajectory in out_PI.dat'
   
   write(*,*) 'NDTS = ', ndts
   write(*,*) '   S = ', new_x(1)
   call store_traj(n_int*2*pi*L_inv/S, ndts)
   write(*,*) 'newton_PI: Done'

   stop

 contains

!*************************************************************************
 END PROGRAM MAIN
!*************************************************************************

!-------------------------------------------------------------------------
!  determine best guess for Lshift
!-------------------------------------------------------------------------
 subroutine getlshift
   use PI_mod
   use newton
   implicit none
   integer, external :: bestshift

   new_x(2) = bestshift(N_x, sqrt(real(rval(:, NTILDE)*conjg(rval(:,NTILDE)))), sqrt(new_x(pos(:, NTILDE, RL))*&
        new_x(pos(:, NTILDE, RL)) + new_x(pos(:, NTILDE, IM))*new_x(pos(:, NTILDE, IM))))*dx

 end subroutine getlshift

 subroutine getphase
   use PI_mod
   use newton
   implicit none
   real, external :: bestphase

   new_x(3) = bestphase(N_x, new_x(pos(:, NTILDE, RL)) + i*new_x(pos(:, NTILDE, IM)), rval(:, NTILDE))

 end subroutine getphase

!-------------------------------------------------------------------------
!  function to be minimised
!-------------------------------------------------------------------------
 subroutine getrhs(n_,x, y)
   use orbit
   use PI_mod
   implicit none
   integer,          intent(in)  :: n_
   double precision, intent(in)  :: x(n)
   double precision, intent(out) :: y(n)
   double precision :: x_(n), y_(n)

   x_ = x
   call steporbit(ndts, x_(2), x_(3), x_, y_)
   y = y_ - x

   y(1) = 0d0					! constraints, rhs=0
   y(2) = 0d0					! constraints, rhs=0
   y(3) = 0d0                                   ! phase constraint, rhs = 0
         
 end subroutine getrhs


!-------------------------------------------------------------------------
!  Jacobian of function + lhs of constraints on update
!-------------------------------------------------------------------------
 subroutine multJ(n_,x, y)
   use newton
   use orbit,    only : n, epsJ, dt__
   use PI_mod
   implicit none
   integer,          intent(in) :: n_
   double precision, intent(in)   :: x(n)
   double precision, intent(out)  :: y(n)
   double precision, external   :: dotprod
   integer :: j
   double precision :: eps  
   double precision :: ss(n)
    				! (F(x0+eps.x)-F(x0))/eps
   eps = dsqrt(dotprod(1,x,x))
   if(eps==0d0)  stop 'multJ: eps=0 (1)'
   eps = epsJ * dsqrt(dotprod(1,new_x,new_x)) / eps
   if(eps==0d0)  stop 'multJ: eps=0 (2)'
   y = new_x + eps*x
   call getrhs(n_,y, ss)
   y = (ss - new_fx) / eps
      				! contstraint,
				! no update in trajectory direction
   call steporbit(1, 0d0, 0d0, new_x, ss)
   ss = (ss - new_x) / dt__
   y(1) = dotprod(-1,ss,x)

   ! Constraint 2 No translation in update
   call init_time
   do j = 1, N_x
      rval(j,     NBAR) = new_x(pos(j,    NBAR, RL))
      rval(j,        E) = new_x(pos(j,       E, RL))
      rval(j, PHITILDE) = new_x(pos(j, PHITILDE, RL))  + i * new_x(pos(j, PHITILDE, IM))
      rval(j,   NTILDE) = new_x(pos(j,   NTILDE, RL))  + i * new_x(pos(j,   NTILDE, IM))
   end do

   call fft(1)

   do j = 1, N_x
      fval(j,     NBAR) = -k(j)*fval(j,     NBAR)
      fval(j,        E) = -k(j)*fval(j,        E)
      fval(j, PHITILDE) = -k(j)*fval(j, PHITILDE)
      fval(j,   NTILDE) = -k(j)*fval(j,   NTILDE)
   end do
   
   call fft(-1)
   
   do j = 1, N_x
      ss(pos(j,     NBAR, RL)) =  dble(rval(j,     NBAR))
      ss(pos(j,        E, RL)) =  dble(rval(j,        E))
      ss(pos(j, PHITILDE, RL)) =  dble(rval(j, PHITILDE))
      ss(pos(j, PHITILDE, IM)) = dimag(rval(j, PHITILDE))
      ss(pos(j,   NTILDE, RL)) =  dble(rval(j,   NTILDE))
      ss(pos(j,   NTILDE, IM)) = dimag(rval(j,   NTILDE))      
   end do

   y(2) = dotprod(-1, ss, x)

   ! Constraint 3, No global phase shift

   do j = 1, N_x
      ss(pos(j,     NBAR, RL)) =   0.0_dp
      ss(pos(j,        E, RL)) =   0.0_dp
      ss(pos(j, PHITILDE, IM)) =   dble(rval(j, PHITILDE))
      ss(pos(j, PHITILDE, RL)) = -dimag(rval(j, PHITILDE))
      ss(pos(j,   NTILDE, IM)) =   dble(rval(j,   NTILDE))
      ss(pos(j,   NTILDE, RL)) = -dimag(rval(j,   NTILDE))      
   end do
   y(3) = dotprod(-1, ss, x)

 end subroutine multJ


!-------------------------------------------------------------------------
!  preconditioner for multJ.  Empty - no preconditioner required
!-------------------------------------------------------------------------
 subroutine multJp(n, x)
   implicit none
   integer,          intent(in)    :: n
   double precision, intent(inout) :: x(n)
 end subroutine multJp


!-------------------------------------------------------------------------
!  called at each newton iteration
!-------------------------------------------------------------------------
 subroutine saveorbit()
   use newton
   use orbit
   use PI_mod
!   use io
   implicit none
   integer :: j1, j2, j
   double precision :: norm_x 
   double precision, external :: dotprod

   do j1 = 1, N_x
      rval(j1,     NBAR) = new_x(pos(j1,    NBAR, RL))
      rval(j1,        E) = new_x(pos(j1,       E, RL))
      rval(j1, PHITILDE) = new_x(pos(j1, PHITILDE, RL))  + i * new_x(pos(j1, PHITILDE, IM))
      rval(j1,   NTILDE) = new_x(pos(j1,   NTILDE, RL))  + i * new_x(pos(j1,   NTILDE, IM))
   end do
   S = new_x(1)
   norm_x = dsqrt(dotprod(-1,new_x,new_x))
   call fft(-1)

   open(99,status='unknown',access='append',file='POPI/Output/newton.dat')
   if(new_nits==0)  write(99,*) ndts, mgmres, n
   write(99,'(2I6,4e13.5)')  &
         new_nits, new_gits, new_tol, new_del, new_tol/norm_x, norm_x
   close(99)
!  newton its completed
!  num gmres its for last newton it
!  current error
!  current trust region
!  relative error
!  size of state
   if (file_stat == 0) then
      open(99,status='unknown',file='POPI/Output/guesses.dat')
   else
      open(99,status='unknown',access='append',file='POPI/Output/guesses.dat')
   end if
   if(new_nits==0)  write(99,*) ndts
   write(99,'(A16, 1I6, 3e26.18)') 'nits, S, L, phi:', new_nits, new_x(1), new_x(2), new_x(3)
   write(99,'(A16, 1I6, 2e26.18)') '   n_int, T, dt:', n_int, n_int*2*pi*L_inv/S, n_int*2*pi*L_inv/(S*ndts)
   close(99)
!  newton its completed
!  current guess for T

!*****
!   SAVE MOST RECENT STATE  new_x(3:)
!*****
   if (file_stat == 0) then
      open(99,status='unknown',file='POPI/Output/states.dat')
      file_stat = 1
   else
      open(99,status='unknown', access = 'append', file='POPI/Output/states.dat')
   end if
   do j1 = 1, 2
      write(99, '(20000e26.18)') real(rval(:, j1))
   end do
   
   do j1 = 3, 4
      write(99, '(20000e26.18)')  real(rval(:, j1))
      write(99, '(20000e26.18)') dimag(rval(:, j1))
   end do
!   write(99,'(20000e26.18)') real(rval(:))
   close(99)

 end subroutine saveorbit

 subroutine store_traj(p, ndts_)
   use PI_mod
   use newton
   use orbit
   implicit none
   integer :: i1, ndts_
   double precision, external :: dotprod
   double precision :: dt_, p, y(n)
   write(*,*) 'Store_traj'
   open(11, file = 'POPI/Output/initFinal.dat')
   do i1 = 1, N_x
      write(11,2) real(rval(i1, NBAR)), real(rval(i1, E)), rval(i1, PHITILDE), rval(i1, NTILDE)
   end do
   close(11)

   do i1 = 1, N_x
      rval(i1,     NBAR) = new_x(pos(i1,    NBAR, RL))
      rval(i1,        E) = new_x(pos(i1,       E, RL))
      rval(i1, PHITILDE) = new_x(pos(i1, PHITILDE, RL))  + i * new_x(pos(i1, PHITILDE, IM))
      rval(i1,   NTILDE) = new_x(pos(i1,   NTILDE, RL))  + i * new_x(pos(i1,   NTILDE, IM))
      rval(i1,   B_PLUS) = rval(i1, PHITILDE) + i * rval(i1, NTILDE)
      rval(i1,  B_MINUS) = rval(i1, PHITILDE) - i * rval(i1, NTILDE)
   end do

   open(11, file = 'POPI/Output/guess.out')
   write(11, *) new_x(1)
   write(11, *) new_x(2)
   write(11, *) new_x(3)
   write(11, *) n_int
   write(11, *) ndts
   write(11, *) info
   write(11, *) nits
   close(11)
   open(11, file = 'POPI/Output/final_traj.dat')
   S = new_x(1)
   dt_ = p/ndts_
   steps_per_frame = int(ndts_/nframes) ! MAY NEED TO BE REVISED
!   call readinit ! Temp here for testing
   call init_time
   call writeout(0)
! Set rval
   call fft(1)
   write(*,*) 'ndts = ', ndts_, nframes, steps_per_frame
   do i1 = 1, ndts_
      call time_step(dt_)
      if (mod(i1, steps_per_frame) == 0) then
         call write_frame(0)
      end if
   end do
   call write_frame(0)
   write(*,*) 'Final state'
   do i1 = 1, N_x
      fval(i1,     NBAR) = fval(i1,     NBAR) * exp(i*(k(i1)*new_x(2))) !+ new_x(3)))
      fval(i1,        E) = fval(i1,        E) * exp(i*(k(i1)*new_x(2))) !+ new_x(3)))
      fval(i1, PHITILDE) = fval(i1, PHITILDE) * exp(i*(k(i1)*new_x(2) + new_x(3)))
      fval(i1,   NTILDE) = fval(i1,   NTILDE) * exp(i*(k(i1)*new_x(2) + new_x(3)))
   end do
   
   call fft(-1)
   write(*,*) '--------------------', rval(85, NBAR)
   call write_frame(0)
   
   do i1 = 1, N_x
      new_fx(pos(i1,     NBAR, RL)) =  dble(rval(i1,     NBAR))
      new_fx(pos(i1,        E, RL)) =  dble(rval(i1,        E))
      new_fx(pos(i1, PHITILDE, RL)) =  dble(rval(i1, PHITILDE))
      new_fx(pos(i1, PHITILDE, IM)) = dimag(rval(i1, PHITILDE))
      new_fx(pos(i1,   NTILDE, RL)) =  dble(rval(i1,   NTILDE))
      new_fx(pos(i1,   NTILDE, IM)) = dimag(rval(i1,   NTILDE))
   end do

   call steporbit(ndts_, new_x(2), new_x(3), new_x, y)
   
   write(*,*) 'store_traj:', dsqrt(dotprod(-1, new_fx, new_fx)), dsqrt(dotprod(-1, y, y))
   

   new_fx = new_fx - new_x
   write(*,*) 'store_traj: ', dsqrt(dotprod(-1, new_fx, new_fx)), dsqrt(dotprod(-1, y-new_x, y -new_x))
   
   dt = dt_
   N_t = ndts
   call writeplotinfo

   close(11)
1  format(100000f30.15)
2  format(6F23.18)
 end subroutine store_traj
   
!-------------------------------------------------------------------------
! dot product.  can flag to exclude parameter T.  Could include weights
!-------------------------------------------------------------------------
 double precision function dotprod(n_,a,b)
   use orbit
   implicit none
   integer,          intent(in) :: n_
   double precision, intent(in) ::  a(n),  b(n)
   double precision :: a1(n), b1(n)
!   double precision :: d,d_
   integer :: n1
   n1 = 1
   a1 = a
   b1 = b
   a1(1:3) = 1.0*a(1:3)
   b1(1:3) = 1.0*b(1:3)
   if(n_==-1) n1 = 4
   dotprod = dot_product(a1(n1:n),b1(n1:n))
 end function dotprod


!-------------------------------------------------------------------------
!  timestep
!-------------------------------------------------------------------------
 subroutine steporbit(ndts_, jump_, phase_, x, y)
   use orbit
   use newton
   use PI_mod
!   use io
   implicit none
   integer,          intent(in)  :: ndts_
   double precision, intent(in)  :: x(n)
   double precision, intent(in)  :: jump_, phase_
   double precision, intent(out) :: y(n)
   double precision, save :: dt_

!   double precision, dimension(1:n-1) :: state,k1,k2,k3,k4
!   double precision :: theta,phi,thedot,phidot
   integer :: j1

   S = x(1)
   if(ndts_/=1) then
      dt_ = n_int*2*pi*L_inv /(S * dble(ndts_)) ! n_int is number of intervals
                                                ! 2pi/LS is the time in one interval
   end if
   call init_time
   
   do j1 = 1, N_x
      rval(j1,     NBAR) = x(pos(j1,    NBAR, RL))
      rval(j1,        E) = x(pos(j1,       E, RL))
      rval(j1, PHITILDE) = x(pos(j1, PHITILDE, RL))  + i * x(pos(j1, PHITILDE, IM))
      rval(j1,   NTILDE) = x(pos(j1,   NTILDE, RL))  + i * x(pos(j1,   NTILDE, IM))
   end do

   call fft(1)
   
   do j1=1,ndts_
      call time_step(dt_)
      call fft(-1)
   enddo

   do j1 = 1, N_x
      fval(j1,     NBAR) = fval(j1,     NBAR) * exp(i*(k(j1)*jump_)) !+ phase_))
      fval(j1,        E) = fval(j1,        E) * exp(i*(k(j1)*jump_)) !+ phase_))
      fval(j1, PHITILDE) = fval(j1, PHITILDE) * exp(i*(k(j1)*jump_ + phase_))
      fval(j1,   NTILDE) = fval(j1,   NTILDE) * exp(i*(k(j1)*jump_ + phase_))
   end do

   call fft(-1)

   do j1 = 1, N_x
      y(pos(j1,     NBAR, RL)) =  dble(rval(j1,     NBAR))
      y(pos(j1,        E, RL)) =  dble(rval(j1,        E))
      y(pos(j1, PHITILDE, RL)) =  dble(rval(j1, PHITILDE))
      y(pos(j1, PHITILDE, IM)) = dimag(rval(j1, PHITILDE))
      y(pos(j1,   NTILDE, RL)) =  dble(rval(j1,   NTILDE))
      y(pos(j1,   NTILDE, IM)) = dimag(rval(j1,   NTILDE))      
   end do

 contains

 end subroutine steporbit

 function bestshift(n_, a, b)
   !Finds the best value of length shift
   implicit none
   integer, intent(in) :: n_
   double precision,    intent(in) :: a(n_), b(n_) ! Initial guess and final position
   integer :: bestshift
   integer             :: i1, i2
   double precision    :: prod(n_), s
   write(*,*) 'bestshift called'
   write(*,*) a
!   write(*,*) b
   prod(:) = 0.0
   do i1 = 1, n_
      do i2 = 1, n_
         !      prod(i1) = prod(i1) + (a(mod(i1 - 1 + i2, n_)) - b(i1 - 1))**2.0
         s = (a(mod(i2+i1-2, n_)+1) - b(i2))
         prod(i1) = prod(i1) + s * s
!         write(*,*) prod(1), a(mod(i2-1, n_) +1)
      end do
   end do
   write(*,*) 'prod = ', prod(:)
   write(*,*) 'minloc(prod) = ', minloc(prod(:)) - 1

!   bestshift = minloc(prod) - 1  
   bestshift = minloc(prod(:), dim = 1)
   bestshift = mod((bestshift - 1 + n_/2), n_) - n_/2
   write(*,*) 'The bestshift integer is: ', bestshift
 end function bestshift

 function bestphase(n_, a, b)
   !Finds the best value of phase shift
   implicit none
   integer, intent(in) :: n_
   double complex, intent(in) :: a(n_), b(n_) ! Initial guess and final position
   real :: bestphase
   integer             :: i1, i2
   double precision    :: prod(629)
   double complex      :: s
   double precision, parameter  :: pi = 4.0*atan(1.0)
   double complex,   parameter  :: i = (0.0, 1.0)
   write(*,*) 'bestphase called'
   prod(:) = 0.0
   do i1 = 0, 628
      do i2 = 1, n_
         !      prod(i1) = prod(i1) + (a(mod(i1 - 1 + i2, n_)) - b(i1 - 1))**2.0
         s = a(i2)*zexp(-i*i1*0.01) - b(i2)
         prod(i1 + 1) = prod(i1 + 1) + s * conjg(s)
!         write(*,*) prod(1), a(mod(i2-1, n_) +1)
      end do
   end do
   write(*,*) 'prod = ', prod(:)
   write(*,*) 'minloc(prod) = ', minloc(prod(:))

!   bestshift = minloc(prod) - 1 
   bestphase = minloc(prod(:) - 1 , dim = 1) * 0.01
   write(*,*) 'bestphase = ', bestphase
 end function bestphase
