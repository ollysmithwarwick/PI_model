module initialise
  implicit none

  integer, parameter :: dp = selected_real_kind(15,307)
  real(dp), parameter :: pi_ = 4*atan(1.0_dp)

  contains

  subroutine init(n)
    implicit none
    integer :: n, k
    real(dp), dimension(n) :: r1, r2
    real(dp),    dimension(n) :: in1, in2
    complex(dp), dimension(n) :: n1, phi1, n0, E
    complex, parameter :: i = (0.0_dp, 1.0_dp)
    integer :: boundL, boundR, snapshot
    real(dp) :: noise = 0.0
    integer, dimension(1) :: seed=1234
    real(dp) :: factor = 1.0
    namelist /frameinfo/ boundL, boundR, snapshot, factor, noise, seed

    n0(:)   = 0.0_dp
    n1(:)   = 0.0_dp
    E(:)    = 0.0_dp
    phi1(:) = 0.0_dp
    open(21, file = 'POPI/Input/frameinfo.in')
    read(21, nml = frameinfo)
    close(21)

!    call random_init
    write(*,*) 'Noise =', noise
    call random_seed(put = seed)
    write(*,*) 'Noise =', noise

    do k = 1, n
!       n1(k)   = 0.5*cos(pi*2.0_dp*real(k-1)/n)
!       n1(k)   = 0.1*exp(-0.25_dp*real(k-n/2)*real(k-n/2))
!       phi1(k) = n1(k)
!       n0(k) = 0.5*cos(pi*2.0_dp*real(k-1)/n)
!       phi1(k) = 0.02
    end do

    open(20, file = './POPI/Input/POPI_in.dat')
    do k = 1, snapshot
       call random_number(r1(:))
       read(20, '(10000f30.15)') in1(:)
       read(20, '(10000f30.15)') in2(:)
       n0(boundL:boundR) = in1(boundL:boundR) * factor 
       n0(:) = n0(:) + noise * (r1(:) - 0.5)
       write(*,*) noise
       call random_number(r1(:))
       read(20, '(10000f30.15)') in1(:)
       read(20, '(10000f30.15)') in2(:) 
       E(boundL:boundR) = in1(boundL:boundR) * factor
       E(:) = E(:) + noise * (r1(:) - 0.5)
       
       call random_number(r1(:))
       call random_number(r2(:))
       read(20, '(10000f30.15)') in1(:)
       read(20, '(10000f30.15)') in2(:)
       phi1(boundL:boundR) = factor * (in1(boundL:boundR) + i * in2(boundL:boundR)) 
       phi1(:) = phi1(:) + noise * (r1(:)-0.5) +&
            i * (noise * (r2(:) - 0.5))
       read(20, '(10000f30.15)') in1(:)
       read(20, '(10000f30.15)') in2(:)
       n1(boundL:boundR) = factor * (in1(boundL:boundR) + i * in2(boundL:boundR))
       n1(:) = n1(:) + noise * (r1(:)-0.5) +&
            i * (noise * (r2(:)-0.5))
    end do

    open(1, file = 'Input/init.dat')

    write(1,2) real(n0(:))
    write(1,2) dimag(n0(:))
    write(1,2) real(E(:))
    write(1,2) dimag(E(:))
    write(1,2) real(phi1(:))
    write(1,2) dimag(phi1(:))
    write(1,2) real(n1(:))
    write(1,2) dimag(n1(:))

    close(1)

2   format (100000f30.15)
  end subroutine init
end module initialise
