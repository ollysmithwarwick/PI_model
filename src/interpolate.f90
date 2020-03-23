program interpolate
  implicit none
  
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: n = 400
  real(dp), dimension(n) :: in1, in2
  REAL(DP) :: d
  complex(dp), dimension(n) :: n0, E, ntilde, phitilde
  complex(dp), parameter :: i = (0.0_dp, 1.0_dp)

  read(*,*) d

  write(*,*) 'read done'
  open(20, file = './out.dat')
  read(20, 1) in1(:)
  read(20, 1) in2(:)
  n0(:) = in1(:) + i * in2(:)
  read(20, 1) in1(:)
  read(20, 1) in2(:)
  E(:) = in1(:) + i * in2(:)
  read(20,1) in1(:)
  read(20, 1) in2(:)
  phitilde(:) = in1(:) + i * in2(:)
  read(20, 1) in1(:)
  read(20, 1) in2(:)
  ntilde(:) = in1(:) + i * in2(:)
  close(20)

  open(20, file = './out_old.dat')
  read(20, 1) in1(:)
  read(20, 1) in2(:)
  n0(:) = n0(:) + d * (n0(:) - (in1(:) + i * in2(:)))
  read(20, 1) in1(:)
  read(20, 1) in2(:)
  E(:) = E(:) + d*(E(:) - (in1(:) + i * in2(:)))
  read(20,1) in1(:)
  read(20, 1) in2(:)
  phitilde(:) = phitilde(:) + d*(phitilde(:) - (in1(:) + i * in2(:)))
  read(20, 1) in1(:)
  read(20, 1) in2(:)
  ntilde(:) = ntilde(:) + d*(ntilde(:) - (in1(:) + i * in2(:)))
  close(20)

  open(21, file = 'POPI_in.dat')
  write(21, 1) real(n0(:))
  write(21, 1) dimag(n0(:))
  write(21, 1) real(E(:))
  write(21, 1) dimag(E(:))
  write(21, 1) real(phitilde(:))
  write(21, 1) dimag(phitilde(:))
  write(21, 1) real(ntilde(:))
  write(21, 1) dimag(ntilde(:))
  write(21, 1) real(n0(:))
  write(21, 1) dimag(n0(:))
  write(21, 1) real(n0(:))
  write(21, 1) dimag(n0(:))
 
 close(20) 

 status = system('cp out.dat out_old.dat')

1 format(100000f30.15)
end program interpolate
