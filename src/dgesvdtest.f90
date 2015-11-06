program main

implicit none

contains

subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests DGESVD.
!
!  Discussion:
!
!    DGESVD computes the singular value decomposition:
!
!      A = U * S * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 6
  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ), parameter :: lwork = 3*min(m,n) + max ( max(m,n), 2*min(m,n) )

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldvt
  character jobu
  character jobvt
  real ( kind = 8 ) s(min(m,n))
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sigma(m,n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) vt(n,n)
  real ( kind = 8 ) work(lwork)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For a double precision real matrix (D)'
  write ( *, '(a)' ) '  in general storage mode (GE):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DGESVD computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V'''
!
!  Set A.
!
  seed = 123456789

  call r8mat_uniform_01 ( m, n, seed, a )

  call r8mat_print ( m, n, a, '  The matrix A:' )
!
!  Compute the singular values and singular vectors.
!
  jobu = 'A'
  jobvt = 'A'
  lda = m
  ldu = m
  ldvt = n

  call dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
    lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DGESVD returned nonzero INFO = ', info
    return
  end if

  call r8vec_print ( min ( m, n ), s, '  Singular values' )

  call r8mat_print ( m, m, u, '  Left singular vectors U:' )
  call r8mat_print ( n, n, vt, '  Right singular vectors V'':' )

  sigma(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    sigma(i,i) = s(i)
  end do

  a(1:m,1:n) = matmul ( u(1:m,1:m), matmul ( sigma(1:m,1:n), vt(1:n,1:n) ) )

  call r8mat_print ( m, n, a, '  The product U * S * V'':' )

  return
end subroutine

end program
