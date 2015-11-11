module lowrankmod

use matrixconverter

implicit none

contains

	subroutine lowrank(rank, epsabs, epsrel)
		double precision, optional :: epsrel, epsabs
		integer, optional :: rank
		integer :: m, n, lda, ldu, ldvt, lwork, info
		character :: jobu, jobvt
		double precision, dimension(:), allocatable :: S
		double precision, dimension(:), allocatable :: work, U
		double precision, dimension(1:1) :: dummy
		type(Matrix), pointer :: A,B

		call matrixReader(A)
		m = size(A%Ut,1)
		n = size(A%Ut,2)
		allocate(S(min(size(A%Ut,1),size(A%Ut,2))))
		lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))
		allocate(work(lwork))
		info = 0
		
		allocate(B)
		B%full = .false.
		if (m<n) then
			jobu = 'O' !first min(m,n) left singular vectors returned in A
			jobvt = 'S' !first min(m,n) right singular vectors returned in Vt
			allocate(B%Vt(n,n))
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, S, dummy, 1, B%Vt, n, work, lwork, info)
		else
			jobu = 'S'
			jobvt = 'O'
			allocate(B%Ut(m,m))
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, S, B%Ut, m, dummy, 1, work, lwork, info)
		endif
		
		deallocate(work)
		deallocate(A,B)
		
	end subroutine

end module

