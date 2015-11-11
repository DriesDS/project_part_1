module lowrankmod

implicit none

contains

	subroutine lowrank(rank, epsabs, epsrel)
		double precision, optional :: epsrel, epsabs
		integer, optional :: rank
		integer :: m, n, lda, ldu, ldvt, lwork, info
		character :: jobu, jobvt
		double precision, dimension(min(size(A%Ut,1),size(A%Ut,2))) :: S
		double precision, allocatable :: lwork, U
		type(Matrix), pointer :: A,B

		call matrixReader(A)
		jobu = 'O' !first min(m,n) left singular vectors returned in A
		jobvt = 'S' !first min(m,n) right singular vectors returned in Vt
		m = size(A%Ut,1)
		n = size(A%Ut,2)
		lwork = 0

		allocate(A%Vt(n,n))

		call dgesvd(jobu, jobvt, m, n, A%Ut, m, S, U, 0, A%Vt, n, 0, lwork, info)

		write(*,*) A%Ut
		write(*,*) A%Vt

		deallocate(A)
		
	end subroutine

end module

