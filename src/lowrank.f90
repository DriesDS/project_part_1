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
		type(Matrix), pointer :: A,B

		write(*,*) 'in lowrank'
		call matrixReader(A)
		write(*,*) 'matrix readed'
		jobu = 'S' !first min(m,n) left singular vectors returned in U
		jobvt = 'S' !first min(m,n) right singular vectors returned in Vt
		m = size(A%Ut,1)
		n = size(A%Ut,2)
		write(*,*) 'variables set, allocating S'
		allocate(S(min(size(A%Ut,1),size(A%Ut,2))))
		lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))
		allocate(work(lwork))
		info = 0
		
		write(*,*) 'allocating B'
		allocate(B)
		B%full = .false.
		allocate(B%Ut(m,m))
		allocate(B%Vt(n,n))
		
		write(*,*) 'calling dgesvd'
		write(*,*) 'jobu: ', jobu, 'jobvt: ', jobvt, 'm: ', m, 'n: ', n
		write(*,*) 
		call dgesvd(jobu, jobvt, m, n, A%Ut, m, S, B%Ut, m, B%Vt, n, work, lwork, info)

		write(*,*) 'info: ', info
		write(*,*) B%Ut
		write(*,*) B%Vt

		deallocate(work)
		deallocate(A%Ut, A%Vt, B%Ut, B%Vt)
		deallocate(A,B)
		
	end subroutine

end module

