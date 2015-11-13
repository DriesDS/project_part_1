module lowrankmod

use matrixconverter

implicit none

contains

	! we berekenen A^t = U,S,Vt
	! we weten dus dat A = Vt^t*S*U^t
	! We moeten dus Vt opslaan in A%Ut en U^t in A%Vt
	subroutine lowrank(rank, epsabs, epsrel)
		double precision, optional :: epsrel, epsabs
		integer, optional :: rank
		integer :: m, n, lda, ldu, ldvt, lwork, info
		character :: jobu, jobvt
		double precision, dimension(:), allocatable :: S
		double precision, dimension(:), allocatable :: work, U
		double precision, dimension(1,1) :: dummy
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
			jobu = 'S' !first min(m,n) left singular vectors returned in A
			jobvt = 'S' !first min(m,n) right singular vectors returned in Vt
			allocate(B%Vt(m,m))
			allocate(B%Ut(m,n))
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, S, B%Vt, m, B%Ut, m, work, lwork, info)
			B%Vt = transpose(B%Vt)
		else
			jobu = 'O' !use O for storing in A
			jobvt = 'S'
			allocate(B%Ut(n,n))
			allocate(B%Vt(n,m))
			write(*,*) size(B%Vt,1), n
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, S, A%Ut, m, B%Ut, n, work, lwork, info)
			write(*,*) size(B%Vt,1), size(B%Vt,2)
			B%Vt = transpose(A%Ut)
			write(*,*) size(B%Vt,1), size(B%Vt,2)
		endif
		
		call matrixWriter(B)
		
		deallocate(work)
		deallocate(A,B)
		
	end subroutine

end module

