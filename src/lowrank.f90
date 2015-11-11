module lowrankmod

implicit none

contains

	subroutine lowrank(rank, epsabs, epsrel)
		double precision, optional :: epsrel, epsabs
		integer, optional :: rank
		integer :: m, n, lda, ldu, ldvt, lwork, info
		character :: jobu, jobvt
<<<<<<< HEAD
		double precision, dimension(:), allocatable :: S
		double precision, dimension(:), allocatable :: work, U
		double precision, dimension(1:1) :: dummy
=======
		double precision, dimension(min(size(A%Ut,1),size(A%Ut,2))) :: S
		double precision, allocatable :: lwork, U
>>>>>>> parent of d2beb5e... working full and matprod
		type(Matrix), pointer :: A,B

		call matrixReader(A)
<<<<<<< HEAD
		write(*,*) 'matrix readed'
		jobu = 'O' !first min(m,n) left singular vectors returned in U
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
		call dgesvd(jobu, jobvt, m, n, A%Ut, m, S, dummy, 1, B%Vt, n, work, lwork, info)

		write(*,*) 'info: ', info
		write(*,*) A%Ut
		write(*,*) 'B:'
		write(*,*) B%Vt(:,1)
		write(*,*) B%Vt(:,2)
		write(*,*) B%Vt(:,3)
		write(*,*) B%Vt(:,4)
		write(*,*) B%Vt(:,5)

		deallocate(work)
		deallocate(A%Ut, A%Vt, B%Ut, B%Vt)
		deallocate(A,B)
=======
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
>>>>>>> parent of d2beb5e... working full and matprod
		
	end subroutine

end module

