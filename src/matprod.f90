module matprodmod

use matrixconverter

implicit none

contains

	subroutine matprod()
		integer :: n
		type(Matrix), pointer :: A, B
		
		call matrixReader(A)
		call matrixReader(B)

		!write(*,*) 'in matprod, matrices readed in'
		if (A%full) then
			!write(*,*) 'going to leftfullprod()'
			call leftfullprod(A, B)
		else
			if (B%full) then
				call rankfullprod(A,B)
			else
				call rankrankprod(A,B)
			endif
		endif

		deallocate(A)
		deallocate(B)
	end subroutine
	
	subroutine leftfullprod(A, B)
		! A*B = C - met A van volle rang dus
		! aangezien we de getransponeerden opslaan: (C^t = C%Vt^t*C%Ut)
		! met B niet van volle rank: C^t = (A*B)^t = B^t*A^t = (B%Ut^t*B%Vt)^t*A%Ut = B%Vt^t*B%Ut*A%Ut
		! met B van volle rank: C^t = (A*B)^t = B^t*A^t = B%Ut*A%Ut
		! in beide gevallen moet enkel het product B%Ut*A%Ut berekend worden
		type(Matrix), pointer :: A, B, C
		double precision :: alpha, beta
		integer :: m,n,k
		
		m = size(B%Ut,1)
		n = size(A%Ut,2)
		k = size(B%Ut,2)
		alpha = 1.0
		beta = 0.0
		allocate(C)
		C%full = B%full
		allocate(C%Ut(m,n))
		allocate(C%Vt(size(B%Vt,1),size(B%Vt,2)))
		C%Vt = B%Vt
		
		!write(*,*) 'DEBUGGING INFORMATION'
		!write(*,*) 'size of A%Ut: ', size(A%Ut,1), ',', size(A%Ut,2)
		!write(*,*) 'size of B%Ut: ', size(B%Ut,1), ',', size(B%Ut,2)
		
		call dgemm('N','N', m, n, k, alpha, B%Ut, m, A%Ut, k, beta, C%Ut, m)

		call matrixWriter(C)
		
		deallocate(C)

	end subroutine
	
	subroutine rankfullprod(A, B)
		! A*B = C - met A niet van volle rank en B wel
		! aangezien we de getransponeerden opslaan: (C^t = C%Vt^t*C%Ut)
		! C^t = (A*B)^t = B^t*A^t = B%Ut*(A%Ut^t*A%Vt)^t = B%Ut*A%Vt^t*A%Ut
		! we berekenen dus C%Ut = A%Ut en C%Vt = B%Ut*A%Vt^t

		type(Matrix), pointer :: A, B, C
		double precision :: alpha, beta
		integer :: m,n,k
		
		m = size(B%Ut,1)
		n = size(A%Vt,1)
		k = size(B%Ut,2)
		alpha = 1.0
		beta = 0.0

		allocate(C)
		C%full = .false.
		allocate(C%Ut(size(A%Ut,1),size(A%Ut,2)))
		C%Ut = A%Ut
		allocate(C%Vt(m,n))
		
		!write(*,*) 'DEBUGGING INFORMATION'
		!write(*,*) 'size of A%Ut: ', size(A%Ut,1), ',', size(A%Ut,2)
		!write(*,*) 'size of B%Ut: ', size(B%Ut,1), ',', size(B%Ut,2)
		
		call dgemm('N','T', m, n, k, alpha, B%Ut, m, A%Vt, n, beta, C%Ut, m)

		call matrixWriter(C)
		
		deallocate(C)

	end subroutine
	
	subroutine rankrankprod(A, B)
		! A*B = C - mat A en B allebei niet van volle rang
		! aangezien we de getransponeerden opslaan: (C^t = C%Vt^t*C%Ut)
		! C^t = (A*B)^t = B^t*A^t = (B%Ut^t*B%Vt)^t*(A%Ut^t*A%Vt)^t
		!						  = B%Vt^t*B%Ut*A%Vt^t*A%Ut
		! we berekenen dus C%Ut = (B%Ut*A%Vt^t)*A%Ut en C%Vt = B%Vt

		type(Matrix), pointer :: A, B, C
		double precision, dimension(size(B%Ut,1),size(A%Vt,1)) :: tempM
		double precision :: alpha, beta
		integer :: m,n,k
		
		m = size(B%Ut,1)
		n = size(A%Vt,1)
		n2 = size(A%Ut,2)
		k = size(B%Ut,2)
		alpha = 1.0
		beta = 0.0

		allocate(C)
		C%full = .false.
		allocate(C%Ut(m,n))
		allocate(C%Vt(size(A%Ut,1),size(A%Ut,2)))
		C%Vt = B%Vt
		
		!write(*,*) 'DEBUGGING INFORMATION'
		!write(*,*) 'size of A%Ut: ', size(A%Ut,1), ',', size(A%Ut,2)
		!write(*,*) 'size of B%Ut: ', size(B%Ut,1), ',', size(B%Ut,2)
		
		call dgemm('N','T', m, n, k, alpha, B%Ut, m, A%Vt, n, beta, tempM, m)
		call dgemm('N','N', m, n2, n, alpha, C, m, A%Ut, n, beta, C%Ut, m)

		call matrixWriter(C)
		
		deallocate(C)

	end subroutine

end module

