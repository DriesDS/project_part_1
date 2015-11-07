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
		endif

		deallocate(A)
		deallocate(B)
	end subroutine
	
	subroutine leftfullprod(A, B)
		! A*B = C - met A van volle rang dus
		! aangezien we de getransponeerden opslaan:
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
		C%full = .false.
		allocate(C%Ut(m,n))
		
		!write(*,*) 'DEBUGGING INFORMATION'
		!write(*,*) 'size of A%Ut: ', size(A%Ut,1), ',', size(A%Ut,2)
		!write(*,*) 'size of B%Ut: ', size(B%Ut,1), ',', size(B%Ut,2)
		
		call dgemm('N','N', m, n, k, alpha, B%Ut, m, A%Ut, k, beta, C%Ut, m)

		if (B%full) then
			call matrixWriter(C)
		else
			call matrixWriter(C,B)
		endif

		deallocate(C)

	end subroutine
	
	subroutine rankfullprod()
	end subroutine
	
	subroutine rankrankprod()
	end subroutine

end module

