module matprodmod

use matrixconverter

implicit none

contains

	subroutine matprod(A, B, C)
		type(Matrix), pointer, intent(in) :: A, B
		type(Matrix), pointer, intent(inout) :: C
		
		if (A%full) then
			call leftfullprod(A, B, C)
		else
			if (B%full) then
				call rankfullprod(A, B, C)
			else
				call rankrankprod(A, B, C)
			endif
		endif

	end subroutine
	
	subroutine leftfullprod(A, B, C)
		! A*B = C - met A van volle rang dus
		! aangezien we de getransponeerden opslaan: (C^t = C%Vt^t*C%Ut)
		! met B niet van volle rank: C^t = (A*B)^t = B^t*A^t = (B%Ut^t*B%Vt)^t*A%Ut = B%Vt^t*B%Ut*A%Ut
		! met B van volle rank: C^t = (A*B)^t = B^t*A^t = B%Ut*A%Ut
		! in beide gevallen moet enkel het product B%Ut*A%Ut berekend worden
		type(Matrix), pointer, intent(in) :: A, B
		type(Matrix), pointer, intent(inout) :: C
		double precision :: alpha, beta
		integer :: m,n,k
		
		m = size(B%Ut,1)
		n = size(A%Ut,2)
		k = size(B%Ut,2)
		alpha = 1.0
		beta = 0.0

		allocate(C)
		C%full = B%full
		C%pointU = .false.
		allocate(C%Ut(m,n))

		if (.not. C%full) then
			C%pointV = .true.
			C%Vt => B%Vt
		endif
		
		call dgemm('N','N', m, n, k, alpha, B%Ut, m, A%Ut, k, beta, C%Ut, m)

	end subroutine
	
	subroutine rankfullprod(A, B, C)
		! A*B = C - met A niet van volle rank en B wel
		! aangezien we de getransponeerden opslaan: (C^t = C%Vt^t*C%Ut)
		! C^t = (A*B)^t = B^t*A^t = B%Ut*(A%Ut^t*A%Vt)^t = B%Ut*A%Vt^t*A%Ut
		! we berekenen dus C%Ut = A%Ut en C%Vt = A%Vt*B%Ut^t

		type(Matrix), pointer, intent(in) :: A, B
		type(Matrix), pointer, intent(inout) :: C
		double precision :: alpha, beta
		integer :: m,n,k
		
		m = size(A%Vt,1)
		n = size(B%Ut,1)
		k = size(B%Ut,2)
		alpha = 1.0
		beta = 0.0

		allocate(C)
		C%full = .false.
		C%pointU = .true.
		C%pointV = .false.
		allocate(C%Vt(m,n))
		C%Ut => A%Ut
				
		call dgemm('N','T', m, n, k, alpha, A%Vt, m, B%Ut, n, beta, C%Vt, m)

	end subroutine
	
	subroutine rankrankprod(A, B, C)
		! A*B = C - mat A en B allebei niet van volle rang
		! aangezien we de getransponeerden opslaan: (C^t = C%Vt^t*C%Ut)
		! C^t = (A*B)^t = B^t*A^t = (B%Ut^t*B%Vt)^t*(A%Ut^t*A%Vt)^t
		!                         = B%Vt^t*B%Ut*A%Vt^t*A%Ut
		! Als matrix B van laagste rang is:
		! we berekenen C%Ut = (B%Ut*A%Vt^t)*A%Ut en C%Vt = B%Vt
		! Als matrix A van laagste rang is:
		! we berekenen C%Ut = A%Ut en C%Vt = A%Vt*B%Ut^t*B%Vt
		type(Matrix), pointer, intent(in) :: A, B
		type(Matrix), pointer, intent(inout) :: C
		double precision, dimension(size(B%Ut,1),size(A%Vt,1)) :: tempM
		double precision :: alpha, beta
		integer :: m, n, n2, k
		logical :: Asmallest

		Asmallest = size(B%Ut,1) > size(A%Ut,1)

		m = size(B%Ut,1)
		n = size(A%Vt,1)
		k = size(B%Ut,2)
		alpha = 1.0
		beta = 0.0
		call dgemm('N','T', m, n, k, alpha, B%Ut, m, A%Vt, n, beta, tempM, m)
		
		allocate(C)
		C%full = .false.
		if (Asmallest) then
			n2 = size(B%Vt,2)
			C%pointU = .true.
			C%pointV = .false.
			allocate(C%Vt(n,size(B%Vt,2)))
			C%Ut => A%Ut
			call dgemm('T','N', n, n2, m, alpha, tempM, m, B%Vt, m, beta, C%Vt, n)
		else
			n2 = size(A%Vt,2)
			allocate(C%Ut(m,n2))
			C%pointU = .false.
			C%pointV = .true.
			C%Vt => B%Vt
			call dgemm('N','N', m, n2, n, alpha, tempM, m, A%Ut, n, beta, C%Ut, m)
		endif
		
	end subroutine

end module
