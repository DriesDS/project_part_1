module fullmod

use matrixconverter

implicit none

contains

	subroutine full()
		integer :: m, n, k
		type(Matrix), pointer :: A, B
		double precision :: alpha, beta

		call matrixReader(A)

		if (A%full) then
			call matrixWriter(A)
		else
			! B^t = A%Vt^t * A%Ut
			m = size(A%Vt,2)
			n = size(A%Ut,2)
			k = size(A%Ut,1)
			alpha = 1.0
			beta = 0.0

			allocate(B)
			B%full = .true.
			allocate(B%Ut(m,n))

			call dgemm('T','N',m,n,k,alpha,A%Vt,k,A%Ut,k,beta,B%Ut)

			call matrixWriter(B)
		endif

		deallocate(A,B)

	end subroutine

end module

