module fullmod

use matrixconverter

implicit none

contains

	subroutine full(A,B)
		integer :: m, n, k
		type(Matrix), pointer, intent(in) :: A
		type(Matrix), pointer, intent(inout) :: B
		double precision :: alpha, beta

		if (A%full) then
			allocate(B)
			B%full = .true.
			B%pointU = .true.
			B%pointV = .false.
			B%Ut => A%Ut
		else
			! B^t = A%Vt^t * A%Ut
			m = size(A%Vt,2)
			n = size(A%Ut,2)
			k = size(A%Ut,1)
			alpha = 1.0
			beta = 0.0

			allocate(B)
			B%full = .true.
			B%pointU = .false.
			allocate(B%Ut(m,n))

			call dgemm('T','N',m,n,k,alpha,A%Vt,k,A%Ut,k,beta,B%Ut,m)

		endif

	end subroutine

end module

