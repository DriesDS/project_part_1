module solveIntFullmod

use makeGFullmod
use matrixconverter

implicit none

contains

	subroutine solveIntFull(N)
		type(Matrix), pointer :: X, U
		double precision :: A(N:N), B(N), IPIV(N)
		integer, intent(in) :: N
		integer :: i, info

		call calculateG(A,N)
		forall (i=1:N) b(i) = cos(2*i*pi/N)
		forall (i=1:N) IPIV(i) = (i)

		call dgetrs('N', N, 1, A, N, IPIV, B, N, info)

		call matrixReader(X)

		allocate(U)
		U%full = .true.
		allocate(U%Ut(size(X%Ut,1),1))

		U%Ut = 0
		do i=1,size(U%Ut,1)
			do j=1,N
				U%Ut(i) = U%Ut(i)-B(j)*1.0/2*log(1+pow(X%Ut(j,1),2)+pow(X%Ut(j,2),2)-2*(X%Ut(j,1)*cos((2*j+1)*pi/N)+X%Ut(j,2)*sin((2*j+1)*pi/N)))
			enddo
		enddo

		call matrixWriter(U)

		deallocate(X,U)

	end subroutine

end module

