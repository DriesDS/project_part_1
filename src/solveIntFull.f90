module solveIntFullmod

use makeGFullmod
use matrixconverter

implicit none

contains

	subroutine solveIntFull(N)
		integer, intent(in) :: N
		type(Matrix), pointer :: X, U
		double precision :: A(N,N), B(N,1)
		integer :: IPIV(N)
		double precision, parameter :: pi=4d0*atan(1d0)
		integer :: i, j, info
		call calculateGt(A,N)
		A = Transpose(A)
		forall (i=1:N) B(i,1) = cos(2d0*i*pi/N)
		IPIV = 0
		call dgetrf(N,N,A,N,IPIV,info)
		call dgetrs('N', N, 1, A, N, IPIV, B, N, info)
		call matrixReader(X)

		allocate(U)
		U%full = .true.
		allocate(U%Ut(size(X%Ut,2),1))

		U%Ut = 0
		do i=1,size(U%Ut,2)
			do j=1,N
				U%Ut(i,1) = U%Ut(i,1)-B(j,1)*1d0/2*log(1d0+X%Ut(j,1)**2+X%Ut(j,2)**2-2d0*(X%Ut(j,1)*cos((2d0*j+1)*pi/N)+X%Ut(j,2)*sin((2d0*j+1)*pi/N)))
			enddo
		enddo

		call matrixWriter(U)

		deallocate(X,U)

	end subroutine

end module

