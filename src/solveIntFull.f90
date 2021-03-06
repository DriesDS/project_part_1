module solveIntFullmod

use makeGFullmod
use matrixconverter

implicit none

contains

	subroutine solveIntFull(U, X, N)
		integer, intent(in) :: N
		type(Matrix), pointer :: X, U
		double precision :: A(N,N), B(N,1)
		integer :: IPIV(N)
		double precision, parameter :: pi=4d0*atan(1d0)
		integer :: i, j, info
		
		call calculateGt(A,N,1,1,N)
		A = Transpose(A)
		forall (i=1:N) B(i,1) = cos(2d0*i*pi/N)
		IPIV = 0
		
		call dgesv(N,1,A,N,IPIV,B,N,info)

		allocate(U)
		U%full = .true.
		U%pointU = .false.
		allocate(U%Ut(size(X%Ut,1),1))

		U%Ut = 0
		do i=1,size(U%Ut,1)
			do j=1,N
				! na het herwerken van de gegeven formule:
				U%Ut(i,1) = U%Ut(i,1)-B(j,1)*1d0/2*log(1d0+X%Ut(i,1)**2+X%Ut(i,2)**2 &
					-2d0*(X%Ut(i,1)*cos((2d0*j-1)*pi/N)+X%Ut(i,2)*sin((2d0*j-1)*pi/N)))
			enddo
		enddo

	end subroutine

end module
