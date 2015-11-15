module solveIntFullmod

use makeGFullmod
use matrixconverter

implicit none

contains

	subroutine solveIntFull(N)
		type(Matrix), pointer :: X, U
		double precision :: A(N:N), B(N), IPIV(N)
		double precision, parameter :: pi=4.0*atan(1.0)
		integer, intent(in) :: N
		integer :: i, j, info
		write(*,*) 'calculating G'
		A = 0
		call calculateG(A,N)
		write(*,*) 'G calculated'
		forall (i=1:N) b(i) = cos(2*i*pi/N)
		forall (i=1:N) IPIV(i) = (i)
		write(*,*) 'b and ipiv set'
		call dgetrs('N', N, 1, A, N, IPIV, B, N, info)
		write(*,*) B

		call matrixReader(X)

		allocate(U)
		U%full = .true.
		allocate(U%Ut(size(X%Ut,1),1))

		U%Ut = 0
		do i=1,size(U%Ut,1)
			do j=1,N
				U%Ut(i,1) = U%Ut(i,1)-B(j)*1.0/2*log(1+X%Ut(j,1)**2+X%Ut(j,2)**2-2*(X%Ut(j,1)*cos((2*j+1)*pi/N)+X%Ut(j,2)*sin((2*j+1)*pi/N)))
			enddo
		enddo

		call matrixWriter(U)

		deallocate(X,U)

	end subroutine

end module

