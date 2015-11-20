module makeGFullmod

use matrixconverter

implicit none

contains

	subroutine makeGFull(G, N)
		type(Matrix), pointer, intent(inout) :: G
		integer, intent(in) :: N

		allocate(G)
		G%full = .true.
		G%pointU = .false.
		allocate(G%Ut(N,N))

		call calculateGt(G%Ut, N)

	end subroutine

	subroutine calculateGt(mat, N)
		integer, intent(in) :: N
		double precision, intent(inout) :: mat(N,N)
		integer :: i,j
		double precision, parameter :: pi=4d0*atan(1d0)

		do i = 1,N
			do j = 0,N-1
				! inverted indices since we are calculating Gt
				mat(j+1,i) = -0.5d0*log( 2*(1-cos((+1d0+2d0*j-2d0*i)*pi/N)) )
			enddo
		enddo
	end subroutine

end module
