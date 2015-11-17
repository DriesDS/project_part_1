module makeGFullmod

use matrixconverter

implicit none

contains

	subroutine makeGFull(G, N)
		type(Matrix), pointer, intent(inout) :: G
		integer :: i,j
		integer, intent(in) :: N

		allocate(G)
		G%full = .true.
		G%pointU = .false.
		allocate(G%Ut(N,N))

		call calculateGt(G%Ut, N)

	end subroutine

	subroutine calculateGt(mat, N)
		double precision, intent(inout) :: mat(N,N)
		integer :: i,j
		double precision, parameter :: pi=4d0*atan(1d0)
		integer, intent(in) :: N

		do i = 1,N
			do j = 1,N
				! inverted indices since we are calculating Gt
				mat(j,i) = -0.5d0*log( 2*(1-cos((-1d0+2d0*j-2d0*i)/N*pi)) )
			enddo
		enddo
	end subroutine

end module
