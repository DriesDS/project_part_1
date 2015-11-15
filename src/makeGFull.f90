module makeGFullmod

use matrixconverter

implicit none

contains

	subroutine makeGFull(N)
		type(Matrix), pointer :: G
		integer :: i,j
		integer, intent(in) :: N
		double precision, parameter :: pi=4d0*atan(1d0)

		allocate(G)
		G%full = .true.
		allocate(G%Ut(N,N))

		call calculateG(G%Ut, N)

		call matrixWriter(G)

		deallocate(G)

	end subroutine

	subroutine calculateG(mat, N)
		double precision, intent(inout) :: mat(N,N)
		integer :: i,j
		double precision, parameter :: pi=4d0*atan(1d0)
		integer, intent(in) :: N

		do j = 1,N
			do i = 1,N
				mat(i,j) = -0.5*log( 2*(1-cos((1.0+2*j-2*i)/N*pi)) )
			enddo
		enddo
	end subroutine

end module
