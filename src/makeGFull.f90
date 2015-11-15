module makeGFullmod

use matrixconverter

implicit none

contains

	subroutine makeGFull(N)
		type(Matrix), pointer :: G
		integer :: i,j
		integer, intent(in) :: N
		double precision, parameter :: pi=4.0*atan(1.0)

		allocate(G)
		G%full = .true.
		allocate(G%Ut(N,N))
		do j = 1,N
			do i = 1,N
				! na uitschrijven van de kwadraten en gebruik te maken van 
				! goniometrische formules.
				G%Ut(i,j) = -0.5*log( 2*(1-cos((1.0+2*j-2*i)/N*pi)) )
			enddo
		enddo

		call matrixWriter(G)

		deallocate(G)

	end subroutine

end module
