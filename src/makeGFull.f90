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

		call calculateG(G%Ut, N)

		call matrixWriter(G)

		deallocate(G)

	end subroutine

	subroutine calculateG(mat, N)
		double precision, intent(inout) :: mat(N,N)
		integer :: i,j
		double precision, parameter :: pi=4.0*atan(1.0)
		integer, intent(in) :: N
		write(*,*) 'calculating G with N= ', N
		write(*,*) 'size matrix mat: ', size(mat,1), size(mat,2)
		do j = 1,N
			do i = 1,N
				write(*,*) 'j=',j,'  i=',i
				! na uitschrijven van de kwadraten en gebruik te maken van 
				! goniometrische formules.
				write(*,*) cos((1.0+2*j-2*i)/N*pi)
				mat(i,j) = -0.5*log( 2*(1-cos((1.0+2*j-2*i)/N*pi)) )
			enddo
		enddo
		write(*,*) 'G calculated', mat
	end subroutine

end module
