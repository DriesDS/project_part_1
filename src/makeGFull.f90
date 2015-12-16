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

		call calculateGt(G%Ut, N, 1, 1, N)

	end subroutine

	subroutine calculateGt(mat, N, beginx, beginy, s)
		integer, intent(in) :: N
		double precision, intent(inout) :: mat(s, s)
		integer :: i,j
		double precision, parameter :: pi=4d0*atan(1d0)

		! aan de onderstaande formule voor de eerste kolom van de matrix zien we
		! dat kolom j dezelfde elementen bevat als kolom 1, enkel j-1 plaatsen
		! naar onder geshift. Of kolom j bevat dezelfde elementen als kolom j-1
		! maar dan 1 plaats naar onder geshift. We berekenen dus enkel de eerste
		! kolom en kopieren die dan telkens in de volgende kolommen.
		! inverted indices since we are calculating Gt
		i = beginx
		do j=beginy-1,beginy+s-2
			mat(j+1,i) = -0.5d0*log( 2*(1-cos((+1d0+2d0*j-2d0*i)*pi/N)) )
		enddo

		do i = beginx+1,beginx+s-1
			mat(1,i) = mat(N,i-1)
			mat(2:N,i) = mat(1:N-1,i-1)
		enddo
	end subroutine

end module
