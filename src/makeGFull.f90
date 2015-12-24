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
		integer, intent(in) :: N, beginx, beginy, s
		double precision, dimension(s,s), intent(inout) :: mat
		integer :: i,j
		double precision, parameter :: pi=4d0*atan(1d0)

		if (.not. size(mat,1) == s) then
			if (.not. size(mat,2) == s) then
				write(0,*) 'this can''t be possible, the matrix should be sxs'
			endif
		endif

		! aan de onderstaande formule voor de eerste kolom van de matrix zien we
		! dat kolom j dezelfde elementen bevat als kolom 1, enkel j-1 plaatsen
		! naar onder geshift. Of kolom j bevat dezelfde elementen als kolom j-1
		! maar dan 1 plaats naar onder geshift. We berekenen dus enkel de eerste
		! kolom en kopieren die dan telkens in de volgende kolommen.
		! inverted indices since we are calculating Gt
		i = beginy
		do j=0,s-1
			mat(j+1,1) = -0.5d0*log( 2*(1-cos((1d0+2d0*(j+beginx)-2d0*i)*pi/N)) )
		enddo

		do i = 2,s
			mat(1,i) = mat(s,i-1)
			mat(2:s,i) = mat(1:s-1,i-1)
		enddo
	end subroutine

end module
