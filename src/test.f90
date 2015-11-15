module testmod

use matrixconverter
use fullmod
use lowrankmod

implicit none

contains

	subroutine test()
		type(Matrix), pointer :: M1, M2, M3, M4, xfield
		integer :: d(22)
		open(7,file='../tests/M1.in')
		call matrixReader(M1, 7)
		close(7)
		open(7,file='../tests/M2.in')
		call matrixReader(M2, 7)
		close(7)
		open(7,file='../tests/M3.in')
		call matrixReader(M3, 7)
		close(7)
		open(7,file='../tests/M4.in')
		call matrixReader(M4, 7)
		close(7)
		open(7,file='../tests/xfield.in')
		call matrixReader(xfield, 7)
		close(7)

		
		
	end subroutine

end module

