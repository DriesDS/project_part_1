module readtestmod

use matrixconverter

implicit none

contains

	subroutine readtest()
		type(Matrix), pointer :: A,B

		call matrixReader(A)
		call matrixReader(B)

		call matrixWriter(A)
		call matrixWriter(B)

		deallocate(A)
		deallocate(B)
	end subroutine

end module 

