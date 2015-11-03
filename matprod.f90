module matprodmod
use matrixconverter

implicit none

contains

	subroutine matprod()
		integer :: n
		type(Matrix), pointer :: mat
		!write(*,*) 'in subroutine matprod'
		call matrixReader(mat)
		call matrixWriter(mat)
		deallocate(mat%Ut)
		if (.not. mat%full) deallocate(mat%Vt)
		deallocate(mat)
	end subroutine

end module

