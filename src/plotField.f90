module plotFieldmod

use matrixconverter
use solveintfullmod

implicit none

contains

	subroutine plotField(X, N)
		Type(Matrix), pointer :: X, U
		integer, intent(in) :: N
		integer :: i

		call solveIntFull(U,X,N)

		open(10,file='matlab/plotfieldpoints.out')
		do i = 1, size(X%Ut,1)
			write(10,*) X%Ut(i,:), U%Ut(i,:)
		enddo
		close(10)

		call SYSTEM('cat ./matlab/plotField.m | matlab -nosplash -noFigureWindows >matlab/log.txt')

		call M_dealloc(U)

	end subroutine

end module

