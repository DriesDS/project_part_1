module helpmod

implicit none

integer, parameter :: maxlines=1000

contains

	subroutine help()
		character(len=128) :: textline
		integer :: i

		open(7, file = 'help.txt')
		do i = 1, maxlines
			read(7,'(a)') textline
			write(*,'(a)') textline
		enddo

	end subroutine

end module

