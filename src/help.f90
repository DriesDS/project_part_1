module helpmod

implicit none

integer, parameter :: maxlines=1000

contains

	subroutine help()
		character(len=128) :: textline
		integer :: i, io

		open(10, file = './src/help.txt')
		do i = 1, maxlines
			read(10,'(a)',iostat=io) textline
			if (io<0) exit
			write(0,'(a)', advance='NO') textline
		enddo

	end subroutine

end module
