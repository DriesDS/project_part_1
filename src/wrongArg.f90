module wrongArgmod

implicit none

contains

	subroutine wrongArg()
		write(*,*) 'wrong argumentlist'
	end subroutine

end module

