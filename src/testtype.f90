program testtype

implicit none

call test_this()

type TT
	logical :: endTT
	type(TT), pointer :: subT
end type

contains

	recursive subroutine test_this()
		integer, parameter :: n=10

		do i = 1,n
			allocate(subT)
			call test_this(testT%subT)
		enddo

	end subroutine

end program