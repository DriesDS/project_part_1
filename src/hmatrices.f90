program hmatrices
	
	use testmod
	use helpmod
	use lowrankmod
	use fullmod
	use matprodmod
	use makeGFullmod
	use solveIntFullmod
	use plotFieldmod
	use wrongArgmod
	use readtestmod
	
	implicit none
	
	integer :: nbarg
	character(len=16) :: cmd
	
	call getarg(1,cmd)
	nbarg = iargc()
	
	select case(cmd)
	case('help')
		call help()
	case('test')
		call test()
	case('full')
		call full()
	case('lowrank')
		call lowrankcall(nbarg)
	case('matprod')
		call matprod()
	case('makeGFull')
		call makeGFull()
	case('solveIntFull')
		call solveIntFull()
	case('plotField')
		call plotField()
	case('readtest')
		call readtest()
	case default
		call wrongArg()
	end select

contains

	subroutine lowrankcall(nbarg)
		character(len=8) :: currarg, arg1, arg2, arg3
		integer, intent(in) :: nbarg
		integer :: i, rank
		double precision :: epsabs, epsrel
		
		!if (mod(nbarg,2) == 0) call wrongArg()
		!
		!do i = 1, (nbarg-1)/2
		!	call getarg(i*2,currarg)
		!	select case(currarg)
		!	case ('rank')
		!		call getarg(i*2+1, currarg)
		!		read(currarg,*) rank
		!	case ('epsabs')
		!		call getarg(i*2+1, currarg)
		!		read(currarg,*) epsabs
		!	case ('epsrel')
		!		call getarg(i*2+1, currarg)
		!		read(currarg,*) epsrel
		!	case default
		!		call wrongArg()
		!	end select
		!enddo

		rank = 0
		epsabs = 0.0
		epsrel = 0.0

		call lowrank(rank, epsabs, epsrel)
	
	end subroutine

end program
