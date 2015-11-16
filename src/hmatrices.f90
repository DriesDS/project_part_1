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
	
	type(Matrix), pointer :: A,B,C
	integer :: nbarg, currarg, N
	character(len=16) :: cmd
	logical :: timing

	currarg = 1
	timing = .false.

	call getarg(currarg,cmd)
	currarg = currarg+1
	nbarg = iargc()

	if (cmd == '-t') then
		!start timing
		call getarg(currarg,cmd)
		currarg = currarg+1
		timing = .true.
	endif
	
	select case(cmd)
	case('help')
		call help()
	case('test')
		call test()
	case('full')
		call matrixReader(A)
		call full(A,B)
		call matrixWriter(B)
		deallocate(A%Ut)
		if (allocated(B%Ut)) deallocate(B%Ut)
		if (allocated(A%Vt)) deallocate(A%Vt)
		if (allocated(B%Vt)) deallocate(B%Vt)
		deallocate(A,B)
	case('lowrank')
		call lowrankcall(nbarg)
	case('matprod')
		call matrixReader(A)
		call matrixReader(B)
		call matprod(A,B,C)
		call matrixWriter(C)
		deallocate(A%Ut, B%Ut, C%Ut)
		if (allocated(A%Vt)) deallocate(A%Vt)
		if (allocated(B%Vt)) deallocate(B%Vt)
		if (allocated(C%Vt)) deallocate(C%Vt)
		deallocate(A,B,C)
	case('makeGFull')
		call getarg(currarg,cmd)
		read(cmd,*) N
		call makeGFull(B,N)
		call matrixWriter(B)
	case('solveIntFull')
		call getarg(2,cmd)
		read(cmd,*) N
		call solveIntFull(N)
	case('plotField')
		call plotField()
	case('readtest')
		call readtest()
	case default
		call wrongArg()
	end select

	if (timing) then
		! stop timing
	endif

contains

	subroutine lowrankcall(nbarg)
		character(len=8) :: currarg, arg1, arg2, arg3
		integer, intent(in) :: nbarg
		integer :: i, rank
		double precision :: epsabs, epsrel
		
		if (mod(nbarg,2) == 0) then
			call wrongArg()
			return
		endif
		
		call getarg(2,currarg)
		select case(currarg)
		case ('rank')
			call getarg(3, currarg)
			read(currarg,*) rank
			call lowrank(rank=rank)
		case ('epsabs')
			call getarg(3, currarg)
			read(currarg,*) epsabs
			call lowrank(epsabs=epsabs)
		case ('epsrel')
			call getarg(3, currarg)
			read(currarg,*) epsrel
			call lowrank(epsrel=epsrel)
		case default
			call wrongArg()
		end select
	
	end subroutine

end program
