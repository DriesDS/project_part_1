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

		call M_dealloc(A)
		call M_dealloc(B)
	case('lowrank')
		call lowrankcall(currarg)
	case('matprod')
		call matrixReader(A)
		call matrixReader(B)

		call matprod(A,B,C)

		call matrixWriter(C)

		call M_dealloc(A)
		call M_dealloc(B)
		call M_dealloc(C)
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

	subroutine lowrankcall(currarg)
		type(Matrix), pointer :: A,B
		character(len=8) :: cmd
		integer :: currarg
		integer :: i, rank
		double precision :: eps
		
		call getarg(currarg,cmd)
		currarg = currarg+1
		select case(cmd)
		case ('rank')
			call getarg(currarg, cmd)
			currarg = currarg+1
			read(cmd,*) rank
			call matrixReader(A)
			call lowrank(A, B, rank=rank)
		case ('epsabs')
			call getarg(currarg, cmd)
			currarg = currarg+1
			read(cmd,*) eps
			call matrixReader(A)
			call lowrank(A, B, rank=rank)
			call lowrank(epsabs=eps)
		case ('epsrel')
			call getarg(currarg, cmd)
			currarg = currarg+1
			read(cmd,*) eps
			call matrixReader(A)
			call lowrank(A, B, rank=rank)
			call lowrank(epsrel=eps)
		case default
			call help()
		end select

		call M_dealloc(A)
		call M_dealloc(B)
	
	end subroutine

	subroutine M_dealloc(A)
		type(Matrix), pointer :: A

		if (A%pointU) then
			nullify(A%Ut)
		else
			deallocate(A%Ut)
		endif

		if (.not.A%full) then
			if (A%pointV) then
				nullify(A%Vt)
			else
				deallocate(A%Vt)
			endif
		endif

		deallocate(A)

	end subroutine

end program
