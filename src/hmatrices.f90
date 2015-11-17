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
		call lowrankcall(A, B, currarg)
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
		call makeGFull(A,N)
		call matrixWriter(A)
		call M_dealloc(A)
	case('solveIntFull')
		call getarg(2,cmd)
		read(cmd,*) N
		call matrixReader(B)
		call solveIntFull(A, B, N)
		call matrixWriter(A)
		call M_dealloc(A)
		call M_dealloc(B)
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

	subroutine lowrankcall(A, B, currarg)
		type(Matrix), pointer :: A,B
		character(len=24) :: cmd, value
		integer :: currarg
		integer :: i, rank, effrank
		double precision :: eps
		
		call getarg(currarg,cmd)
		currarg = currarg+1
		call getarg(currarg,value)
		currarg = currarg+1

		call matrixReader(A)
		select case(cmd)
		case ('rank')
			read(value,*) rank
			call lowrank(A, B, effrank, rank=rank)
		case ('epsabs')
			read(value,*) eps
			call lowrank(A, B, effrank, epsabs=eps)
		case ('epsrel')
			read(value,*) eps
			call lowrank(A, B, effrank, epsrel=eps)
		case default
			call help()
		end select

		call matrixWriter(B, effrank)

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
