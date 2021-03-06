program hmatrices
	
	use testmod
	use helpmod
	use lowrankmod
	use fullmod
	use matprodmod
	use makeGFullmod
	use makeGHmatmod
	use solveIntFullmod
	use plotFieldmod
	use vecProdHmatmod

	implicit none

	type(Matrix), pointer :: A,B,C
	type(HMatrix), pointer :: AH
	integer :: currarg, N, elems
	character(len=16) :: cmd
	logical :: timing
	double precision :: start, total_time, y

	currarg = 1
	timing = .false.

	call getarg(currarg,cmd)
	currarg = currarg+1

	if (cmd == '-t') then
		call CPU_TIME(start)
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
		call getarg(currarg,cmd)
		currarg = currarg+1
		read(cmd,*) N
		call matrixReader(B)
		call solveIntFull(A, B, N)
		call matrixWriter(A)
		call M_dealloc(A)
		call M_dealloc(B)
	case('plotField')
		call getarg(currarg,cmd)
		currarg = currarg+1
		read(cmd,*) N
		call matrixReader(B)
		call plotField(B,N)
		call M_dealloc(B)
	case('makeGHmat')
		allocate(AH)
		call getarg(currarg,cmd)
		currarg = currarg+1
		read(cmd,*) N
		call getarg(currarg,cmd)
		currarg = currarg+1
		read(cmd,*) y
		call makeGHmat(AH, N, y, elems)
		write(*,*) elems
		call Hm_dealloc(AH)
	case('vecProdHmat')
		allocate(AH)
		call getarg(currarg,cmd)
		currarg = currarg+1
		read(cmd,*) N
		call getarg(currarg,cmd)
		currarg = currarg+1
		read(cmd,*) y
		call matrixReader(A)
		call vecProdHmat(AH, N, y, A, B, elems)
		call Hm_dealloc(AH)
		call matrixWriter(B)
		call M_dealloc(A)
		call M_dealloc(B)
		write(*,*) elems
	case default
		call help()
	end select

	if (timing) then
		call CPU_TIME(total_time)
		total_time = total_time-start
		write(0,*) 'CPU_TIME: ', total_time
	endif

contains

	subroutine lowrankcall(A, B, currarg)
		type(Matrix), pointer :: A,B
		character(len=24) :: cmd, value
		integer :: currarg
		integer :: rank
		double precision :: eps
		
		call getarg(currarg,cmd)
		currarg = currarg+1
		call getarg(currarg,value)
		currarg = currarg+1

		call matrixReader(A)
		select case(cmd)
		case ('rank')
			read(value,*) rank
			call lowrank(A, B, rank=rank)
		case ('epsabs')
			read(value,*) eps
			call lowrank(A, B, epsabs=eps)
		case ('epsrel')
			read(value,*) eps
			call lowrank(A, B, epsrel=eps)
		case default
			call lowrank(A, B)
		end select

		call matrixWriter(B)

		call M_dealloc(A)
		call M_dealloc(B)
	
	end subroutine

end program
