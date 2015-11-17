module testmod

use matrixconverter
use fullmod
use lowrankmod
use matprodmod
use makeGFullmod
use solveIntFullmod

implicit none

contains

	subroutine test()
		call test1()
		call test2()
		call test3()
	end subroutine

	subroutine test1()
		type(Matrix), pointer :: M1, M2, M3, M4, xfield, M1u, M2u, M3u, M4u, onesten, onestwo, onetillten
		integer, parameter :: nbtests=22
		integer :: i, d(nbtests), effrank

		call loadmatrices(M1, M2, M3, M4, xfield)
		allocate(onesten)
		allocate(onestwo)
		onesten%full = .true.
		onesten%pointU = .false.
		allocate(onesten%Ut(1,10))
		onesten%Ut = 1
		onestwo%full = .true.
		onestwo%pointU = .false.
		allocate(onestwo%Ut(2,1))
		onestwo%Ut = 1
		allocate(onetillten)
		onetillten%full = .true.
		allocate(onetillten%Ut(1,10))
		forall (i=1:10) onetillten%Ut(1,i) = i

		! test 1
		call full(M1, M1u)
		call unit_digit(d(1), M1u, 2, 1)
		call Mu_dealloc(M1u)

		! test 2
		call full(M1, M1u)
		call lowrank(M1u, M2u, effrank, rank=2)
		call full(M2u, M3u)
!		call matrixWriter(M3u)
		call unit_digit(d(2), M3u, 1,6)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)
		call Mu_dealloc(M3u)

		! test 3
		call full(M1, M1u)
		call lowrank(M1u, M2u, effrank, epsrel=0.7d0)
		call full(M2u, M3u)
		call unit_digit(d(3), M3u, 7,5)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)
		call Mu_dealloc(M3u)

		! test 4
		call full(M1, M1u)
		call lowrank(M1u, M2u, effrank, epsabs=60.0d0)
		call full(M2u, M3u)
		call unit_digit(d(4), M3u, 1,4)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)
		call Mu_dealloc(M3u)

		! test 5
		call full(M1, M1u)
		call matprod(M1u, onesten, M2u)
		call unit_digit(d(5), M2u, 3,1)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)

		! test 6
		call lowrank(M1, M1u, effrank, rank=3)
		call full(M1u, M2u)
		call unit_digit(d(6), M2u, 2, 2)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)

		! test 7
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M1, M1u, effrank, epsrel=0.5d0)
		call full(M1u, M2u)
		call unit_digit(d(7), M2u, 8, 10)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)

		! test 8
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M1, M1u, effrank, epsabs=80.0d0)
		call full(M1u, M2u)
		call unit_digit(d(8), M2u, 5, 7)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)

		! test 9
		call loadmatrices(M1, M2, M3, M4, xfield)
		call matprod(M1, M2, M1u)
		call lowrank(M1u, M2u, effrank, rank=4)
		call full(M2u, M3u)
		call unit_digit(d(9), M3u, 5, 2)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)
		call Mu_dealloc(M3u)

		! test 10
		call loadmatrices(M1, M2, M3, M4, xfield)
		call matprod(M1, M2, M1u)
		call full(M1u, M2u)
		call unit_digit(d(10), M2u, 2, 10)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)

		! test 11
		call matprod(onestwo, xfield, M1u)
		call unit_digit(d(11), M1u, 1,1)
		call Mu_dealloc(M1u)

! 		! ----------------- next sequence
		
		! test 1
		call full(M3, M3u)
		call unit_digit(d(12), M3u, 10, 10)
		call Mu_dealloc(M3u)

		! test 2
		call full(M3, M3u)
		call lowrank(M3u, M4u, effrank, rank=5)
		call full(M4u, M1u)
		call unit_digit(d(13), M1u, 7, 1)
		
		call Mu_dealloc(M3u)
		call Mu_dealloc(M4u)
		call Mu_dealloc(M1u)

		! test 3
		call full(M3, M3u)
		call lowrank(M3u, M4u, effrank, epsrel=0.2d0)
		call full(M4u, M1u)
		call unit_digit(d(14), M1u, 1, 7)
		call Mu_dealloc(M3u)
		call Mu_dealloc(M4u)
		call Mu_dealloc(M1u)

		! test 4
		call full(M3, M3u)
		call lowrank(M3u, M4u, effrank, epsabs=51.0d0)
		call full(M4u, M1u)
		call unit_digit(d(15), M1u, 2, 8)
		call Mu_dealloc(M3u)
		call Mu_dealloc(M4u)
		call Mu_dealloc(M1u)
		
		! test 5
		call matprod(M3, onetillten, M4u)
		call full(M4u, M1u)
		call unit_digit(d(16), M1u, 9, 1)
		call Mu_dealloc(M4u)
		call Mu_dealloc(M1u)

		! test 6
		call lowrank(M3, M3u, effrank, rank=2)
		call full(M3u, M1u)
		call unit_digit(d(17), M1u, 3, 2)
		call Mu_dealloc(M3u)
		call Mu_dealloc(M1u)

		! test 7
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M3, M3u, effrank, epsrel=0.8d0)
		call full(M3u, M1u)
		call unit_digit(d(18), M1u, 5, 5)
		call Mu_dealloc(M3u)
		call Mu_dealloc(M1u)

		! test 8
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M3, M4u, effrank, epsabs=61.125d0)
		call full(M4u, M1u)
		call unit_digit(d(19), M1u, 7, 6)
		call Mu_dealloc(M4u)
		call Mu_dealloc(M1u)

		! test 9
		call loadmatrices(M1, M2, M3, M4, xfield)
		call matprod(M3, M4, M4u)
		call lowrank(M4u, M1u, effrank, rank=10)
		call full(M1u, M2u)
		call unit_digit(d(20), M2u, 2, 1)
		call Mu_dealloc(M4u)
		call Mu_dealloc(M1u)
		call Mu_dealloc(M2u)

		! test 10
		call loadmatrices(M1, M2, M3, M4, xfield)
		call matprod(M3, M4, M4u)
		call full(M4u, M1u)
		call unit_digit(d(21), M1u, 2, 4)
		call Mu_dealloc(M4u)
		call Mu_dealloc(M1u)

		! test 11
		call full(xfield, M1u)
		call unit_digit(d(22), M1u, 2, 10)
		call Mu_dealloc(M1u)

		call Mu_dealloc(M1)
		call Mu_dealloc(M2)
		call Mu_dealloc(M3)
		call Mu_dealloc(M4)
		call Mu_dealloc(xfield)
		call Mu_dealloc(onetillten)
		call Mu_dealloc(onesten)
		call Mu_dealloc(onestwo)


		write(*,'(t10, i0,",",10(i0),t50,"  <-- This should be Pi")') d(1:11) 
		write(*,'(t10, 7(i0,x),2(2(i0),x),t50,"  <-- This should be the row of Fibonacci")') d(12:22) 
		write(*,*)
		write(*,*)

	end subroutine

	subroutine test2()
		type(Matrix), pointer :: G, X

		call makeGFull(G, 5)

		write(*,'(e15.3,"   This should be equal to 0.481")') G%Ut(1,1)
		write(*,'(e15.4,"   This should be equal to -0.4812")') G%Ut(1,2)
		write(*,'(e15.3,"   This should be equal to -0.693")') G%Ut(1,3)

		call Mu_dealloc(G)

		allocate(X)
		X%full = .true.
		X%pointU = .false.
		allocate(X%Ut(2,1))
		X%Ut(1,1) = 0.7
		X%Ut(2,1) = -0.9

		call solveIntFull(G,X,40)
		write(*,*) G%Ut(1,1)

		call Mu_dealloc(G)
		call Mu_dealloc(X)
	end subroutine

	subroutine test3()

	end subroutine

	subroutine loadmatrices(M1, M2, M3, M4, xfield)
		type(Matrix), pointer :: M1, M2, M3, M4, xfield		

		open(10,file='../tests/M1.in')
		call matrixReader(M1, 10)
		close(10)
		open(10,file='../tests/M2.in')
		call matrixReader(M2, 10)
		close(10)
		open(10,file='../tests/M3.in')
		call matrixReader(M3, 10)
		close(10)
		open(10,file='../tests/M4.in')
		call matrixReader(M4, 10)
		close(10)
		open(10,file='../tests/xfield.in')
		call matrixReader(xfield, 10)
		close(10)

	end subroutine


	subroutine unit_digit(digit, A, row, col)
		type(Matrix), pointer :: A
		integer :: row, col, digit

		digit = mod(abs(int(A%Ut(col, row))),10)
	end subroutine

	subroutine Mu_dealloc(A)
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

end module

