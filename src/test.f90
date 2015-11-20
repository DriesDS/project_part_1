module testmod

use matrixconverter
use fullmod
use lowrankmod
use matprodmod
use makeGFullmod
use solveIntFullmod
use plotfieldmod

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
		integer :: i, d(nbtests)

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
		call M_dealloc(M1u)

		! test 2
		call full(M1, M1u)
		call lowrank(M1u, M2u, rank=2)
		call full(M2u, M3u)
		call unit_digit(d(2), M3u, 1,6)
		call M_dealloc(M1u)
		call M_dealloc(M2u)
		call M_dealloc(M3u)

		! test 3
		call full(M1, M1u)
		call lowrank(M1u, M2u, epsrel=0.7d0)
		call full(M2u, M3u)
		call unit_digit(d(3), M3u, 7,5)
		call M_dealloc(M1u)
		call M_dealloc(M2u)
		call M_dealloc(M3u)

		! test 4
		call full(M1, M1u)
		call lowrank(M1u, M2u, epsabs=60.0d0)
		call full(M2u, M3u)
		call unit_digit(d(4), M3u, 1,4)
		call M_dealloc(M1u)
		call M_dealloc(M2u)
		call M_dealloc(M3u)

		! test 5
		call full(M1, M1u)
		call matprod(M1u, onesten, M2u)
		call unit_digit(d(5), M2u, 3,1)
		call M_dealloc(M1u)
		call M_dealloc(M2u)

		! test 6
		call lowrank(M1, M1u, rank=3)
		call full(M1u, M2u)
		call unit_digit(d(6), M2u, 2, 2)
		call M_dealloc(M1u)
		call M_dealloc(M2u)

		! test 7
		!call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M1, M1u, epsrel=0.5d0)
		call full(M1u, M2u)
		call unit_digit(d(7), M2u, 8, 10)
		call M_dealloc(M1u)
		call M_dealloc(M2u)

		! test 8
		!call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M1, M1u, epsabs=80.0d0)
		call full(M1u, M2u)
		call unit_digit(d(8), M2u, 5, 7)
		call M_dealloc(M1u)
		call M_dealloc(M2u)

		! test 9
		!call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call loadmatrices(M1, M2, M3, M4, xfield)
		call matprod(M1, M2, M1u)
		call lowrank(M1u, M2u, rank=4)
		call full(M2u, M3u)
		call unit_digit(d(9), M3u, 5, 2)
		call M_dealloc(M1u)
		call M_dealloc(M2u)
		call M_dealloc(M3u)

		! test 10
		!call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call matprod(M1, M2, M1u)
		call full(M1u, M2u)
		call unit_digit(d(10), M2u, 2, 10)
		call M_dealloc(M1u)
		call M_dealloc(M2u)

		! test 11
		call matprod(onestwo, xfield, M1u)
		call unit_digit(d(11), M1u, 1,1)
		call M_dealloc(M1u)

		! ----------------- next sequence
		
		! test 1
		call full(M3, M3u)
		call unit_digit(d(12), M3u, 10, 10)
		call M_dealloc(M3u)

		! test 2
		call full(M3, M3u)
		call lowrank(M3u, M4u, rank=5)
		call full(M4u, M1u)
		call unit_digit(d(13), M1u, 7, 1)
		
		call M_dealloc(M3u)
		call M_dealloc(M4u)
		call M_dealloc(M1u)

		! test 3
		call full(M3, M3u)
		call lowrank(M3u, M4u, epsrel=0.2d0)
		call full(M4u, M1u)
		call unit_digit(d(14), M1u, 1, 7)
		call M_dealloc(M3u)
		call M_dealloc(M4u)
		call M_dealloc(M1u)

		! test 4
		call full(M3, M3u)
		call lowrank(M3u, M4u, epsabs=51.0d0)
		call full(M4u, M1u)
		call unit_digit(d(15), M1u, 2, 8)
		call M_dealloc(M3u)
		call M_dealloc(M4u)
		call M_dealloc(M1u)
		
		! test 5
		call matprod(M3, onetillten, M4u)
		call full(M4u, M1u)
		call unit_digit(d(16), M1u, 9, 1)
		call M_dealloc(M4u)
		call M_dealloc(M1u)

		! test 6
		call lowrank(M3, M3u, rank=2)
		call full(M3u, M1u)
		call unit_digit(d(17), M1u, 3, 2)
		call M_dealloc(M3u)
		call M_dealloc(M1u)

		! test 7
! 		call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M3, M3u, epsrel=0.8d0)
		call full(M3u, M1u)
		call unit_digit(d(18), M1u, 5, 5)
		call M_dealloc(M3u)
		call M_dealloc(M1u)

		! test 8
! 		call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call loadmatrices(M1, M2, M3, M4, xfield)
		call lowrank(M3, M4u, epsabs=61.125d0)
		call full(M4u, M1u)
		call unit_digit(d(19), M1u, 7, 6)
		call M_dealloc(M4u)
		call M_dealloc(M1u)

		! test 9
! 		call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call loadmatrices(M1, M2, M3, M4, xfield)
		call matprod(M3, M4, M4u)
		call lowrank(M4u, M1u, rank=10)
		call full(M1u, M2u)
		call unit_digit(d(20), M2u, 2, 1)
		call M_dealloc(M4u)
		call M_dealloc(M1u)
		call M_dealloc(M2u)

		! test 10
! 		call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call loadmatrices(M1, M2, M3, M4, xfield)
		call matprod(M3, M4, M4u)
		call full(M4u, M1u)
		call unit_digit(d(21), M1u, 2, 4)
		call M_dealloc(M4u)
		call M_dealloc(M1u)

		! test 11
		call full(xfield, M1u)
		call unit_digit(d(22), M1u, 2, 10)
		call M_dealloc(M1u)

		call dealloc_test_matrices(M1, M2, M3, M4, xfield)
		call M_dealloc(onetillten)
		call M_dealloc(onesten)
		call M_dealloc(onestwo)


		write(*,'(t10, i0,",",10(i0),t40,"  <-- This should be Pi")') d(1:11) 
		write(*,'(t10, 7(i0,1x),2(2(i0),1x),t40,"  <-- This should be the row of Fibonacci")') d(12:22)
		write(*,*)
		write(*,*)

	end subroutine

	subroutine test2()
		type(Matrix), pointer :: G, X, xfield
		integer, parameter :: N=40

		call makeGFull(G, 5)

		write(*,'(t10, e15.3,t40,"  <-- This should be equal to 0.481")') G%Ut(1,1)
		write(*,'(t10, e15.4,t40,"  <-- This should be equal to -0.4812")') G%Ut(1,2)
		write(*,'(t10, e15.3,t40,"  <-- This should be equal to -0.693")') G%Ut(1,3)

		call M_dealloc(G)

		allocate(X)
		X%full = .true.
		X%pointU = .false.
		allocate(X%Ut(2,1))
		X%Ut(1,1) = 0.7
		X%Ut(2,1) = -0.9

		call solveIntFull(G,X,N)
		write(*,'(t10, e15.3,t40,"  <-- This should be around 0.712")') G%Ut(1,1)

		call M_dealloc(G)
		call M_dealloc(X)

		open(10,file='./tests/xfield.in')
		call matrixReader(xfield, 10)
		close(10)
		call plotField(xfield, N)
		call M_dealloc(xfield)
	end subroutine

	subroutine test3()
		integer, parameter :: N=5
		integer :: i
		type(Matrix), pointer :: A, B, C, D, E, F
		double precision :: start, gestopt, time(N,2)

		call random_seed()
! 		do i = 5,10
! 			allocate(A)
! 			A%full = .true.
! 			A%pointU = .false.
! 			allocate(A%Ut(2**i,2**i))
! 			call random_number(A%Ut)
! 			allocate(B)
! 			B%full = .true.
! 			B%pointU = .false.
! 			allocate(B%Ut(2**i,2**i))
! 			call random_number(B%Ut)
			
! 			call cpu_time(start)
! 			call matprod(A,B,C)
! 			call cpu_time(gestopt)
! 			time(i,1) = gestopt-start

! !			deallocate(A%Ut,B%Ut,C%Ut)
! !			deallocate(C, A, B)

! 			allocate(D)
! 			D%full = .false.
! 			D%pointU = .false.
! 			D%pointV = .false.
! 			allocate(D%Ut(10,2**i))
! 			allocate(D%Vt(10,2**i))
! 			call random_number(D%Ut)
! 			call random_number(D%Vt)
! 			allocate(E)
! 			E%full = .false.
! 			E%pointU = .false.
! 			E%pointV = .false.
! 			allocate(E%Ut(10,2**i))
! 			allocate(E%Vt(10,2**i))
! 			call random_number(E%Ut)
! 			call random_number(E%Vt)
! 			call cpu_time(start)
! 			call matprod(D,E,F)
! 			call cpu_time(gestopt)
! 			time(i,2) = gestopt-start
! 			write(*,*) 2**i, time(i,1), time(i,2)
! 			write(0,*) 2**i, time(i,1), time(i,2)

! 		enddo


! 		do i = 5,10

! 			allocate(D)
! 			D%full = .false.
! 			D%pointU = .false.
! 			D%pointV = .false.
! 			allocate(D%Ut(2**i,512))
! 			allocate(D%Vt(2**i,512))
! 			call random_number(D%Ut)
! 			call random_number(D%Vt)
! 			allocate(E)
! 			E%full = .false.
! 			E%pointU = .false.
! 			E%pointV = .false.
! 			allocate(E%Ut(2**i,512))
! 			allocate(E%Vt(2**i,512))
! 			call random_number(E%Ut)
! 			call random_number(E%Vt)
! 			call cpu_time(start)
! 			call matprod(D,E,F)
! 			call cpu_time(gestopt)
! 			time(i,2) = gestopt-start
! 			write(*,*) 2**i, time(i,2)
! 			write(0,*) 2**i, time(i,2)

! 		enddo


	end subroutine

	subroutine loadmatrices(M1, M2, M3, M4, xfield)
		type(Matrix), pointer :: M1, M2, M3, M4, xfield		

		open(10,file='./tests/M1.in')
		call matrixReader(M1, 10)
		close(10)
		open(10,file='./tests/M2.in')
		call matrixReader(M2, 10)
		close(10)
		open(10,file='./tests/M3.in')
		call matrixReader(M3, 10)
		close(10)
		open(10,file='./tests/M4.in')
		call matrixReader(M4, 10)
		close(10)
		open(10,file='./tests/xfield.in')
		call matrixReader(xfield, 10)
		close(10)

	end subroutine

	subroutine dealloc_test_matrices(M1, M2, M3, M4, xfield)
		type(Matrix), pointer :: M1, M2, M3, M4, xfield		

		write(*,*) associated(M1), associated(M2), associated(M3), associated(M4), associated(xfield)
		call M_dealloc(M1)
		call M_dealloc(M2)
		call M_dealloc(M3)
		call M_dealloc(M4)
		call M_dealloc(xfield)

	end subroutine

	subroutine unit_digit(digit, A, row, col)
		type(Matrix), pointer :: A
		integer :: row, col, digit

		digit = mod(abs(int(A%Ut(col, row))),10)
	end subroutine

end module
