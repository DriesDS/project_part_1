module testmod

use matrixconverter
use fullmod
use lowrankmod
use matprodmod
use makeGFullmod
use solveIntFullmod
use plotfieldmod
use vecProdHmatmod
use makeGHmatmod

implicit none

contains

	subroutine test()
		call test1()
	!	call test2()
	!	call test3()
		call test4()
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
		call loadmatrices(M1, M2, M3, M4, xfield)
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
		allocate(X%Ut(1,2))
		X%Ut(1,1) = 0.9
		X%Ut(1,2) = -0.7

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
		integer, parameter :: N=20, Nfull=50, stepprod=50, steplowrank=20, stepfull=5
		integer :: i
		type(Matrix), pointer :: A, B, C, D, E, F
		double precision :: start, gestopt, timef, timep(2,N), timel(N), timefull(N)

		call random_seed()
		allocate(A)
		A%full = .true.
		A%pointU = .false.
		allocate(A%Ut(stepprod*N,stepprod*N))
		call random_number(A%Ut)
		allocate(B)
		B%full = .true.
		B%pointU = .false.
		allocate(B%Ut(stepprod*N,stepprod*N))
		call random_number(B%Ut)
		
		call cpu_time(start)
		call matprod(A,B,C)
		call cpu_time(gestopt)
		timef = gestopt-start

		write(*,'(a)') '     rank:          rank*rank:          rank*full:'
		
		! testen voor tijdsmeting matprod
		do i = 1,N

			allocate(D)
			D%full = .false.
			D%pointU = .false.
			D%pointV = .false.
			allocate(D%Ut(stepprod*i,stepprod*N))
			allocate(D%Vt(stepprod*i,stepprod*N))
			call random_number(D%Ut)
			call random_number(D%Vt)
			allocate(E)
			E%full = .false.
			E%pointU = .false.
			E%pointV = .false.
			allocate(E%Ut(stepprod*i,stepprod*N))
			allocate(E%Vt(stepprod*i,stepprod*N))
			call random_number(E%Ut)
			call random_number(E%Vt)
			call cpu_time(start)
			call matprod(D,E,F)
			call cpu_time(gestopt)
			timep(i,1) = gestopt-start

			call cpu_time(start)
			call matprod(A,D,F)
			call cpu_time(gestopt)

			timep(i,2) = gestopt-start
			write(*,'(i10,3(e20.10))') stepprod*i, timep(i,1), timep(i,2), timef

 		enddo

 		! testen voor tijdsmeting lowrank
 		do i=1,N

			allocate(A)
			A%full = .true.
			A%pointU = .false.
			allocate(A%Ut(steplowrank*N,steplowrank*N))
			call random_number(A%Ut)
 			call cpu_time(start)
 			call lowrank(A,B,steplowrank*i)
 			call cpu_time(gestopt)
 			timel(i) = gestopt-start
 			write(*,*) steplowrank*i, timel(i) 			

 		enddo

 		! testen voor tijdmeting full
 		do i = 1,Nfull
			allocate(D)
			D%full = .false.
			D%pointU = .false.
			D%pointV = .false.
			allocate(D%Ut(stepfull*i,stepfull*Nfull))
			allocate(D%Vt(stepfull*i,stepfull*Nfull))
			call random_number(D%Ut)
			call random_number(D%Vt)
			call cpu_time(start)
			call full(D,E)
			call cpu_time(gestopt)
			timefull(i) = gestopt-start
 			write(*,*) stepfull*i, timefull(i) 			
		enddo

	end subroutine

	subroutine test4()
		type(HMatrix), pointer :: AH
		type(Matrix), pointer :: x1,x2
		integer, parameter :: begin=5, eind=14, beginy=1, eindy=10, gammatestN=2**10
		integer :: elems, i, j, elemslist(eind-begin+1), elemslisty(eindy-beginy+1)
		integer :: nbflops(7,5), curnbflops
		double precision, dimension(5), parameter :: gamma= (/ 1d0, 2d0, 5d0, 1d1, 2d1 /)
		double precision :: y, procent(eind-begin+1), procenty(eindy-beginy+1)
		double precision :: xnorm(1), fnorm, diff(512)

		do i = 1,7
			do j = 1,5
				call flops_amount(2**(i+4), gamma(j), curnbflops)
				nbflops(i,j) = curnbflops
			enddo
		enddo
		do i = 1,7
			write(0,*) nbflops(i,:)
		enddo

! 		call SYSTEM('./hmatrices makeGFull 512 >G.out')
! 		call SYSTEM('cat tests/randn512.in | ./hmatrices vecProdHmat 512 5 >x1.out')
! 		call SYSTEM('cat G.out tests/randn512.in | ./hmatrices matprod >x2.out')

! 		open(10,file='x1.out')
! 		call matrixReader(x1, 10)
! 		close(10)
! 		open(10,file='x2.out')
! 		call matrixReader(x2, 10)
! 		close(10)
! 		diff = (x1%Ut(1,1:512)-x2%Ut(1,1:512))**2
! 		fnorm = sum(diff)
! 		xnorm = sum(x1%Ut,2)
! 		write(*,*) xnorm(1), fnorm

! 		do i = begin, eind
! 			allocate(AH)
! 			call makeGHmat(AH, 2**i, gamma)
! 			call elemsinHmat(AH, elems)
! 			elemslist(1+i-begin) = elems
! 			procent(1+i-begin) = elems*1d0/4**i
! 			call HM_dealloc(AH)
! 		enddo

! 		do i = beginy, eindy
! 			y = i
! 			allocate(AH)
! 			call makeGHmat(AH, gammatestN, y)
! 			call elemsinHmat(AH, elems)
! 			elemslisty(1+i-beginy) = elems
! 			procenty(1+i-beginy) = elems*1d0/(gammatestN**2)
! 			call HM_dealloc(AH)
! 		enddo

! 		do i = 1,eind-begin+1
! 			write(*,'(i10, e12.4)') elemslist(i), procent(i)
! 		enddo
! 		do i = 1,eindy-beginy+1
! 			write(*,'(i10, e12.4)') elemslisty(i), procenty(i)
! 		enddo

	end subroutine

	subroutine flops_amount(N,y, nbfloats)
		integer, intent(in) :: N
		double precision, intent(in) :: y
		integer :: nbfloats
		character(len=128) :: command
		character(len=16) :: nbfloatsstr

		write(command,'(a,i0,x,e12.4,a)') './hmatrices makeGFull ', N, y, ' >flops.out'
		call SYSTEM(command)
		open(10,file='flops.out')
		read(10,*) nbfloatsstr
		close(10)
		read(nbfloats,*) nbfloats

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
