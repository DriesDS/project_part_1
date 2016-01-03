module makeGHmatmod

	use matrixconverter
	use makegfullmod
	use lowrankmod

implicit none

integer, parameter :: k=3, smalls=8
double precision, parameter :: pi=4d0*atan(1d0)

type HMatrix
	logical :: subtypeH
	type(Matrix), pointer :: endmat
	type(HMatrix), pointer :: leftup
	type(HMatrix), pointer :: rightup
	type(HMatrix), pointer :: leftdown
	type(HMatrix), pointer :: rightdown
end type

contains

	subroutine makeGHmat(AH, N, gamma, elems)
		type(HMatrix), pointer :: AH
		integer, intent(in) :: N
		integer, intent(inout) :: elems
		double precision, intent(in) :: gamma
		double precision :: NN

		NN = N

		if (min(mod(log(NN)/log(2d0),1d0), abs(mod(log(NN)/log(2d0),1d0))) > 1d-5) then
			write(0,*) "N has to be a power of 2"
			return
		elseif (gamma<=0) then
			write(0,*) "gamma has to be real positive greater than zero."
			return
		endif

		AH%subtypeH = .true.
		call makeGrec(AH, N, 1, 1, N, gamma)

		call elemsinHmat(AH, elems)

	end subroutine

	recursive subroutine makeGrec(GH, N, beginx, beginy, s, gamma)
	! preconditie: N is power of 2
		type(HMatrix), pointer :: GH
		type(Matrix), pointer :: fullapprox, rankapprox
		integer, intent(in) :: N, beginx, beginy, s
		double precision, intent(in) :: gamma
		double precision :: ister, jster, Dist

		ister = beginx + s*1d0/2
		jster = beginy + s*1d0/2
		Dist = min(abs(ister-jster),abs(ister-jster-N),abs(ister-jster+N))

		! makeing a rank-k approximation of the matrix when the distance to the
		! diagonal is wide enough
		if (Dist>s*gamma) then
			GH%subtypeH = .false.

			!!! uncomment for approximation by formula
			allocate(GH%endmat)
			GH%endmat%full = .false.
			GH%endmat%pointU = .false.
			GH%endmat%pointV = .false.
			allocate(GH%endmat%Ut(k,s),GH%endmat%Vt(k,s))
			call calculateGapprox(GH%endmat, N, beginx, beginy, s)
			!!! \uncomment for approximation by formula

! 			!!! uncomment for approximation by svd
! 			allocate(fullapprox)
! 			allocate(GH%endmat)
! 			fullapprox%full = .true.
! 			fullapprox%pointU = .false.
! 			GH%endmat%full = .false.
! 			GH%endmat%pointU = .false.
! 			GH%endmat%pointV = .false.
! 			allocate(fullapprox%Ut(s,s))
! 			allocate(GH%endmat%Ut(k,s),GH%endmat%Vt(k,s))

! 			call calculateGt(fullapprox%Ut, N, beginx, beginy, s)
! 			call lowrank(fullapprox, rankapprox, rank=k)

! 			GH%endmat%Ut(1:k, 1:s) = rankapprox%Ut(1:k, 1:s)
! 			GH%endmat%Vt(1:k, 1:s) = rankapprox%vt(1:k, 1:s)
! 			call M_dealloc(rankapprox)
! 			call M_dealloc(fullapprox)
! 			!!! \uncomment for approximation by svd

			return
		endif

		! making a full matrix because we reached the end of recursion due to
		! too small submatrices. Here the matrix is full.
		if (s <= smalls) then
			GH%subtypeH = .false.
			allocate(GH%endmat)
			GH%endmat%full = .true.
			GH%endmat%pointU = .false.
			allocate(GH%endmat%Ut(s, s))
 			call calculateGt(GH%endmat%Ut, N, beginx, beginy, s)
 			return
		endif

		! recursie
		GH%subtypeH = .true.
		allocate(GH%leftup)
		allocate(GH%rightup)
		allocate(GH%leftdown)
		allocate(GH%rightdown)
		call makeGrec(GH%leftup, N, beginx, beginy, s/2, gamma)
		call makeGrec(GH%rightup, N, beginx+s/2, beginy, s/2, gamma)
		call makeGrec(GH%leftdown, N, beginx, beginy+s/2, s/2, gamma)
		call makeGrec(GH%rightdown, N, beginx+s/2, beginy+s/2, s/2, gamma)
 
	end subroutine

	subroutine calculateGapprox(mat, N, beginx, beginy, s)
		type(Matrix), pointer :: mat
		integer, intent(in) :: N, beginx, beginy, s
		integer :: ister, i, j
		double precision :: arg, d

		ister = beginy + s*1d0/2 - 0.5d0

		mat%Ut(1,1:s) = 1d0
		mat%Ut(2,1:s) = (/ (2*pi*(beginy+i - ister)/N, i=0,s-1) /)
		mat%Ut(3,1:s) = (mat%Ut(2,1:s)**2)/2d0

		do j = -1,s-2 ! cause indices of j go from 0 to N-1
			arg = pi*(2d0*ister-2d0*(j+beginx)-1d0)/N
			d = 2d0-2d0*cos((1d0+2d0*(j+beginx)-2d0*ister)*pi/N)
			mat%Vt(1,j+2) = -0.5d0*log(d)
			mat%Vt(2,j+2) = -1d0*d**(-1d0)*sin(arg)
			mat%Vt(3,j+2) = 2d0*d**(-2d0)*sin(arg)**2-d**(-1d0)*cos(arg)
		enddo

	end subroutine

	recursive subroutine elemsinHmat(AH, elems)
		type(HMatrix), pointer, intent(in) :: AH
		integer :: elems, elemssub

		! trivial case
		if (.not. AH%subtypeH) then
			elems = size(AH%endmat%Ut,1) * size(AH%endmat%Ut,2)
			if (.not. AH%endmat%full) then
				elems = elems + size(AH%endmat%Vt,1) * size(AH%endmat%Vt,2)
			endif
			return
		endif

		! recursion
		call elemsinHmat(AH%leftup, elemssub)
		elems = elemssub
		call elemsinHmat(AH%leftdown, elemssub)
		elems = elems+elemssub
		call elemsinHmat(AH%rightup, elemssub)
		elems = elems+elemssub
		call elemsinHmat(AH%rightdown, elemssub)
		elems = elems+elemssub

	end subroutine

	recursive subroutine Hm_dealloc(AH)
		type(HMatrix), pointer :: AH

		if (.not. AH%subtypeH) then
			call M_dealloc(AH%endmat)
		else
			call Hm_dealloc(AH%leftup)
			call Hm_dealloc(AH%rightup)
			call Hm_dealloc(AH%leftdown)
			call Hm_dealloc(AH%rightdown)
		endif

		deallocate(AH)

	end subroutine

end module
