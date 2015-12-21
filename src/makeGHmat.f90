module makeGHmatmod

	use matrixconverter
	use makegfullmod

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

	subroutine makeGHmat(AH, N, gamma)
		type(HMatrix), pointer :: AH
		integer, intent(in) :: N
		integer :: elems
		double precision, intent(in) :: gamma
		double precision :: NN, delems

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
		delems = elems
		write(0,'(a,i0)') 'number of elements in the matrix: ', elems
		write(0,'(a,e10.3,a)') 'This is a percentage of ', delems/N/N*100, ' of the elements of the matrix.'

	end subroutine

	recursive subroutine makeGrec(GH, N, beginx, beginy, s, gamma)
	! preconditie: N is power of 2
		type(HMatrix), pointer :: GH
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
			allocate(GH%endmat)
			GH%endmat%full = .false.
			GH%endmat%pointU = .false.
			GH%endmat%pointV = .false.
			allocate(GH%endmat%Ut(k,s),GH%endmat%Vt(k,s))
			call calculateGapprox(GH%endmat, N, beginx, beginy, s)
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

		do i = 1,s
			mat%Ut(1,i) = 1d0
			mat%Ut(2,i) = 2*pi*(beginy -1 + i - ister)/N
			mat%Ut(3,i) = mat%Ut(2,i)**2 / 2
		enddo

		do j = 1,s
			arg = pi*(2d0*ister-2d0*j-3d0)/N
			d = 2*(1-cos((+1d0+2d0*j-2d0*ister)*pi/N))
			mat%Vt(1,j) = -0.5d0*log(d)
			mat%Vt(2,j) = d**(-3)*sin(arg)
			mat%Vt(3,j) = 3d0*d**(-5)*sin(arg)**2-d**(-3)*cos(arg)
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
