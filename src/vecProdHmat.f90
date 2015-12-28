module vecProdHmatmod

	use matrixconverter
	use makeGHmatmod
	use matprodmod

contains

	subroutine vecProdHmat(AH, N, gamma, vec, prodmat)
		type(HMatrix), pointer :: AH
		type(Matrix), pointer :: vec, prodmat
		double precision, intent(in) :: gamma
		integer, intent(in) :: N
		integer :: flops
		double precision :: NN

		NN = N
		if (min(mod(log(NN)/log(2d0),1d0), abs(mod(log(NN)/log(2d0),1d0))) > 1d-5) then
			write(0,*) "N has to be a power of 2"
			return
		endif

		if (.not. vec%full .or. size(vec%Ut,2) .ne. N) then
			write(0,*) 'Invalid size of matrix.'
			return
		endif

		call makeGHmat(AH, N, gamma)

		call recvecProdHmat(AH, N, vec, prodmat, flops)

		write(0,'(a,i0,a,i0)') "het aantal gebruikte flops = ", flops, ", met een volledige matrix zou dit ", 2*N**2-N, " zijn."

	end subroutine

	recursive subroutine recvecProdHmat(AH, s, vec, prodmat, flops)
		type(HMatrix), pointer :: AH
		type(Matrix), pointer, intent(in) :: vec
		type(Matrix), pointer :: prodmat
		type(Matrix), pointer :: vec1, vec2, prodmatsub
		integer, intent(in) :: s
		integer :: flops, subflops

		! trivial case: AH is a full or rank-k matrix, not an HMatrix
		if (.not. AH%subtypeH) then
			call matprod(AH%endmat, vec, prodmat)
			if (AH%endmat%full) then
				flops = 2*size(AH%endmat%Ut,1)**2-size(AH%endmat%Ut,1)
			else
				flops = 11*size(AH%endmat%Ut,1)-1
			endif
			return
		endif

		! recursie:
		allocate(vec1, vec2, prodmat)
		vec1%full = .true.
		vec2%full = .true.
		prodmat%full = .true.
		vec1%pointU = .false.
		vec2%pointU = .false.
		prodmat%pointU = .false.
		allocate(vec1%Ut(1,s/2), vec2%Ut(1,s/2), prodmat%Ut(1,s))
		vec1%Ut(1,1:s/2) = vec%Ut(1,1:s/2)
		vec2%Ut(1,1:s/2) = vec%Ut(1,s/2+1:s)

		call recvecProdHmat(AH%leftup, s/2, vec1, prodmatsub, subflops)
		prodmat%Ut(1,1:s/2) = prodmatsub%Ut(1,1:s/2)
		call M_dealloc(prodmatsub)
		flops = subflops

		call recvecProdHmat(AH%rightup, s/2, vec2, prodmatsub, subflops)
		prodmat%Ut(1,1:s/2) = prodmat%Ut(1,1:s/2)+prodmatsub%Ut(1,1:s/2)
		call M_dealloc(prodmatsub)
		flops = flops+subflops+s/2

		call recvecProdHmat(AH%leftdown, s/2, vec1, prodmatsub, subflops)
		prodmat%Ut(1,s/2+1:s) = prodmatsub%Ut(1,1:s/2)
		call M_dealloc(prodmatsub)
		flops = flops+subflops

		call recvecProdHmat(AH%rightdown, s/2, vec2, prodmatsub, subflops)
		prodmat%Ut(1,s/2+1:s) = prodmat%Ut(1,s/2+1:s)+prodmatsub%Ut(1,1:s/2)
		call M_dealloc(prodmatsub)
		flops = flops+subflops+s/2

		call M_dealloc(vec1)
		call M_dealloc(vec2)

	end subroutine

end module
