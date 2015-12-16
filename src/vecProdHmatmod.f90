module vecProdHmatmod

	use matrixconverter
	use makeGHmatmod

contains

	subroutine vecProdHmat(AH, N, gamma, vec)
		type(HMatrix), pointer :: AH
		type(Matrix), pointer :: vec, prodmat
		double precision, intent(in) :: gamma
		integer, intent(in) :: N

		NN = N
		if (mod(log(NN)/log(2d0),1d0) == 1) then
			write(0,*) "N has to be a power of 2"
			return
		endif

		call MatrixReader(vec)
		allocate(AH)
		call makeGHmat(AH, N, gamma)

		allocate(prodmat)
		prodmat%full = .true.
		prodmat%pointU = .false.
		allocate(prodmat%Ut(1,N))


	end subroutine

	recursive subroutine recvecProdHmat(AH, s, prodmat)

	end subroutine

end module