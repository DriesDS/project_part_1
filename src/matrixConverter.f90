module matrixconverter

implicit none

integer, parameter :: defaultin=5, defaultout=6
character(len=16), parameter :: separationstr='---------'
integer, parameter :: wp = kind(0.d0)

type Matrix
	! when the matrix is full, only Ut will be used
	logical :: full
	logical :: pointU, pointV
	double precision, pointer :: Ut(:,:), Vt(:,:)
end type

contains

	subroutine matrixReader(matp, optinunit)
		type(Matrix), pointer :: matp
		integer, optional, intent(in) :: optinunit
		integer :: inunit, rows, cols, rank
		character(len=8) :: matkind, mat, r, x, c
		
		inunit = defaultin
		if (present(optinunit)) inunit = optinunit
		
		read(inunit,*) matkind, mat, r, x, c
		read(r(2:),*) rows
		read(c(:len_trim(c)-1),*) cols
		
		if (index(matkind,'rank') == 1) then
			read(matkind(6:),*) rank
			call rankMatrixReader(matp, inunit, rows, cols, rank)
		elseif(index(matkind,'full') == 1) then
			call fullMatrixReader(matp, inunit, rows, cols)
		endif

	end subroutine

	subroutine fullMatrixReader(matp, inunit, rows, cols)
		type(Matrix), pointer :: matp
		integer, intent(in) :: inunit, rows, cols
		double precision :: rowelems(cols)
		integer :: row
		
		allocate(matp)
		matp%pointU = .false.
		matp%pointV = .false.
		matp%full = .true.
		allocate(matp%Ut(cols, rows))
		do row = 1,rows
			read(inunit,*) rowelems
			matp%Ut(:,row) = rowelems
		enddo

	end subroutine

	subroutine rankMatrixReader(matp, inunit, rows, cols, rank)
		type(Matrix), pointer :: matp
		integer, intent(in) :: inunit, rows, cols, rank
		real(wp) :: elems(rank)
		integer :: row, col
		
		allocate(matp)
		matp%pointU = .false.
		matp%pointV = .false.
		matp%full = .false.
		allocate(matp%Ut(rank, rows))
		allocate(matp%Vt(rank, cols))
		do row = 1,rows
			read(inunit,*) elems
			matp%Ut(:,row) = elems
		enddo
		read(inunit,*)
		do col = 1,cols
			read(inunit,*) elems
			matp%Vt(:,col) = elems
		enddo 
	end subroutine

	subroutine matrixWriter(matpu, rank, optoutunit)
		type(Matrix), pointer :: matpu
		integer, optional, intent(in) :: rank, optoutunit
		integer :: row, col, rows, cols, rank_k, outunit, dig, prec
		character(len=32) :: header
		character(len=32) :: fmtstr

		dig = 24
		prec = 16
		outunit = defaultout
		if (present(optoutunit)) then
			outunit = optoutunit
		endif

		if (matpu%full) then
			rows = size(matpu%Ut,2)
			cols = size(matpu%Ut,1)
			rank_k = cols
			write(header,'(a,i0,a,i0,a)') 'full matrix [', rows, ' x ', cols, ']'
		else
			rows = size(matpu%Ut,2)
			cols = size(matpu%Vt,2)
			rank_k = size(matpu%Ut,1)
			if(present(rank)) rank_k = rank
			write(header,'(a,i0,a,i0,a,i0,a)') 'rank-', rank_k, ' matrix [', rows, ' x ', cols, ']'
		endif

		write(outunit,'(a)') trim(header)
		write(fmtstr,'(a,i0,a,i0,a,i0,a)') '(', rank_k, '(e', dig, '.', prec, '))'
		do row = 1,rows
			write(outunit,fmtstr) matpu%Ut(:rank_k,row)
		enddo
		if (.not. matpu%full) then
			write(outunit,*) '----------'
			do row = 1,cols
				write(outunit,fmtstr) matpu%Vt(:rank_k,row)
			enddo
		endif
		write(outunit,*)
	end subroutine

	subroutine lowlevelwriter(U, jobU, V, jobV, rank, optoutunit)
		double precision, intent(in) :: U(:,:)
		double precision, optional, intent(in) :: V(:,:)
		character :: jobU
		character, optional :: jobV
		character(len=32) :: header
		integer, optional, intent(in) :: optoutunit, rank
		integer :: outunit, rows, cols, row, col, r
		logical :: full
		
		full = .true.
		if (present(V) .and. present(jobV)) full = .false.
		outunit = defaultout
		if (present(optoutunit)) outunit = optoutunit
		
		if (full) then
			if (jobU == 'T' .or. jobU == 't') then
				rows = size(U,2)
				cols = size(U,1)
			else
				rows = size(U,1)
				cols = size(U,2)
			endif
			write(header,'(a,i0,a,i0,a)') 'full matrix [', rows, ' x ', cols, ']'
		else
			if (jobU == 'T' .or. jobV == 't') then
				rows = size(U,2)
				r = size(U,1)
			else
				rows = size(U,1)
				r = size(U,2)
			endif
			if (present(rank)) r = rank
			if (jobV == 'T' .or. jobV == 't') then
				cols = size(V,1)
			else
				cols = size(V,2)
			endif

			write(header,'(a,i0,a,i0,a,i0,a)') 'rank-', r, ' matrix [', rows, ' x ', cols, ']'
		endif

		write(outunit,'(a)') trim(header)
		
		if (jobU == 'T' .or. jobU == 't') then
			do row = 1,rows
				write(outunit,*) U(:r,row)
			enddo
		else
			do row = 1,rows
				write(outunit,*) U(row,:r)
			enddo
		endif
		if (.not. full) then
			write(outunit,*) separationstr
			if (jobV == 'T' .or. jobV == 't') then
				do col = 1,cols
					write(outunit,*) V(:r,col)
				enddo
			else
				do col = 1,cols
					write(outunit,*) V(col,:r)
				enddo
			endif
		endif
		write(outunit,*)

	end subroutine

end module

