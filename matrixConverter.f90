module matrixconverter

implicit none

integer, parameter :: defaultin=5, defaultout=6
integer, parameter :: wp = kind(0.d0)

type Matrix
	! when the matrix is full, only Ut will be used
	logical :: full
	real(wp), allocatable :: Ut(:,:), Vt(:,:)
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
		real(wp) :: rowelems(cols)
		integer :: row
		
		allocate(matp)
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
		matp%full = .false.
		allocate(matp%Ut(rank, rows))
		allocate(matp%Vt(rank, cols))
		do row = 1,rows
			read(inunit,*) elems
			matp%Ut(:,row) = elems
		enddo
		read(*,*)
		do col = 1,cols
			read(inunit,*) elems
			matp%Vt(:,col) = elems
		enddo 
	end subroutine

	subroutine matrixWriter(matp, optoutunit)
		type(Matrix), pointer :: matp
		integer, optional, intent(in) :: optoutunit
		integer :: row, col, rows, cols, rank, outunit
		character(len=16) :: rowstr, colstr, rankstr
		character(len=100) :: header
		
		outunit = defaultout
		if (present(optoutunit)) outunit = optoutunit
		
		if (matp%full) then
			rows = size(matp%Ut,2)
			cols = size(matp%Ut,1)
			write(header,'(a,i0,a,i0,a)') 'full matrix [', rows, ' x ', cols, ']'
		else
			rows = size(matp%Ut,2)
			cols = size(matp%Vt,2)
			rank = size(matp%Ut,1)
			write(header,'(a,i0,a,i0,a,i0,a)') 'rank-', rank, ' matrix [', rows, ' x ', cols, ']'
		endif

		write(outunit,'(a)') header
		do row = 1,rows
			write(outunit,*) matp%Ut(:,row)
		enddo
		if (.not. matp%full) then
			write(outunit,*) '----------'
			do row = 1,cols
				write(outunit,*) matp%Vt(:,row)
			enddo
		endif
	end subroutine

end module

