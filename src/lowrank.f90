module lowrankmod

use matrixconverter
use fullmod

implicit none

contains

	! we berekenen A^t = U,S,Vt
	! we weten dus dat A = Vt^t*S*U^t
	! We moeten dus Vt opslaan in A%Ut en U^t in A%Vt
	
	subroutine lowrank(A, B, rank, epsabs, epsrel)
		double precision, optional :: epsrel, epsabs
		integer, optional :: rank
		integer :: m, n, lwork, info, s, r
		character :: jobu, jobvt
		double precision, dimension(:), allocatable :: Sigma, work
		double precision :: rc
		type(Matrix), pointer, intent(inout) :: A, B
		type(Matrix), pointer :: dump

		if (.not.A%full) then
			call full(A,dump)
			deallocate(A%Ut, A%Vt)
			A%full = .true.
			allocate(A%Ut(size(dump%Ut,1), size(dump%Ut,2)))
			A%Ut = dump%Ut
			deallocate(dump%Ut)
			deallocate(dump)
		endif

		m = size(A%Ut,1)
		n = size(A%Ut,2)
		allocate(Sigma(min(m,n)))
		lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))
		allocate(work(lwork))
		info = 0
		
		allocate(B)
		B%full = .false.
		if (m<n) then 
			jobu = 'S' !first min(m,n) left singular vectors returned in Ut
			jobvt = 'O' !first min(m,n) right singular vectors returned in A
			allocate(A%Vt(m,m))
			A%full = .false.
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, Sigma, A%Vt, m, A%Ut, m, work, lwork, info)
			A%Vt = transpose(A%Vt)
		else
			jobu = 'O' !use O for storing in A
			jobvt = 'S'
			allocate(A%Vt(n,n))
			A%full = .false.
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, Sigma, A%Ut, m, A%Vt, n, work, lwork, info)
		endif

		r = 0
		if (present(epsrel)) then
			do s = 1,min(m,n)
				if (Sigma(s)>Sigma(1)*epsrel) r=r+1
			enddo
		elseif (present(epsabs)) then
			do s = 1,min(m,n)
				if (Sigma(s)>epsabs) r = r+1
			enddo
		elseif (present(rank)) then
			r = rank
			do s = 1,rank
				if (Sigma(s)/Sigma(1) < 1d-14) then
					r = s-1
					exit
				endif
			enddo
		else
			! wanneer er geen argument is meegegeven om de rank van de benadering
			! te bepalen, nemen we ofwel de helft van de kleinste dimensie van de matrix
			! (op deze manier gebruiken we voor dezelfde hoeveelheid geheugen)
			! ofwel enkel alle van nul verschillende waarden.
			rc = min(m,n)/2
			r = floor(rc)
			do s = 1,floor(rc)
				if (Sigma(s)/Sigma(1) < 1d-14) then
					r = s-1
					exit
				endif
			enddo
		endif

		if (m<n) then
			B%Ut => A%Ut(:r, :)
			B%Vt => A%Vt(:r, :)
			B%pointU = .true.
			B%pointV = .true.
		else
			B%Ut => A%Vt(:r, :)
			allocate(B%Vt(r,n))
			B%Vt = transpose(A%Ut(:, :r))
			B%pointU = .true.
			B%pointV = .false.
		endif
		
		do s = 1,r
			B%Ut(s,:) = B%Ut(s,:)*Sigma(s)
		enddo
				
		deallocate(work, Sigma)

	end subroutine

end module
