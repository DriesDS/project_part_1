module lowrankmod

use matrixconverter

implicit none

contains

	! we berekenen A^t = U,S,Vt
	! we weten dus dat A = Vt^t*S*U^t
	! We moeten dus Vt opslaan in A%Ut en U^t in A%Vt
	subroutine lowrank(rank, epsabs, epsrel)
		use matrixconverter
	
		double precision, optional :: epsrel, epsabs
		integer, optional :: rank
		integer :: m, n, lda, ldu, ldvt, lwork, info, k, s, r
		character :: jobu, jobvt
		double precision, dimension(:), allocatable :: Sigma
		double precision, dimension(:), allocatable :: work, U
		double precision, dimension(1,1) :: dummy
		double precision :: rc
		type(Matrix), pointer :: A,B

		call matrixReader(A)
		m = size(A%Ut,1)
		n = size(A%Ut,2)
		allocate(Sigma(min(m,n)))
		lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))
		allocate(work(lwork))
		info = 0
		
		allocate(B)
		B%full = .false.
		if (m<n) then 
			jobu = 'S' !first min(m,n) left singular vectors returned in A
			jobvt = 'S' !first min(m,n) right singular vectors returned in Vt
			allocate(B%Vt(m,m))
			allocate(B%Ut(m,n))
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, Sigma, B%Vt, m, B%Ut, m, work, lwork, info)
			B%Vt = transpose(B%Vt)
		else
			jobu = 'O' !use O for storing in A
			jobvt = 'S'
			allocate(B%Ut(n,n))
			allocate(B%Vt(n,m))
			call dgesvd(jobu, jobvt, m, n, A%Ut, m, Sigma, A%Ut, m, B%Ut, n, work, lwork, info)
			B%Vt = transpose(A%Ut)
		endif
		
		if (present(rank)) write(*,*) 'rank: ', rank
		if (present(epsrel)) write(*,*) 'epsrel: ', epsrel
		if (present(epsabs)) write(*,*) 'epsabs: ', epsabs
		r = 0
		if (present(epsrel)) then
			do s = 1,min(m,n)
				if (Sigma(s)>Sigma(1)*epsrel) r=r+1
			enddo
		elseif (present(epsabs)) then
			write(*,*) 'min(m,n): ', min(m,n)
			write(*,*) 'r:',r,'Sigma: ',Sigma
			do s = 1,min(m,n)
				if (Sigma(s)>epsabs) r = r+1
			enddo
		elseif (present(rank)) then
			r = rank
			do s = 1,rank
				if (Sigma(s) == 0) then
					r = s
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
				if (Sigma(s) == 0) then
					r = s
					exit
				endif
			enddo
		endif
		
		do s = 1,r
			B%Ut(s,:) = B%Ut(s,:)*Sigma(s)
		enddo
		
		call matrixWriter(B,r)
		
		deallocate(work)
		deallocate(A,B)
		
	end subroutine

end module

