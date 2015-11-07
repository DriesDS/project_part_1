program main

implicit none

call test2()

write(*,'(/,a,/)') 'test with bigger matrices: '

call test()

contains

	subroutine test()
		integer :: i
		double precision :: alpha, beta
		integer :: m,n,k
		double precision, dimension(5,3) :: A = reshape((/&
			1.0, 1.0, 2.0, &
			2.0, 1.0, 3.0, &
			3.0, 1.0, 4.0, &
			4.0, 1.0, 0.0, &
			5.0, 1.0, 0.0  /),(/5,3/),order=(/2,1/))
		double precision, dimension(5,3) :: B = reshape((/&
			5.0, 0.0, 1.0, &
			4.0, 1.0, 0.0, &
			3.0, 0.0, 1.0, &
			2.0, 1.0, 0.0, &
			1.0, 0.0, 1.0  /),(/5,3/),order=(/2,1/))
		double precision, dimension(5,5) :: C

		do i = 1,5
			write(*,'(2(3(e12.3),5x))') A(i,:), B(i,:)
		enddo

		alpha = 1.0
		beta = 0.0
		m = 5
		n = 5
		k = 3
		call dgemm('N','T',m,n,k,alpha,A,m,B,m,beta,C,m)

		write(*,*) '\n', '\n'
		do i = 1,5
			write(*,'(5(e12.3))') C(i,:)
		enddo

	end subroutine

	subroutine test2()
		integer :: i
		double precision,dimension(2,2) :: A = reshape((/1.0, 2.0, 2.0, 3.0/),(/2,2/))
		double precision,dimension(2,2) :: B = reshape((/2.0, 1.0, 4.0, 1.0/),(/2,2/))
		double precision,dimension(2,2) :: C = reshape((/0.0, 0.0, 0.0, 0.0/),(/2,2/))

		do i = 1,2
			write(*,'(2(2(e12.3),5x))') A(i,:), B(i,:)
		enddo

		call dgemm('N','N',2,2,2,1.0,A,2,B,2,0,C,2)

		do i = 1,2
			write(*,'(2(e12.3))') C(i,:)
		enddo

	end subroutine

end program
