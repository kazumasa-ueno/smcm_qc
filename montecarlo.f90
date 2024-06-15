program mc
	implicit none
	
	integer, parameter :: N=50000
	integer :: nt, i, j, ios, X(N)
	real(8) :: sigma(0:3), M(0:3,0:3), R(0:3,0:3), tau(0:3,0:3), tmp(0:3)
	real(8), parameter :: C=0.25d0, D=0.75d0, dt=0.05d0

	tau(:,:) = 0
	tau(0,1) = 1
	tau(1,0) = 5
	tau(1,2) = 1
	tau(0,2) = 2
	tau(2,3) = 3
	tau(2,0) = 5
	tau(3,0) = 5
	R(:,:) = 0.d0
	R(0,1) = Gamma(C)*Gamma(D)/tau(0,1)
	R(0,2) = Gamma(C)*(1-Gamma(D))/tau(0,2)
	R(1,0) = Gamma(D)/tau(1,0)
	R(1,2) = Gamma(C)*(1-Gamma(D))/tau(1,2)
	R(2,0) = (1-Gamma(C))/tau(2,0)
	R(2,3) = 1.d0/tau(2,3)
	R(3,0) = 1.d0/tau(3,0)
	M(:,0) = (/1-R(0,1)*dt-R(0,2)*dt, R(0,1)*dt, R(0,2)*dt, 0.d0/)
	M(:,1) = (/R(1,0)*dt, 1-R(1,0)*dt-R(1,2)*dt, R(1,2)*dt, 0.d0/)
	M(:,2) = (/R(2,0)*dt, 0.d0, 1-R(2,0)*dt-R(2,3)*dt, R(2,3)*dt/)
	M(:,3) = (/R(3,0)*dt, 0.d0, 0.d0, 1-R(3,0)*dt/)

	sigma(:) = (/0.25, 0.25, 0.25, 0.25/)
	X(1:N/4) = 0
	X(N/4+1:N/2) = 1
	X(N/2+1:N/4*3) = 2
	X(N/4*3+1:N) = 3

	open(unit=10, file="mcout50000.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file digout.txt"
	
	do nt = 1, 2000
		do i = 1, N
			X(i) = sample_from_distribution(M(:,X(i)))
		enddo
		do i = 0, 3
			sigma(i) = 1.d0*count(X==i)/N
		enddo
		write(10,*) sigma
	enddo

contains
	real function Gamma(x)
		implicit none
		real(8), intent(in) :: x
		if(x>0) then
			Gamma = 1-exp(-x)
		else
			Gamma = 0
		endif
	end function Gamma

	function sample_from_distribution(p) result(sample)
		implicit none
		real(8), dimension(:), intent(in) :: p
		real(8) :: random_value, cumulative_prob
		integer :: i, sample

		call random_seed()  ! 乱数のシードを設定
		call random_number(random_value)

		cumulative_prob = 0.0d0
		do i = 1, size(p)
			cumulative_prob = cumulative_prob + p(i)
			if (random_value <= cumulative_prob) then
				sample = i - 1  ! 0-indexed
				exit
			end if
		end do
	end function sample_from_distribution

end program mc