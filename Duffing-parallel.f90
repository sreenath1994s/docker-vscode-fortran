!######### Duffing oscillator program for calculating the safe basin ##########!

program Duffing

    use omp_lib         ! Use Open MP Library for parallel processing

    implicit none

    !########## Constants and Variable ############!
    
    double precision, parameter :: pi =  4.0d0*ATAN(1.0d0)
    
    double precision :: k, beta, gamma               ! damping parameters
    double precision :: alp1, alp3, alp5, alp7, alp9 ! Restoring parameters
    double precision :: B, omega                     ! Excitation parameters 
    double precision :: t_start, t_end, T, Ndim, dt, phi_v                     ! Simulation parameters 
    double precision :: phi_final, phi_dot_final                               !
    double precision, dimension(:), allocatable :: phi_bounds, phi_dot_bounds  !
    double precision, dimension(:), allocatable :: B_list                      !
	logical, dimension(:,:,:), allocatable :: state_list                       !
	double precision, dimension(:,:,:), allocatable :: phi_final_list          !
	double precision, dimension(:,:,:), allocatable :: phi_dot_final_list      !
    integer :: grid_size, B_size, i, j, h                                      !
    logical :: state = .false.                                                 !
       
    double precision, dimension(2) :: y, a
    double precision :: T1, T2
    
    T1 = wtime( )
    
    !######### File ########!
    
    open(unit = 1, file="Safe_List.csv")
    open(unit = 2, file="Safe_Final_List.csv")
    
    write(1,*) "Phi_initial", ",", "Phi_dot_initial"
    write(2,*) "B", ",", "Phi_initial", ",", "Phi_dot_initial", ",", "Phi_final", ",", "Phi_dot_final"
    

    !########## Variables assignment ############!
    
    omega   = 0.31260d0
    k       = 0.00788d0
    beta    = 0.51914d0
    gamma   = 0.0d0
    phi_v   = 1.22173d0
    
    alp1    = 0.09772d0
    alp3    = 0.54952d0
    alp5    = -1.62631d0
    alp7    = 1.31168d0
    alp9    = -0.35322d0
    
    T       = 2.0d0 * pi / omega
    Ndim    = 200.0
    dt      = T/Ndim
	
	t_start = 0.0d0
    t_end   = T*100.0d0  ! Simulation is run for 100 cycles
    
    grid_size = 100
    B_size    = 200
    
    !############# Main Logic ##############!
    
    allocate(B_list(B_size))                                        ! Allocate memory for arrays
	allocate(State_list(grid_size, grid_size, B_size))
	allocate(Phi_final_list(grid_size, grid_size, B_size))
	allocate(Phi_dot_final_list(grid_size, grid_size, B_size))
	
	
    call linspace(from=0.0000d0, to_=0.0450d0, array=B_list)         ! Create arrays for B_list

        
    allocate(phi_bounds(grid_size))           ! Allocate memory for arrays based on the grid size           
    allocate(phi_dot_bounds(grid_size))       !
        
    call linspace(from=-1.5d0,   to_=1.5d0, array=phi_bounds)       ! Create arrays for phi and phi_dot
    call linspace(from=-1.5d0,   to_=1.5d0, array=phi_dot_bounds)   !
	
	! Parallel directives for OpenMP
    !$OMP PARALLEL  
	!$OMP DO
    do h = 1, B_size
        
        print *, h, " B = ", B_list(h) 
        
        do i = 1, grid_size     
            
            do j = 1, grid_size
                
                    state = .true.
                
                    call solve_duffing(phi_bounds(i), phi_dot_bounds(j), B_list(h), phi_final_list(j,i,h), &
					phi_dot_final_list(j,i,h), state_list(j,i,h), t_start, t_end, dt)
						
            end do
            
        end do
        
    end do
    !$OMP END DO
    !$OMP END PARALLEL
	! End parallel directives for OpenMP
	
    T2 = wtime( )  ! Simulation time
	
	print *, "Simulation finished, writing output...."
	
	do h = 1, B_size
        
        B  =  B_list(h)
        
        print *, h, " B = ", B  
        
        do i = 1, grid_size     
            
            do j = 1, grid_size							
                    
                    if (state_list(j,i,h)) then
                    
                        !print *,"Safe for phi = ", phi_bounds(i), " and phi_dot = ", phi_dot_bounds(j), " B = ", B
                        
                        write(1,*) phi_bounds(i), ",", phi_dot_bounds(j)
                        
                        write(2,"(1X, F10.6, A1, F10.6, A1, F10.6, A1, F10.6, A1, F10.6)") B, ",",  phi_bounds(i), ",",&
						& phi_dot_bounds(j), ",", phi_final_list(j,i,h), ",", phi_dot_final_list(j,i,h)
                        
                    else
                        
                        !print *,"Capsize for phi = ", phi_bounds(i), " and phi_dot = ", phi_dot_bounds(j)
                    
                    endif					
            end do
            
        end do
        
    end do
    
    print *,"Calculation time =", (T2-T1)
	print *,"Total time =", ( wtime( )-T1)
    
    !######### Function and subroutine defenitions ##########!
    
    contains

    ! 1) Subroutine duffing
    
    subroutine solve_duffing(phi, phi_dot, B, phi_final, phi_dot_final, state, t_start, t_end, dt)
    
        implicit none
        
        double precision,    intent (in)    :: phi, phi_dot, B, t_start, t_end, dt
        double precision,    intent (out)   :: phi_final, phi_dot_final
        logical, intent (out)   :: state
        
        double precision, dimension(2) :: y, k1, k2, k3, k4
        double precision :: t
        
        y = (/phi,phi_dot/)
        
        t = t_start
        
        do while (t <= t_end)
            
            !call runge
            
            k1 = calc_derivative(y, t, B)
            
            k2 = calc_derivative(y + dt/2.0d0 * k1 , t + dt/2.0d0, B)
            
            k3 = calc_derivative(y + dt/2.0d0 * k2 , t + dt/2.0d0, B)
            
            k4 = calc_derivative(y + dt * k3, t + dt, B)
            
            y = y + 1.0d0/6.0d0 * dt * (k1+ 2.0d0*k2 + 2.0d0*k3 + k4)
            
            if (abs(y(1)) >= phi_v) then
                state = .false.
                exit      
            else  
                state = .true.      
            endif
            
            t = t + dt
            
        end do
        
        phi_final = y(1)
        
        phi_dot_final = y(2)
                
    end subroutine solve_duffing
    
    
    ! 2) Function to calculate derivatives
    
    function calc_derivative(y, t, B) 
    
        implicit none
        
        double precision, dimension(2), intent(in) :: y
        double precision, intent(in) :: t, B
        double precision, dimension(2) :: calc_derivative
        
        calc_derivative(1) = y(2) ! Phi_dot                                             
        calc_derivative(2) = B * COS(omega*t) - k*y(2) - beta*y(2)*ABS(y(2)) - gamma*(y(2)**3) - &
                            ( alp1*y(1) + alp3*y(1)**3 + alp5*y(1)**5 + alp7*y(1)**7 + alp9*y(1)**9 ) ! Phi_dot_dot
        
    end function calc_derivative
	
	! 3) Subroutine to generate grid of points
	
	subroutine linspace(from, to_, array)
    
        implicit none
    
        double precision, intent(in)  :: from, to_
        double precision, intent(out) :: array(:)
        double precision    :: range      
        integer :: n, i
        
        n     = size(array)
        range = to_ - from

        if (n == 0) return

        if (n == 1) then
            
          array(1) = from
          
          return
        
        end if

        do i=1, n
            
            array(i) = from + range * (i - 1) / (n - 1)
            
        end do
    
    end subroutine
	
	! 3) Function to calculate wall clock time
	
	function wtime ( )
	
      implicit none
    
      integer ( kind = 4 ) clock_max
      integer ( kind = 4 ) clock_rate
      integer ( kind = 4 ) clock_reading
      real ( kind = 8 ) wtime
    
      call system_clock ( clock_reading, clock_rate, clock_max )
    
      wtime = real ( clock_reading, kind = 8 ) &
            / real ( clock_rate, kind = 8 )
    
      
    end function wtime
    
end program Duffing
    
    
        
    
    
    
    
    
    
    
    
