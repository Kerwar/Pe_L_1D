module variables_m
   
    implicit none

    double precision, ALLOCATABLE, dimension(:,:,:), SAVE :: A 
    double precision, ALLOCATABLE, dimension(:,:),  SAVE :: b , x , xp, &
    rho, omega_c, omega_b
    double precision, allocatable, dimension(:), save :: xgrid
    double precision, SAVE :: Le_f, Le_z, beta, gamma, q, bK, relax, &
    tol, xmin, xmax, h, x1change, x2change, time_start, time_finish, fixt,&
    relaxPe
    DOUBLE PRECISION, DIMENSION(2), SAVE ::  Pe 
    INTEGER, SAVE :: n, show, maxiter, nmin, nmax, initial_data, & 
    ctedens,index_fixt, q_or_Pe 
    CHARACTER(400), SAVE :: file_output = ''

    
end module variables_m