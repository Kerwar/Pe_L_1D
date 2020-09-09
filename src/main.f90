program main
    
    use functions
    use initial
    use variables_m  
    
    implicit none
    
    double precision :: Pe_mean(1000)
    integer :: i,j 
    CALL CPU_TIME(time_start)
    
    call readparam()
    
    allocate(A(2,3*n,3*n), b(2,3*n), x(2,3*n), xp(2,3*n), &
    rho(2,n), omega_b(2,n), omega_c(2,n), xgrid(n))
    ! read input constants file
    
    gamma = q/(q+1d0)
    h = DBLE(xmax-xmin)/DBLE(n-1)
    do i = 1,n
        if(abs(xmin+(i-1)*h-x1change)<=10d-6) then
            nmin = i
            exit
        endif
        if (i ==n) write(*,*) "ERROR: COGE EL NUMERO DE PUNTOS MULTIPLO " &
        // "DE LA LONGITUD DEL DOMINIO"
    enddo
    do i =1,n
        if(abs(xmin+(i-1)*h-x2change)<=10d-6) then
            nmax = i
            exit
        endif
    enddo
    
    index_fixt = floor(0.4*n)
    if(fixt <= xmax) then 
        xgrid = (/(xmin+h*(i-1),i=1,n)/)
        index_fixt = minloc(abs(xgrid-fixt), dim =1)
    endif
    if (initial_data == 0) then
        CALL initialT(x(1,1:n),q,index_fixt)
        CALL initialT(x(2,n:1:-1),q,index_fixt)
        
        CALL initialY(x(1,n+1:2*n),index_fixt) 
        CALL initialY(x(2,2*n:n+1:-1),index_fixt)
        
        CALL initialZ(x(1,2*n+1:3*n),index_fixt) 
        CALL initialZ(x(2,3*n:2*n+1:-1),index_fixt)
    else
        CALL read_data(initial_data)
    endif
    
    
    rho = 1d0/(1d0+(gamma*x(:,1:n)/(1d0-gamma)))
    if (ctedens == 1) rho = 1d0
    xp = x
    write(*,'(A8, 7A25)') 'Iter', 'Error', 'Pe', 'q'               , & 
    'Maximum Temperature 1', 'Maximum Temperature 2', 'Maximum Z 1', &
    'Maximum Z 2'
    print '(A190)', repeat('-',190)
    OPEN(301,file = 'resultstemp.txt')
    do j = 1,n
        WRITE(301,'(9F20.8)') xmin+(j-1)*h, x(1,j), x(1,n+j), x(1,2*n+j), &
        x(2,j), x(2,n+j), x(2,2*n+j)   
    enddo
    CLOSE(301)
    
    Pe_mean = Pe(1)
    do i = 1,maxiter  
        !Pe(1) = 3.3d0
        !Pe(2) = -Pe(1)
        !if ( i <= 11000) then
        !    Pe(1) = 3.75d0
        !    Pe(2) = -Pe(1)
        !endif
        call temperature(A(:,1:n,1:n), b(:,1:n), h, Pe, bK, q, rho, &
        x(:,2*n+1:3*n), x(:,1:n), nmin, nmax)
        CALL onestepgaussseidel(A(1,1:n,1:n), b(1,1:n), x(1,1:n), xp(1,1:n))
        CALL onestepgaussseidel(A(2,1:n,1:n), b(2,1:n), x(2,1:n), xp(2,1:n))
        if(fixt <= xmax .and. fixt>=xmin) call fixingt(x(:,1:n),x(:,2*n+1:3*n)&
        ,q, bK, rho, h, Pe, index_fixt, relaxPe,i,maxiter, q_or_Pe)
        x(:,1:n) = (1d0-relax)*xp(:,1:n) + relax*x(:,1:n)
        !if ( i >= 200000 .and. i <= 300000) then
        !    Pe_mean(mod(i,show)+1) = Pe(1) 
        !    Pe(1) = sum(Pe_mean)/show
        !    Pe(2) = -Pe(1)
        !endif
        if (ctedens == 0) rho = 1d0/(1d0 + gamma*x(:,1:n)/(1d0-gamma))
        
        CALL massfractionF(A(:,n+1:2*n,n+1:2*n), b(:,n+1:2*n) , h, Pe, Le_f, &
        beta, gamma, rho, x(:,2*n+1:3*n), x(:,1:n))
        CALL onestepgaussseidel(A(1,n+1:2*n,n+1:2*n), b(1,n+1:2*n), &
        x(1,n+1:2*n), xp(1,n+1:2*n))
        CALL onestepgaussseidel(A(2,n+1:2*n,n+1:2*n), b(2,n+1:2*n), &
        x(2,n+1:2*n), xp(2,n+1:2*n))
        x(:,n+1:2*n) = (1d0-relax)*xp(:,n+1:2*n) + relax*x(:,n+1:2*n)
        
        CALL massfractionZ(A(:,2*n+1:3*n,2*n+1:3*n), h, Pe, Le_z, beta, gamma, rho, x(:,n+1:2*n), x(:,1:n))
        CALL onestepgaussseidel(A(1,2*n+1:3*n,2*n+1:3*n), b(1,2*n+1:3*n), x(1,2*n+1:3*n), xp(1,2*n+1:3*n))
        CALL onestepgaussseidel(A(2,2*n+1:3*n,2*n+1:3*n), b(2,2*n+1:3*n), x(2,2*n+1:3*n), xp(2,2*n+1:3*n))
        x(:,2*n+1:3*n) = (1d0-relax)*xp(:,2*n+1:3*n) + relax*x(:,2*n+1:3*n)
        
        
        if (mod(i,show) == 0) then
            if (mod(i,show*10) == 0) then
                print '(A190)', repeat('-',190)
                write(*,'(A8, 7A25)') 'Iter', 'Error', 'Pe', 'q'               , & 
                'Maximum Temperature 1', 'Maximum Temperature 2', 'Maximum Z 1', &
                'Maximum Z 2'
                print '(A190)', repeat('-',190)
            endif
            print '(i8, 7F25.13)', i, maxval(x-xp), Pe(1), q, maxval(x(1,1:n)), &
            maxval(x(2,1:n)), maxval(x(1,2*n+1:3*n)), maxval(x(2,2*n+1:3*n)) 
            OPEN(301,file = 'resultstemp.txt')
            do j = 1,n
                WRITE(301,'(9F20.8)') xmin+(j-1)*h, x(1,j), x(1,n+j), x(1,2*n+j), x(2,j), x(2,n+j), x(2,2*n+j)   
            enddo
            CLOSE(301)
        end if
        if (maxval(x-xp) <= tol) exit
        
        xp = x
    enddo
    
    omega_b = rho**2*beta*beta*dexp((beta*(x(:,1:n)-1d0))/(1d0+gamma*(x(:,1:n)-1d0)))*x(:,2*n+1:3*n)*x(:,n+1:2*n)
    omega_c = rho**2*x(:,2*n+1:3*n)
    
    write(file_output,17) n,xmax, x2change, ctedens, Le_f, Le_z, q, beta, gamma, bK, Pe(1) 
    
    
    if (isnan(Pe(1)) .eqv. .false.) then
        OPEN(201,file = file_output,status='replace')
        OPEN(401,file = 'initial_profiles.txt')
        write(201, '(A3,9E22.12)') "#", xmax, x2change, Le_F, Le_Z, q, beta, &
        gamma, bK, Pe(1)
        write(201, '(13A22)') '"X"', '"T_1"', '"F_1"', '"Z_1"', '"rho1"', &
        '"omega_b1"', '"omega_c1"', '"T_2"', '"F_2"', '"Z_2"', '"rho2"', &
        '"omega_b2"', '"omega_c2"'
        do i = 1,n
            WRITE(201,'(13E22.12)') xmin+(i-1)*h, x(1,i), x(1,n+i), x(1,2*n+i), &
            1d0/rho(1,i), omega_b(1,i), omega_c(1,i), x(2,i), x(2,n+i),         &
            x(2,2*n+i), 1d0/rho(2,i), omega_b(2,i), omega_c(2,i)
            WRITE(401,'(9E22.12)') xmin+(i-1)*h, x(1,i), x(1,n+i), x(1,2*n+i), x(2,i), x(2,n+i), x(2,2*n+i)    
        enddo
        CLOSE(201) 
        CLOSE(401)
    endif
    
    
    deallocate(A,b,x,xp,rho,omega_b,omega_c,xgrid)
    CALL CPU_TIME(time_finish)
    
    print *, int(time_finish-time_start)/3600,':' ,mod(int(time_finish-time_start)/60,60), ':',mod(int(time_finish-time_start), 60) 
    print *, Pe
    print*, q
    17 format('./Results/results_n-',I0,'_range','-',F0.0,'-',F0.0,'_cterho-',i1,'_lef-',F0.1,'_lez-',F0.1,'_q-',F0.1,&
    '_beta-',F0.1,'_gamma',F0.1,'_b-',F0.2,'_Pe-',F0.4,'.txt')
end program main
