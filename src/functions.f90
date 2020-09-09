module functions

  use variables_m, only: sparse_t
  implicit none
   
  contains
  
  subroutine temperature(A, b, h, Pe, bK, q, rho, Z, T, nmin, nmax)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: rho, Z, T
    DOUBLE PRECISION, INTENT(IN) :: h, bK, q
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: Pe
    INTEGER, INTENT(IN) :: nmin, nmax
    type(sparse_t), INTENT(INOUT), DIMENSION(:,:)  :: A
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:) ::  b
    INTEGER :: i,n
    
    n = size(b, dim = 2)

    A(:,:)%l = 0.0
    A(:,:)%c = 0.0
    A(:,:)%r = 0.0

    A(1,1)%c = h*h
    
    A(2,1)%c = 2d0 
    A(2,1)%r = -2d0
    
     
    do i = 2, nmin-1
      A(:,i)%l = -Pe(:)*h/2d0-1d0
      A(:,i)%c = 2d0 
      A(:,i)%r = Pe(:)*h/2d0-1d0
  
    enddo
    do i = nmin,nmax
        A(:,i)%l = -Pe(:)*h/2d0-1d0
        A(:,i)%c = 2d0 + bK*h*h
        A(:,i)%r = Pe(:)*h/2d0-1d0
  
    enddo
    do i = nmax+1,n-1
      A(:,i)%l = -Pe(:)*h/2d0-1d0
      A(:,i)%c = 2d0 
      A(:,i)%r = Pe(:)*h/2d0-1d0
  
    enddo
  
    A(1,n)%l = -2d0
    A(1,n)%c = 2d0 
    
    A(2,n)%c = h*h
    
    A%l = A%l /(h*h)
    A%c = A%c /(h*h)
    A%r = A%r /(h*h)
  
    b = q*Z*rho*rho
    b(:,nmin:nmax) = b(:,nmin:nmax) +bK*T(2:1:-1,nmin:nmax)
    
    
    !DIRICHLET BOUNDARY CONDITIONS
    b(1,1) = 0d0
    b(2,n) = 0d0
  
    
  end subroutine temperature
  
  subroutine massfractionF(A, b, h, Pe, Le_f, beta, gamma, rho, Z, T)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: rho, Z, T
    DOUBLE PRECISION, INTENT(IN) :: h,  Le_f, beta, gamma
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: Pe
    type(sparse_t), INTENT(INOUT), DIMENSION(:,:)  :: A
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:) :: b
    INTEGER :: i,n
    
    n = size(b, dim = 2)
  
    A(:,:)%l = 0.0
    A(:,:)%c = 0.0
    A(:,:)%r = 0.0

    A(1,1)%c = h*h*Le_f
    
    A(2,1)%c = 2d0 + rho(2,1)*rho(2,1) * beta*beta * &
    dexp((beta*(T(2,1)-1d0))/(1d0+gamma*(T(2,1)-1d0)))*Z(2,1)*h*h*Le_f
    A(2,1)%r = -2d0
  
    do i = 2,n-1
        A(:,i)%l = -Pe(:)*h*Le_f/2-1d0
        A(:,i)%c = 2d0 + rho(:,i)**2 * beta*beta * &
        dexp((beta*(T(:,i)-1d0))/(1d0+gamma*(T(:,i)-1d0)))*Z(:,i)*h*h*Le_f
        A(:,i)%r = Pe(:)*h*Le_f/2-1d0
        
    enddo
  
    A(1,n)%l = -2d0
    A(1,n)%c = 2d0 + rho(1,n)*rho(1,n) * beta*beta * &
    dexp((beta*(T(1,n)-1d0))/(1d0+gamma*(T(1,n)-1d0)))*Z(1,n)*h*h*Le_f
    
    A(2,n)%c = h*h*Le_f
    
    A(:,:)%l = A(:,:)%l / (h*h*Le_f)
    A(:,:)%c = A(:,:)%c / (h*h*Le_f)
    A(:,:)%r = A(:,:)%r / (h*h*Le_f)

    b(1,1) = 1d0
    b(2,n) = 1d0
  end subroutine massfractionF
  
  subroutine massfractionZ(A, h, Pe, Le_z, beta, gamma, rho, F, T)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: rho, F, T
    DOUBLE PRECISION, INTENT(IN) :: h, Le_z, beta, gamma
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: Pe
    type(sparse_t)  , INTENT(INOUT), DIMENSION(:,:)  :: A
    INTEGER :: i,n
    
    n = size(F, dim = 2)
  
    A(:,:)%l = 0.0
    A(:,:)%c = 0.0
    A(:,:)%r = 0.0

    A(1,1)%c = h*h*Le_z
     
    A(2,1)%c = 2d0 + rho(2,1)*rho(2,1)*h*h*Le_z - rho(2,1)*rho(2,1) * &
    beta*beta*dexp((beta*(T(2,1)-1d0))/(1d0+gamma*(T(2,1)-1d0)))*F(2,1)*h*h*Le_z
    A(2,1)%r = -2d0
    
    do i = 2,n-1
        A(:,i)%l = -Pe(:)*h*Le_z/2d0-1d0
        A(:,i)%c = 2d0 +rho(:,i)*rho(:,i)*h*h*Le_z - rho(:,i)*rho(:,i) * &
        beta*beta*dexp((beta*(T(:,i)-1d0))/(1d0+gamma*(T(:,i)-1d0))) * &
        F(:,i)*h*h*Le_z
        A(:,i)%r = Pe(:)*h*Le_z/2d0-1d0
  
    enddo
    
    A(1,n)%l = -2d0
    A(1,n)%c = 2d0 + rho(1,n)*rho(1,n)*h*h*Le_z - rho(1,n)*rho(1,n) * &
    beta*beta*dexp((beta*(T(1,n)-1d0))/(1d0+gamma*(T(1,n)-1d0)))*F(1,n)*h*h*Le_z
    
    A(2,n)%c = h*h*Le_z
    
    A(:,:)%l = A(:,:)%l / (h*h*Le_z)
    A(:,:)%c = A(:,:)%c / (h*h*Le_z)
    A(:,:)%r = A(:,:)%r / (h*h*Le_z)
  
  end subroutine massfractionZ
  
  subroutine onestepgaussseidel(A, b, x, xp)
    type(sparse_t), INTENT(IN), DIMENSION(:)  :: A
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) ::  b
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) ::  x, xp
    INTEGER :: i,n
    
    n = size(b)
  
    x(1) = (1d0/A(1)%c)*(-A(1)%r*xp(2)+b(1))
    do i= 2, n-1
      x(i) = (1d0/A(i)%c)*&
      (-A(i)%l * x(i-1) -A(i)%r * xp(i+1)+b(i))
    enddo
  
    x(n) = (1d0/A(n)%c)*(-A(n)%l * x(n-1)+b(n))
  
  end subroutine onestepgaussseidel
  
  subroutine fixingt(T, Z, q, bK, rho, h, Pe, index_fixt, relaxPe,i, maxiter, q_or_Pe)
    double precision, intent(in)  :: Z(:,:), bK, rho(:,:), h, relaxPe
    integer, intent(in) :: index_fixt,i,maxiter, q_or_Pe 
    double precision, intent(inout) :: T(:,:), Pe(:),q
    double precision :: t_x,t_xx, reaction_heat, heat_ex
    integer :: n
     
    n = size(T, dim = 2)
    T(1, index_fixt) = 1d0
  
    t_x = (T(1,index_fixt+1)-T(1,index_fixt-1))/(2d0*h)
    t_xx = (T(1,index_fixt+1)-2d0*T(1, index_fixt) &
    +T(1,index_fixt-1)) /(h*h)
    heat_ex = bK*(T(2,index_fixt)-T(1, index_fixt))
    reaction_heat = rho(1,index_fixt)*rho(1,index_fixt)* Z(1,index_fixt)
   
    if (q_or_Pe == 1) then
      Pe(1) = Pe(1)*(1d0-relaxPe)+&
      (relaxPe)*1d0/t_x*(t_xx + q*reaction_heat + heat_ex )
  
      Pe(2) = -Pe(1)
    else
      q = q*(1d0-relaxPe*i/maxiter)+(relaxPe*i/maxiter-0.0)*1d0/reaction_heat*(Pe(1)*t_x-t_xx)
    endif
       
  end subroutine fixingt
 
  !subroutine writetemp()
  !  print '(i7, 7F25.13)', i, maxval(x-xp), Pe(1), q, maxval(x(1,1:n)), &
  !  maxval(x(2,1:n)), maxval(x(1,2*n+1:3*n)), maxval(x(2,2*n+1:3*n)) 
  !  OPEN(301,file = 'resultstemp.txt')
  !  write(301, '(7A20)') "X", "T_1", "F_1", "Z_1", "T_2", "F_2", "Z_2"
  !  do j = 1,n
  !      WRITE(301,'(7F20.8)') xmin+(j-1)*h, x(1,j), x(1,n+j), x(1,2*n+j),&
  !       x(2,j), x(2,n+j), x(2,2*n+j)   
  !  enddo
  !  CLOSE(301)
  !end subroutine writetemp
  end module functions