module functions

  use type_m, only: DP
  use param_m, only: param_t
  use variables_m, only: sparse_t
  implicit none
   
  contains
  
  subroutine set_t(A, b, rho, Z, T, param)
    type(sparse_t), intent(inout), dimension(:,:) :: A
    real(DP)      , intent(inout), dimension(:,:) :: b
    real(DP)      , intent(in   ), dimension(:,:) :: rho, Z, T
    type(param_t) , intent(in   ) :: param
    integer :: i,n
    
    n = size(b, dim = 2)

    A(:,:)%l = 0.0_DP
    A(:,:)%c = 0.0_DP
    A(:,:)%r = 0.0_DP

    A(1,1)%c = 1.0_DP
    
    A(2,1)%c =  2.0_DP * param%h2  
    A(2,1)%r = -2.0_DP * param%h2
    
     
    do i = 2, param%nmin-1
      A(:,i)%l = -param%Pe(:) * param%h1 * 0.5_DP - param%h2
      A(:,i)%c = 2.0_DP * param%h2 
      A(:,i)%r =  param%Pe(:) * param%h1 * 0.5_DP - param%h2
  
    enddo
    do i = param%nmin, param%nmax
        A(:,i)%l = -param%Pe(:) * param%h1 * 0.5_DP - param%h2
        A(:,i)%c = 2.0_DP * param%h2 + param%bK
        A(:,i)%r =  param%Pe(:) * param%h1 * 0.5_DP - param%h2
  
    enddo
    do i = param%nmax + 1, n-1
      A(:,i)%l = -param%Pe(:) * param%h1 * 0.5_DP - param%h2
      A(:,i)%c = 2.0_DP * param%h2 
      A(:,i)%r =  param%Pe(:) * param%h1 * 0.5_DP - param%h2
  
    enddo
  
    A(1,n)%l = -2.0_DP * param%h2
    A(1,n)%c =  2.0_DP * param%h2 
    
    A(2,n)%c = 1.0_DP
    
    b = param%q * Z * rho * rho
    b(:, param%nmin : param%nmax) = b(:, param%nmin : param%nmax) + &
    param%bK * T(2:1:-1, param%nmin : param%nmax)
    
    
    !DIRICHLET BOUNDARY CONDITIONS
    b(1,1) = 0.0_DP
    b(2,n) = 0.0_DP
  
    
  end subroutine set_t
  
  subroutine act_t(A, b, rho, Z, T, param)
    type(sparse_t), intent(inout), dimension(:,:) :: A
    real(DP)      , intent(inout), dimension(:,:) :: b
    real(DP)      , intent(in   ), dimension(:,:) :: rho, Z, T
    type(param_t) , intent(in   ) :: param
    integer :: n
    
    n = size(b, dim = 2)
  
    b = param%q * Z * rho * rho
    b(:, param%nmin : param%nmax) = b(:, param%nmin : param%nmax) + &
    param%bK * T(2:1:-1, param%nmin : param%nmax)
        
    !DIRICHLET BOUNDARY CONDITIONS
    b(1,1) = 0.0_DP
    b(2,n) = 0.0_DP
  
  end subroutine act_t

  subroutine set_F(A, b, rho, Z, T, param)
    type(sparse_t), intent(inout), dimension(:,:) :: A
    real(DP)      , intent(inout), dimension(:,:) :: b
    real(DP)      , intent(in   ), dimension(:,:) :: rho, Z, T
    type(param_t) , intent(in   ) :: param
    INTEGER :: i,n
    
    n = size(b, dim = 2)
  
    A(:,:)%l = 0.0_DP
    A(:,:)%c = 0.0_DP
    A(:,:)%r = 0.0_DP

    A(1,1)%c = 1.0_DP
    
    A(2,1)%c =  2.0_DP * param%h2 / param%Le_f + rho(2,1)** 2 * param%beta** 2 * &
    dexp((param%beta*(T(2,1)-1.0_DP))/(1.0_DP+param%gamma*(T(2,1) - 1.0_DP))) * Z(2,1)
    A(2,1)%r = -2.0_DP * param%h2 / param%Le_f
  
    do i = 2,n-1
        A(:,i)%l = -param%Pe(:) * param%h1 * 0.5_DP - param%h2 / param%Le_f
        A(:,i)%c = 2.0_DP * param%h2 / param%Le_f + rho(:,i)**2 * param%beta*param%beta * &
        dexp((param%beta*(T(:,i)-1.0_DP))/(1.0_DP+param%gamma*(T(:,i)-1.0_DP)))*Z(:,i)
        A(:,i)%r =  param%Pe(:) * param%h1 * 0.5_DP - param%h2 / param%Le_f
        
    enddo
  
    A(1,n)%l = -2.0_DP * param%h2 / param%Le_f
    A(1,n)%c = 2.0_DP * param%h2 / param%Le_f + rho(1,n)*rho(1,n) * param%beta*param%beta * &
    dexp((param%beta*(T(1,n)-1.0_DP))/(1.0_DP+param%gamma*(T(1,n)-1.0_DP)))*Z(1,n)
    
    A(2,n)%c = 1.0_DP
    
    b(:,:) = 0.0_DP
    b(1,1) = 1.0_DP
    b(2,n) = 1.0_DP
  end subroutine set_F
  
  subroutine act_F(A, b, rho, Z, T, param)
    type(sparse_t), intent(inout), dimension(:,:) :: A
    real(DP)      , intent(inout), dimension(:,:) :: b
    real(DP)      , intent(in   ), dimension(:,:) :: rho, Z, T
    type(param_t) , intent(in   ) :: param
    INTEGER :: i,n
    
    n = size(b, dim = 2)
  
    A(2,1)%c =  2.0_DP * param%h2 / param%Le_f + rho(2,1)** 2 * param%beta** 2 * &
    dexp((param%beta*(T(2,1)-1.0_DP))/(1.0_DP+param%gamma*(T(2,1) - 1.0_DP))) * Z(2,1)
  
    do i = 2,n-1
        A(:,i)%c = 2.0_DP * param%h2 / param%Le_f + rho(:,i)**2 * param%beta*param%beta * &
        dexp((param%beta*(T(:,i)-1.0_DP))/(1.0_DP+param%gamma*(T(:,i)-1.0_DP)))*Z(:,i)
    enddo
  
    A(1,n)%c = 2.0_DP * param%h2 / param%Le_f + rho(1,n)*rho(1,n) * param%beta*param%beta * &
    dexp((param%beta*(T(1,n)-1.0_DP))/(1.0_DP+param%gamma*(T(1,n)-1.0_DP)))*Z(1,n)
    
  end subroutine act_F

  subroutine set_Z(A, b, rho, F, T, param)
    type(sparse_t), intent(inout), dimension(:,:) :: A
    real(DP)      , intent(inout), dimension(:,:) :: b
    real(DP)      , intent(in   ), dimension(:,:) :: rho, F, T
    type(param_t) , intent(in   ) :: param
    INTEGER :: i,n
    
    n = size(F, dim = 2)
  
    A(:,:)%l = 0.0
    A(:,:)%c = 0.0
    A(:,:)%r = 0.0

    A(1,1)%c = 1.0_DP
     
    A(2,1)%c = 2.0_DP * param%h2 / param%Le_z + rho(2,1)** 2 - rho(2,1)** 2 * param%beta** 2 * &
    dexp((param%beta*(T(2,1)-1.0_DP))/(1.0_DP+param%gamma*(T(2,1)-1.0_DP))) * F(2,1)
    A(2,1)%r = -2.0_DP * param%h2 / param%Le_z
    
    do i = 2,n-1
        A(:,i)%l = -param%Pe(:) * param%h1 * 0.5_DP - param%h2 / param%Le_z
        A(:,i)%c = 2.0_DP * param%h2 / param%Le_z + rho(:,i)** 2 - rho(:,i)** 2 * param%beta** 2 * &
        dexp((param%beta*(T(:,i)-1.0_DP))/(1.0_DP+param%gamma*(T(:,i)-1d0))) * F(:,i)
        A(:,i)%r =  param%Pe(:) * param%h1 * 0.5_DP - param%h2 / param%Le_z
  
    enddo
    
    A(1,n)%l = -2.0_DP * param%h2 / param%Le_z
    A(1,n)%c = 2.0_DP * param%h2 / param%Le_z + rho(1,n)** 2  - rho(1,n)** 2 * param%beta** 2 * &
    dexp((param%beta*(T(1,n)-1.0_DP))/(1.0_DP+param%gamma*(T(1,n)-1.0_DP)))*F(1,n)
    
    A(2,n)%c = 1.0_DP
    
    b(:,:) = 0.0_DP
  end subroutine set_Z
  
  subroutine act_Z(A, b, rho, F, T, param)
    type(sparse_t), intent(inout), dimension(:,:) :: A
    real(DP)      , intent(inout), dimension(:,:) :: b
    real(DP)      , intent(in   ), dimension(:,:) :: rho, F, T
    type(param_t) , intent(in   ) :: param
    INTEGER :: i,n
    
    n = size(F, dim = 2)
  
    A(2,1)%c = 2.0_DP * param%h2 / param%Le_z + rho(2,1)** 2 - rho(2,1)** 2 * param%beta** 2 * &
    dexp((param%beta*(T(2,1)-1.0_DP))/(1.0_DP+param%gamma*(T(2,1)-1.0_DP))) * F(2,1)

    do i = 2,n-1
        A(:,i)%c = 2.0_DP * param%h2 / param%Le_z + rho(:,i)** 2 - rho(:,i)** 2 * param%beta** 2 * &
        dexp((param%beta*(T(:,i)-1.0_DP))/(1.0_DP+param%gamma*(T(:,i)-1d0))) * F(:,i)
    enddo
    
    A(1,n)%c = 2.0_DP * param%h2 / param%Le_z + rho(1,n)** 2  - rho(1,n)** 2 * param%beta** 2 * &
    dexp((param%beta*(T(1,n)-1.0_DP))/(1.0_DP+param%gamma*(T(1,n)-1.0_DP)))*F(1,n)
    
  end subroutine act_Z

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
  
  subroutine fixingt(T, Z, rho, i, param)
    real(DP)     , intent(inout) :: T(:,:)
    real(DP)     , intent(in   ) :: Z(:,:), rho(:,:)
    integer      , intent(in   ) :: i
    type(param_t), intent(inout) :: param
    double precision :: t_x,t_xx, reaction_heat, heat_ex
    integer :: ift
    
    ift = param%i_fixt
    T(1, ift) = 1d0
  
    t_x = (T(1,ift+1) - T(1,ift-1)) * 0.5_DP * param%h1
    t_xx = (T(1,ift+1) - 2.0_DP * T(1, ift) + T(1,ift-1)) * param%h2
    heat_ex = param%bK * (T(2,ift)-T(1, ift))
    reaction_heat = rho(1,ift)*rho(1,ift)* Z(1,ift)
    
    if (param%q_or_Pe == 1) then
      param%Pe(1) = param%Pe(1)*(1.0_DP-param%relaxPe)+&
      (param%relaxPe) * 1.0_DP / t_x * &
      (t_xx + param%q * reaction_heat + heat_ex )
  
      param%Pe(2) = -param%Pe(1)
    else
      param%q = param%q*(1.0_DP - param%relaxPe*i/param%m_iter) + &
      (param%relaxPe * i/param%m_iter - 0.0) * 1.0_DP /           &
      reaction_heat * (param%Pe(1) * t_x - t_xx - heat_ex)
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