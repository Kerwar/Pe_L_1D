program main
  
  use type_m, only: DP
  use param_m, only: param_t
  use functions
  use initial
  use variables_m  
  
  implicit none
  
  type(sparse_t), allocatable, dimension(:,:) :: A_T, A_F, A_Z 
  type(param_t ) :: param 
  
  real(DP), ALLOCATABLE, dimension(:,:),  SAVE :: bT, bF, bZ, T, F, Z, Tp, Fp, &
  Zp, rho, omega_c, omega_b

  real(DP) :: t_start, t_finish, Pe_mean(1000)
  
  CHARACTER(400) :: file_output = ''
  
  integer :: i,j 
  
  call cpu_time(t_start)
  
  call param%read()
  
  allocate(A_T(2,param%n), A_F(2,param%n), A_Z(2,param%n), &
  bT(2,param%n), bF(2,param%n), bZ(2,param%n), T(2,param%n), &
  F(2,param%n), Z(2,param%n), Tp(2,param%n), Fp(2,param%n),  &
  Zp(2,param%n), rho(2,param%n), omega_b(2,param%n), omega_c(2,param%n))
  
  
  write(*,'(A8, 7A25)') 'Iter', 'Error', 'Pe', 'q'               , & 
  'Maximum Temperature 1', 'Maximum Temperature 2', 'Maximum Z 1', &
  'Maximum Z 2'
  print '(A190)', repeat('-',190)
  
  
  if (param%initial_data == 0) then
    CALL initialT(T(1, :), param)
    CALL initialT(T(2, param%n:1:-1), param)
    
    CALL initialY(F(1, :), param)
    CALL initialY(F(2, param%n:1:-1), param)
    
    CALL initialZ(Z(1, :), param)
    CALL initialZ(Z(2, param%n:1:-1), param)
  else
    CALL read_data(T, F, Z, param)
  endif
  
  OPEN(301,file = 'resultstemp.txt')
  do j = 1, param%n
    WRITE(301,'(9F20.8)') param%xmin + (j-1) * param%h, T(1,j), F(1,j), Z(1,j), &
    T(2,j), F(2,j), Z(2,j)   
  enddo
  CLOSE(301)
  
  rho = 1.0_DP / (1.0_DP + (param%gamma * T(:,:)/(1.0_DP- param%gamma)))
  if (param%cte_rho == 1) rho = 1.0_DP
  
  Tp = T
  Fp = F
  Zp = Z

  Pe_mean = param%Pe(1)
  
  call set_t(A_T(:,:), bT(:,:), rho, Z(:,:), T(:,:), param)
  call set_F(A_F(:,:), bF(:,:), rho, Z(:,:), T(:,:), param)
  call set_Z(A_Z(:,:), bZ(:,:), rho, F(:,:), T(:,:), param)

  do i = 1, param%m_iter  

    if (abs(param%fixt) < param%xmax ) &
      call set_t(A_T(:,:), bT(:,:), rho, Z(:,:), T(:,:), param)
    call act_t(A_T(:,:), bT(:,:), rho, Z(:,:), T(:,:),param)
    CALL onestepgaussseidel(A_T(1,:), bT(1,:), T(1,:), Tp(1,:))
    CALL onestepgaussseidel(A_T(2,:), bT(2,:), T(2,:), Tp(2,:))

    T(:,:) = (1.0_DP-param%relax)*Tp(:,:) + param%relax*T(:,:)
    !T(:,:) = max(0.0_DP, T(:,:))

    if (param%cte_rho == 0) rho = 1.0_DP / &
    (1.0_DP + param%gamma * T(:,:)/(1.0_DP - param%gamma))
    
    if (abs(param%fixt) < param%xmax ) &
      call set_F(A_F(:,:), bF(:,:), rho, Z(:,:), T(:,:), param)
    call act_F(A_F(:,:), bF(:,:), rho, Z(:,:), T(:,:), param)
    CALL onestepgaussseidel(A_F(1,:), bF(1,:), F(1,:), Fp(1,:))
    CALL onestepgaussseidel(A_F(2,:), bF(2,:), F(2,:), Fp(2,:))
    F(:,:) = (1.0_DP-param%relax) * Fp(:,:) + param%relax * F(:,:)
    !F(:,:) = max(0.0_DP, F(:,:))

    if (abs(param%fixt) < param%xmax ) &
      call set_Z(A_Z(:,:), bZ(:,:), rho, F(:,:), T(:,:), param)
    call act_Z(A_Z(:,:), bZ(:,:), rho, F(:,:), T(:,:), param)
    CALL onestepgaussseidel(A_Z(1,:), bZ(1,:), Z(1,:), Zp(1,:))
    CALL onestepgaussseidel(A_Z(2,:), bZ(2,:), Z(2,:), Zp(2,:))
    Z(:,:) = (1.0_DP-param%relax)*Zp(:,:) + param%relax * Z(:,:)
    !Z(:,:) = max(0.0_DP, Z(:,:))

    if(abs(param%fixt) <= param%xmax ) &
      call fixingt(T(:,:), Z(:,:), rho(:,:), i, param)

    if (mod(i, param%show) == 0) then
      !print *, maxval(T-Tp), maxval(F-Fp), maxval(Z-Zp)
      if (mod(i, param%show*10) == 0) then
        print '(A190)', repeat('-',190)
        write(*,'(A8, 7A25)') 'Iter', 'Error', 'Pe', 'q'               , & 
        'Maximum Temperature 1', 'Maximum Temperature 2', 'Maximum Z 1', &
        'Maximum Z 2'
        print '(A190)', repeat('-',190)
      endif
      print '(i8, 7F25.13)', i, max(maxval(T-Tp), maxval(F-Fp), maxval(Z-Zp)), &
       param%Pe(1), param%q, maxval(T(1,:)), maxval(T(2,:)), maxval(Z(1,:)),   &
       maxval(Z(2,:)) 

      OPEN(301,file = 'resultstemp.txt')
      do j = 1, param%n
        WRITE(301,'(9F20.8)') param%xmin + (j-1) * param%h, T(1,j), F(1,j), Z(1,j), &
        T(2,j), F(2,j), Z(2,j)   
      enddo
      CLOSE(301)
    end if
    if (max(maxval(T-Tp), maxval(F-Fp), maxval(Z-Zp)) <= param%tol) exit
    
    Tp = T
    Fp = F
    Zp = Z
  enddo
  
  omega_b = rho**2 * param%beta** 2 * dexp((param%beta*(T(:,:)-1.0_DP)) &
  /(1.0_DP+param%gamma*(T(:,:)-1.0_DP)))*Z(:,:)*F(:,:)
  omega_c = rho**2 * Z(:,:)
  
  write(file_output,17) param%n, param%xmax, param%x2c, param%cte_rho, &
  param%Le_f, param%Le_z, param%q, param%beta, param%gamma, param%bK, param%Pe(1)
  
  
  if (isnan(param%Pe(1)) .eqv. .false.) then
    OPEN(201,file = file_output,status='replace')
    OPEN(401,file = 'initial_profiles.txt')
    write(201, '(A3,9E22.12)') "#", param%xmax, param%x2c, param%Le_F, &
    param%Le_Z, param%q, param%beta, param%gamma, param%bK, param%Pe(1)
    write(201, '(13A22)') '"X"', '"T_1"', '"F_1"', '"Z_1"', '"rho1"', &
    '"omega_b1"', '"omega_c1"', '"T_2"', '"F_2"', '"Z_2"', '"rho2"', &
    '"omega_b2"', '"omega_c2"'
    do i = 1, param%n
      WRITE(201,'(13E22.12)') param%xmin + (i-1) * param%h, T(1,i), F(1,i), Z(1,i), &
      1.0_DP/rho(1,i), omega_b(1,i), omega_c(1,i), T(2,i), F(2,i),         &
      Z(2,i), 1.0_DP/rho(2,i), omega_b(2,i), omega_c(2,i)
      WRITE(401,'(9E22.12)') param%xmin + (i-1) * param%h, T(1,i), F(1,i), Z(1,i), &
      T(2,i), F(2,i), Z(2,i)    
    enddo
    CLOSE(201) 
    CLOSE(401)
  endif
  
  
  deallocate(A_T, A_F, A_Z, bT, bF, bZ, T, F, Z, Tp, Fp,  &
  Zp, rho, omega_b, omega_c)
  
  CALL CPU_TIME(t_finish)
  
  print *, int(t_finish-t_start)/3600,':' ,mod(int(t_finish-t_start)/60,60), &
  ':',mod(int(t_finish-t_start), 60) 
  
  17 format('./Results/results_n-',I0,'_range','-',F0.0,'-',F0.0,'_cterho-',i1,&
  '_lef-',F0.1,'_lez-',F0.1,'_q-',F0.1, '_beta-',F0.1,'_gamma',F0.1,'_b-',F0.2,&
  '_Pe-',F0.4,'.txt')
end program main
