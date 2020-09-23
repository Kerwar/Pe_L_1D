module initial
  
  use type_m , only: DP
  use param_m, only: param_t
  implicit none
  
  contains
  
  
  
  subroutine initialT(x, param) 
    real(DP), dimension(:), intent(inout) :: x
    type(param_t)         , intent(in   ) :: param
    integer :: i,j
    
    x = param%q
    
    do i = 1, floor(param%i_fixt - param%n * 0.1)
      x(i) = 0d0
    enddo
    
    j = floor(param%n * 0.75)
    
    do i = floor(param%i_fixt - param%n * 0.1), param%n
      
      x(i) = x(i-1) + 1d0/(param%i_fixt -floor(param%i_fixt - param%n*0.1)+1)
      if (x(i) .ge. param%q*1.2) then
        x(i) = param%q * 1.2
        j = i
        exit
      endif
    enddo
    
    do i = j, floor(j + param%n*0.2)
      x(i) = param%q*1.2
    enddo

    do i = floor(j + param%n * 0.2), &
      min(floor(j + param%n * 0.3), floor(param%n * 0.95))
      
      x(i) = x(i-1) - 0.2 * param%q / (floor(param%i_fixt + param%n * 0.3) - &
      floor(param%i_fixt + param%n * 0.2) + 1)
    enddo
  endsubroutine
  
  subroutine initialY(x, param) 
    real(DP), dimension(:), intent(inout) :: x
    type(param_t)         , intent(in   ) :: param
    integer :: i
    
    x = 0.0_DP
    
    do i = 1, floor(param%i_fixt - param%n * 0.05_DP)
      x(i) = 1.0_DP
    enddo

    do i = floor(param%i_fixt - param%n * 0.05_DP), &
      floor(param%i_fixt + param%n * 0.05_DP)
      
      x(i) = x(i-1)-1.0_DP/(floor(param%n*0.1_DP)+1)
    enddo
  endsubroutine
  
  subroutine initialZ(x, param) 
    real(DP), dimension(:), intent(inout) :: x
    type(param_t)         , intent(in   ) :: param
    integer :: i
    
    x = 0.0_DP
    
    do i = floor(param%i_fixt - param%n*0.05), param%i_fixt
      x(i) = x(i-1) + 0.5_DP/(floor(param%n * 0.05)+1)
    enddo
    do i = param%i_fixt, floor(param%i_fixt + param%n * 0.05)
      x(i) = x(i-1) - 0.5_DP/(floor(param%n * 0.05)+1)
    enddo
  endsubroutine
  
  subroutine read_data(T, F, Z, param)
    real(DP), dimension(:,:), intent(inout) :: T, F, Z
    type(param_t)         , intent(in   ) :: param
    real(DP) :: xcoord
    real(DP) :: dummy, dummy2, dummy3, &
    dummy4, dummy5, dummy6
    integer  :: i
    
    OPEN(101, file = 'initial_profiles.txt')
    
    if(param%initial_data == 1) then
      
      do i = 1, param%n
        read(101, *) xcoord, T(1,i), F(1,i), Z(1,i), T(2,i), F(2,i), Z(2,i) 
      enddo
    else
      read(101, *)
      read(101, *)
      do i = 1, param%n
        read(101, *) xcoord, T(1,i), F(1,i), Z(1,i), dummy, dummy2, &
        dummy3, T(2,i), F(2,i), Z(2,i), dummy4, dummy5, dummy6
      enddo
    endif
    
    CLOSE(101)
  end subroutine read_data
  
  
end module
