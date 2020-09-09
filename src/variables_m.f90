module variables_m
  
  implicit none
  
  integer , parameter :: DP = selected_real_kind(15 , 307)
  
  type, public :: sparse_t
  ! left, central, right diagonal
  real(DP) :: l, c, r

  ! contains
  ! procedure, pass(self) :: sparse_mult
  ! procedure, pass(self) :: sparse_div
  ! generic :: operator(*) => sparse_mult!, sparse_mult_reverse
  ! generic :: operator(/) => sparse_div!, sparse_div_reverse
  end type sparse_t
  
  ! interface operator (*)
  ! include sparse_t
  ! function sparse_mult(A, r) result (B)
  !   type(sparse_t), intent(in   ) :: A
  !   real(DP)      , intent(in   ) :: r
  !   type(sparse_t) :: B
  ! end function sparse_mult
  
  ! function sparse_mult_reverse(r, A) result (B)
  !   type(sparse_t), intent(in   ) :: A
  !   real(DP)      , intent(in   ) :: r
  !   type(sparse_t), intent(  out) :: B
  ! end function sparse_mult_reverse
  ! end interface
  
  ! interface operator (/)
  ! function sparse_mult(A, r) result (B)
  !   type(sparse_t), intent(in   ) :: A
  !   real(DP)      , intent(in   ) :: r
  !   type(sparse_t), intent(  out) :: B
  ! end function sparse_mult
  
  ! function sparse_mult_reverse(r, A) result (B)
  !   type(sparse_t), intent(in   ) :: A
  !   real(DP)      , intent(in   ) :: r
  !   type(sparse_t), intent(  out) :: B
  ! end function sparse_mult_reverse
  ! end interface
  type(sparse_t), ALLOCATABLE, dimension(:,:), SAVE :: A 
  real(DP), ALLOCATABLE, dimension(:,:),  SAVE :: b , x , xp, &
  rho, omega_c, omega_b
  real(DP), allocatable, dimension(:), save :: xgrid
  real(DP), SAVE :: Le_f, Le_z, beta, gamma, q, bK, relax, &
  tol, xmin, xmax, h, x1change, x2change, time_start, time_finish, fixt,&
  relaxPe
  real(DP), DIMENSION(2), SAVE ::  Pe 
  INTEGER, SAVE :: n, show, maxiter, nmin, nmax, initial_data, & 
  ctedens,index_fixt, q_or_Pe 
  CHARACTER(400), SAVE :: file_output = ''
  
  contains

!   pure function sparse_mult(self, r) result (B_spr)
!       type(sparse_t), intent(in   ) :: self
!       real(DP)      , intent(in   ) :: r
!       type(sparse_t) :: B_spr

!       B_spr%l = self%l * r
!       B_spr%c = self%c * r
!       B_spr%r = self%r * r
!   end function sparse_mult

! !   function sparse_mult_reverse(r, self) result (B_spr)
! !     type(sparse_t), intent(in   ) :: self
! !     real(DP)      , intent(in   ) :: r
! !     type(sparse_t) :: B_spr

! !     B_spr%l = self%l * r
! !     B_spr%c = self%c * r
! !     B_spr%r = self%r * r
! ! end function sparse_mult_reverse

! pure function sparse_div(self, r) result (B_spr)
!   type(sparse_t), intent(in   ) :: self
!   real(DP)      , intent(in   ) :: r
!   type(sparse_t) :: B_spr

!   B_spr%l = self%l / r
!   B_spr%c = self%c / r
!   B_spr%r = self%r / r
! end function sparse_div

! ! function sparse_div_reverse(r, self) result (B_spr)
! ! type(sparse_t), intent(in   ) :: self
! ! real(DP)      , intent(in   ) :: r
! ! type(sparse_t) :: B_spr

! ! B_spr%l = self%l / r
! ! B_spr%c = self%c / r
! ! B_spr%r = self%r / r
! ! end function sparse_div_reverse
end module variables_m