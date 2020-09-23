module param_m

  use type_m, only: DP
  
  implicit none

  type, public :: param_t 

    ! Parameters of the mixture
    real(DP) :: Pe(2), Le_f, Le_z, beta, gamma, q, bK

    ! Parameters of the domain
    real(DP) :: xmin, xmax, x1c, x2c, fixt, h, h1, h2
    
    real(DP) :: relax, relaxPe, tol

    integer :: n, nmin, nmax, m_iter, show
    
    integer :: initial_data, cte_rho, i_fixt, q_or_Pe

    contains 

    procedure :: read => read_param 
  end type param_t

  contains

  subroutine read_param(this)

    class(param_t), intent(inout) :: this
    real(DP), allocatable, dimension(:) :: xgrid
    integer :: i

    open(101,file='datainput.dat',status='old')
    read(101,*) this%cte_rho              !is it cte the density?
    read(101,*) this%fixt
    read(101,*) this%q_or_Pe
    read(101,*) this%Le_f, this%Le_z        !lewis numbers
    read(101,*) this%beta, this%gamma     !zeldovich and gamma numbers
    read(101,*) this%q              !dimensionless heat release
    read(101,*) this%Pe(1), this%Pe(2)              !flow rate
    read(101,*) this%bK              !b
    read(101,*) this%relax, this%relaxPe  !relaxation factor 
    read(101,*) this%tol            !tolerance of the finite differences
    read(101,*) this%xmin, this%xmax      !domain
    read(101,*) this%m_iter        !maximum number of iterations
    read(101,*) this%show           !show the iteration
    read(101,*) this%x1c, this%x2c  !percentage of the distance that has heat exchange
    read(101,*) this%initial_data   !does it have initial data
    read(101,*) this%n
    close(101) 
    
    ! read input constants file
    
    this%gamma = this%q / (this%q + 1.0_DP)
    
    this%h = DBLE(this%xmax-this%xmin)/DBLE(this%n-1)
    
    do i = 1, this%n
      if( abs(this%xmin + (i-1) * this%h - this%x1c) <= 10d-6 ) then
        this%nmin = i
        exit
      endif
      if (i == this%n) write(*,*) "ERROR: COGE EL NUMERO DE PUNTOS MULTIPLO " &
      // "DE LA LONGITUD DEL DOMINIO"
    enddo
    do i =1, this%n
      if( abs(this%xmin + (i-1) * this%h - this%x2c) <= 10d-6 ) then
        this%nmax = i
        exit
      endif
    enddo
    
    this%i_fixt = floor(0.4 * this%n)
    if(this%fixt <= this%xmax) then 
      allocate(xgrid(this%n))
      xgrid = (/(this%xmin + this%h*(i-1), i=1, this%n)/)
      this%i_fixt = minloc(abs(xgrid - this%fixt), dim =1)
      deallocate(xgrid)
    endif
    
    this%h1 = 1.0_DP / (this%h)
    this%h2 = 1.0_DP / (this%h ** 2)
    
  end subroutine read_param
end module param_m