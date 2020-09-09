module initial

  use variables_m 
  implicit none

contains

subroutine readparam()
    open(101,file='datainput.dat',status='old')
    read(101,*) ctedens              !is it cte the density?
    read(101,*) fixt
    read(101,*) q_or_Pe
    read(101,*) Le_f,Le_z        !lewis numbers
    read(101,*) beta,gamma     !zeldovich and gamma numbers
    read(101,*) q              !dimensionless heat release
    read(101,*) Pe(1), Pe(2)              !flow rate
    read(101,*) bK              !b
    read(101,*) relax, relaxPe  !relaxation factor 
    read(101,*) tol            !tolerance of the finite differences
    read(101,*) xmin,xmax      !domain
    read(101,*) maxiter        !maximum number of iterations
    read(101,*) show           !show the iteration
    read(101,*) x1change, x2change  !percentage of the distance that has heat exchange
    read(101,*) initial_data   !does it have initial data
    read(101,*) n
    close(101) 

end subroutine readparam

subroutine initialT(x, q, index_fixt) 
    real(kind=8), dimension(:), intent(inout) :: x
    DOUBLE PRECISION, INTENT(IN):: q
    integer, intent(in) :: index_fixt
    integer :: i,j
    x = q
    
    do i = 1, floor(index_fixt-n*0.1)
        x(i) = 0d0
    enddo
    j=floor(n*0.75)
    do i = floor(index_fixt-n*0.1), n
        x(i) = x(i-1)+1d0/(index_fixt-floor(index_fixt-n*0.1)+1)
        if (x(i) .ge. q*1.2) then
            x(i) = q *1.2
            j = i
            exit
        endif
    enddo
    do i = j, floor(j+n*0.2)
        x(i) = q*1.2
    enddo
    do i = floor(j+n*0.2), min(floor(j+n*0.3), floor(n*0.95))
        x(i) = x(i-1)-0.2*q/(floor(index_fixt+n*0.3)-floor(index_fixt+n*0.2)+1)
    enddo
endsubroutine

subroutine initialY(x,index_fixt) 
    real(kind=8), dimension(:), intent(inout) :: x
    integer, intent(in) :: index_fixt
    integer :: i
    x = 0d0
    do i = 1, floor(index_fixt-n*0.05)
        x(i) = 1d0

    enddo
    do i = floor(index_fixt-n*0.05), floor(index_fixt+n*0.05)
        x(i) = x(i-1)-1d0/(floor(n*0.1)+1)
    enddo
endsubroutine

subroutine initialZ(x, index_fixt) 
    real(kind=8), dimension(:), intent(inout) :: x
    integer, intent(in) :: index_fixt
    integer :: i
    x = 0d0
    
    do i = floor(index_fixt-n*0.05), index_fixt
        x(i) = x(i-1)+5d-1/(floor(n*0.05)+1)
    enddo
    do i = index_fixt, floor(index_fixt+n*0.05)
        x(i) = x(i-1)-5d-1/(floor(n*0.05)+1)
    enddo
endsubroutine

subroutine read_data(read_from)
    INTEGER, INTENT(IN) :: read_from
    INTEGER :: i
    DOUBLE PRECISION, DIMENSION(n) :: xcoord
    double precision               :: dummy, dummy2, dummy3, &
    dummy4, dummy5, dummy6
    
    OPEN(101, file = 'initial_profiles.txt')
    
    if(read_from == 1) then

    do i = 1,n
        read(101, *) xcoord(i), x(1,i), x(1,n+i), x(1,2*n+i), x(2,i),&
         x(2,n+i), x(2,2*n+i) 
    enddo
    else
    read(101, *)
    read(101, *)
    do i = 1,n
        read(101, *) xcoord(i), x(1,i), x(1,n+i), x(1,2*n+i), dummy, &
        dummy2, dummy3, x(2,i), x(2,n+i), x(2,2*n+i), dummy4, dummy5,&
        dummy6
    enddo
    endif
    
    CLOSE(101)
end subroutine read_data


end module
