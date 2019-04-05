module Soln
!===============================================================================!
!==============================----START----====================================!
!===============================================================================!

    implicit none

!===============================================================================!
!==============================----TYPE----=====================================!
!--                NAME    sol                                                --!
!--    Contains parameter array (vars), velocity array (velo), current        --!
!--    fitting function value (fValue) and the best values for the parameter  --!
!--    and fitting value (bvars and bValue).                                  --!
!===============================================================================!
    type  sol
        private
        real*8, public                            :: fValue, bValue
        integer, public                           :: SwarmNo
        real*8, public, allocatable, dimension(:) :: vars, bvars, velo, vmax, vmin
    end type

    contains

        subroutine INITIATE( cSol, vNum )
            type(sol), intent(inout) :: cSol
            integer, intent(in)      :: vNum

            allocate(cSol%vars(  vNum  )); allocate(cSol%bvars( vNum  )); allocate(cSol%velo(  vNum  ))
            allocate(cSol%vmax(  vNum  )); allocate(cSol%vmin(  vNum  ));
        end subroutine

        subroutine INITIATEA( cSol, vNum )
            type(sol), dimension(:), intent(inout) :: cSol
            integer, intent(in)                    :: vNum
            integer                                :: i

            do i=1, size(cSol)
                call INITIATE(cSol(i),vNum)
            end do
        end subroutine

!===============================================================================!
!==============================-----END-----====================================!
!===============================================================================!

end module
