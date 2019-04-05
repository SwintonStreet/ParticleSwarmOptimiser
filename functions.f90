module functions

!===============================================================================!
!==============================----START----====================================!
!===============================================================================!

!<<----------------------------------- Module imports
    use    Soln
!<<----------------------------------- End imports

    implicit none

    contains

!===============================================================================!
!========================----Function----=======================================!
!--    NAME func                                                              --!
!--    Calculates the fitting function for a given solution.                  --!
!===============================================================================!

        function fitfunc(CurSol, Darray, Tarray) result(ffNum)    !The fitting function
            Type (sol), intent(in)           :: CurSol
            real*8, dimension(:), intent(in) :: Darray, Tarray
            real*8                           :: ffNum
            integer                          :: i

            ffNum=0

            do i=1, size(Darray)
                ffNum = ffNum + ( testfunc(CurSol,Darray(i)) - Tarray(i) )**2
            end do

        end function

!===============================================================================!
!========================----Function----=======================================!
!--    NAME    testfunc                                                       --!
!--    Calculates the value of the test function for the current solution.    --!
!--                                                                           --!
!--    Implement the function being tested here.                              --!
!===============================================================================!

        function testfunc(CurSol, Dependent) result(tfNum)
            Type (sol), intent(in) :: CurSol
            real*8, intent(in)     :: Dependent
            real*8                 :: tfNum

            ! dummy model being calculated
            tfNum =   (CurSol%vars(1) / 10 ) * Cos(2*Dependent + 5*CurSol%vars(1)) &
                    + (CurSol%vars(2) / 10 ) * Sin(3*Dependent + 5*CurSol%vars(2)) &
                    - (CurSol%vars(3) / 5  ) * Cos(5*Dependent + 5*CurSol%vars(3))

        end function

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME    UPDATEF                                                        --!
!--    Subroutine which updates a solutions fitting function value.           --!
!===============================================================================!

        subroutine UPDATEF(Sols, Darray, Tarray, solBest)
            Type (sol), dimension(:), intent(inout) :: Sols
            Type (sol), dimension(:), intent(inout) :: SolBest
            Type (sol)                              :: SSolBest
            real*8, dimension(:), intent(in)        :: Darray, Tarray
            integer                                 :: i
            real*8                                  :: dummy

            do i=1, size(Sols)
                dummy          = fitfunc(Sols(i), Darray, Tarray)
                Sols(i)%fValue = dummy
                SSolBest       = SolBest(Sols(i)%SwarmNo)

                !---- Update the swarm best value
                if ( dummy .lt. SSolBest%fValue .or. SSolBest%fValue .eq. -1) then
                    SolBest(Sols(i)%SwarmNo) = Sols(i)
                end if

                !---- Update the global best value
                if ( dummy .lt. SolBest(size(SolBest))%fValue .or. SSolBest%fValue .eq. -1) then
                    SolBest(size(SolBest)) = Sols(i)
                end if

                !---- Update the particle best value
                if ( Sols(i)%fValue .lt. Sols(i)%bValue .or. Sols(i)%bValue .eq. -1 ) then

                    Sols(i)%bValue = dummy
                    Sols(i)%bvars  = Sols(i)%vars
                end if

            end do

        end subroutine

!===============================================================================!
!==============================-----END-----====================================!
!===============================================================================!
end module
