module velocity

!===============================================================================!
!==============================----START----====================================!
!===============================================================================!

!<<----------------------------------- Module imports
    use    Soln
    use    functions
!<<----------------------------------- End imports

    implicit none

    contains

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME LOOP                                                              --!
!--    Performs the position and velocity iterations to try and improve a     --!
!--    solution.                                                              --!
!===============================================================================!

        subroutine LOOP(SolArray,    &
                        Darray,      &
                        Tarray,      &
                        solBest,     &
                        solMin,      &
                        solMax,      &
                        NoLoops,     &
                        fileOut,     &
                        ScalingOMax, &
                        ScalingOMin, &
                        LoopbofB)
            Type (sol), dimension(:), intent(inout) :: SolArray, SolBest
            Type (sol), dimension(:), intent(in)    :: solMin, solMax
            Type (sol), intent(in)                  :: LoopbofB
            real*8, dimension(:), intent(in)        :: Darray, Tarray
            integer, intent(in)                     :: NoLoops
            integer                                 :: i
            real*8                                  :: dummy
            character*80                            :: fileOut
            logical                                 :: xyzOut
            real*8, intent(in)                      :: ScalingOMax, ScalingOMin


            !open(unit=20,file=trim(fileOut), status='unknown')

            do i=1, NoLoops

                call UPDATEF(SolArray, Darray, Tarray, solBest) !Updates fitting value (functions module)

                call UPDATEV(SolArray, solBest, i, NoLoops, solMin, solMax, &
                    ScalingOMax, ScalingOMin, LoopbofB) !Updates velocity

                call UPDATEP(SolArray, solMin, solMax) !Updates position

                !call PRINTV(SolArray, SolBest, fileOut, i) !Print out all data

                !Please only use these if the problem is a 2 dimensional parameter space.
                !call XYZPRINT(SolArray, i) !Print out xyz data for plotting
                !call SXYZPRINT(solBest, i) !Prints swarm centres

                !if (i == NoLoops) call SOPTPRINT(solBest)


            end do

            !close(20)

        end subroutine

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME UPDATEV                                                           --!
!--    Updates the parameter space velocity.                                  --!
!===============================================================================!

        subroutine UPDATEV(cSol,        &
                           solBest,     &
                           it,          &
                           imax,        &
                           solMin,      &
                           solMax,      &
                           ScalingOMax, &
                           ScalingOMin, &
                           LoopbofB)
            Type (sol), intent(inout), dimension(:) :: cSol
            Type (sol), intent(in), dimension(:)    :: SolBest, solMin, solMax
            Type (sol), intent(in)                  :: LoopbofB
            Type (sol)                              :: SSolBest
            real*8, intent(in)                      :: ScalingOMax, ScalingOMin
            real*8                                  :: newV, phi1, phi2, phi3, input, lvmax, philocal, curVmax
            integer                                 :: time_array(8), i, inint, j, sNo
            integer, intent(in)                     :: it, imax

            phi1 = 2.0; phi2 = 2.0 ; phi3  = 1.0; lvmax = 0; philocal = 0

            !call date_and_time(values=time_array)
            !inint = (1000 * time_array(7)) + time_array(8)
            !input = rand(inint)

            do j=1, size(cSol)

                sNo      = cSol(j)%SwarmNo
                SSolBest = SolBest(sNo)

                do i=1, size(cSol(j)%vars)

                    lvmax = scalefunc( it, imax, solMax(sNo)%vmin(i), solMax(sNo)%vmax(i) )

                    ! you can play around with this, at the moment it has an attraction to
                    ! the best solution the swarm has found and the best solution found by all the
                    ! swarms
                    !
                    ! you can add more complex potentials or make the swarms repluse each other
                    ! or anything else you can play aorund with
                    newV = cSol(j)%velo(i) * scalefunc( it, imax, ScalingOMax, ScalingOMin ) + &
                           phi1 * rand()  * ( cSol(j)%bvars(i) -  cSol(j)%vars(i)  )     + &
                           phi2 * rand()  * ( SSolBest%bvars(i) - cSol(j)%vars(i)  )  !+ &
                        !phi3 * rand()  * ( LoopbofB%vars(i) -  cSol(j)%vars(i)  )  * heavy(it, 2*imax/3)! + &
                        !phi3 * rand() * repulsive(cSol(j)%SwarmNo, SolBest, i, cSol(j)%vars(i))

                    if (newV .gt. lvmax ) then
                        newV = lvmax
                    else if ( newV .lt. -lvmax ) then
                        newV = -lvmax
                    end if

                    cSol(j)%velo(i) = newV
                end do
            end do



        end subroutine

!===============================================================================!
!========================----Function----=======================================!
!--    NAME heavy                                                             --!
!--    Function for switching on and off interactions.                        --!
!===============================================================================!

        function heavy(itc, ito) result(gNumber)
            integer, intent(in) :: itc, ito
            real*8              :: gNumber

            gNumber = 0

            if (itc >= ito) gNumber = 1

        end function

!===============================================================================!
!========================----Function----=======================================!
!--    NAME scalefunc                                                         --!
!--    Scales a number between Nmin and imax in the range of 1 to imax.       --!
!===============================================================================!

        function scalefunc(it, imax, Nmin, Nmax) result(ScNumber)
            integer, intent(in) :: it, imax
            real*8, intent(in)  :: Nmin, Nmax
            real*8              :: ScNumber

            ScNumber = -( (real(it) - 1 ) * ( (Nmax - Nmin) / real(imax) ) ) + Nmax

        end function

!===============================================================================!
!========================----Function----=======================================!
!--    NAME repulsive                                                         --!
!--    A coulomb-like repulsive interaction between the particle and other    --!
!--    swarm's optimal position.                                              --!
!===============================================================================!

        function repulsive(sNum, BSols, pNum, variable) result(rNum)
            Type (sol), intent(in), dimension(:) :: BSols
            integer, intent(in)                  :: sNum, pNum
            real *8, intent(in)                  :: variable
            integer                              :: i
            real*8                               :: rNum

            rNum = 0

            do i=1, (size(BSols) - 1)
                if (i == sNum) cycle
                if ((variable - BSols(i)%bvars(pNum)) .ne. 0 ) then
                    rNum = rNum +  ( ( variable - BSols(i)%bvars(pNum) ) &
                            / ( abs(( variable - BSols(i)%bvars(pNum) )**2 )) )
                end if
            end do

        end function

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME UPDATEP                                                           --!
!--    Updates the parameter space positions                                  --!
!===============================================================================!

        subroutine UPDATEP(cSol, solMin, solMax)
            Type (sol), intent(inout),dimension(:) :: cSol
            Type (sol), intent(in),dimension(:)    :: solMin, solMax
            real*8                                 :: vmin  , vmax
            real*8                                 :: dummy
            integer                                :: i, j

            do j=1, size(cSol)
                do i=1, size(cSol(j)%vars)
                    dummy = cSol(j)%vars(i) + cSol(j)%velo(i)

                    !vmin = solMin%vars(i)
                    !vmax = solMax%vars(i)

!                    if (dummy .gt. vmax) then !Box size
!
!                        dummy = vmax
!                        cSol(j)%velo(i) = 0
!
!                    else if (dummy .lt. vmin) then
!
!                        dummy = vmin
!                        cSol(j)%velo(i) = 0
!
!                    end if

                !    print *, "j=", j, "i=", i
                !    print *, "old value", cSol(j)%vars(i)
                !    print *, "new value", dummy

                    cSol(j)%vars(i) = dummy

                end do
            end do



        end subroutine

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME PRINTV                                                            --!
!--    prints out all the data stored by the program in a formatted manor.    --!
!===============================================================================!

        subroutine PRINTV(cSol, SolBest, fileOut, loop)
            Type (sol), intent(in), dimension(:) :: cSol
            character*80                         :: fileOut
            Type (sol), intent(in), dimension(:) :: SolBest
            Type (sol)                           :: SSolBest
            integer                              :: i, imax
            integer                              :: loop

            imax = size(cSol)
            SSolBest = SolBest(size(SolBest))


            write(20,46)
            write(20,47) loop
            write(20,40)
            write(20,41)

            do i=1, size(cSol)

                write(20,42) i,  cSol(i)%vars( 1), cSol(i)%vars( 2), cSol(i)%fValue, &
                         cSol(i)%velo( 1), cSol(i)%velo( 2),                         &
                         cSol(i)%bvars(1), cSol(i)%bvars(2), cSol(i)%bValue,         &
                         cSol(i)%SwarmNo

            end do

            write(20,40)
            write(20,43) "Best  Solution"
            write(20,44)
            write(20,45) SSolBest%vars(1), SSolBest%vars(2), SSolBest%fValue

            write(20,40)
            write(20,46)



            40 format(85("-"))
            41 format(3("-"), "Soln.", 7x,  "a", 7x, "b",   &
                              9x,  "f", 6x, "av", 5x, "bv", &
                              7x, "ab", 6x, "bb", 8x, "fb", 4x, "SN", 3("-"))
            42 format(3("-"),  i5, 2x, f6.2, 2x, f6.2, 2x,  &
                               ES8.2E2, 3x, f5.2, 2x, f5.2, &
                               3x, f6.2, 2x, f6.2, 2x, ES8.2E2, 1x, i5, 3("-"))
            43 format(33("-"),3x,a14,2x,33("-") )
            44 format(23("-"),12x,"a",6x,"b",9x,"f",9x,23("-") )
            45 format(23("-"),7x,f6.2,1x,f6.2,2x,ES8.2E2,9x,23("-") )
            46 format(85("*"))
            47 format(33("*"),3x,"ITERATION", i5,3x,32("*"))


        end subroutine

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME    XYZPRINT                                                       --!
!--    Prints the parameter space position and the current fitted value for   --!
!--    the iteration. Please only use for 2D parameter spaces.                --!
!===============================================================================!

        subroutine XYZPRINT(cSol, iteration)
            Type (sol), intent(in), dimension(:) :: cSol
            integer, intent(in)                  :: iteration
            character*80                         :: fileName, dummy
            integer                              :: i

            write(dummy,*) iteration

            fileName = "/home/jake/PSOFortran/xyz/xyz" // trim(adjustl(dummy)) // ".dat"

            open(unit=21,file=trim(fileName), status='unknown')

            do i=1, size(cSol)

                write(21,48) cSol(i)%vars( 1), cSol(i)%vars( 2), cSol(i)%fValue

            end do

            48 format(f8.2, 8x, f8.2, 8x, f8.2)

            close(21)

        end subroutine

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME    SXYZPRINT                                                      --!
!--    Prints the optimal parameter space position of each swarm and the      --!
!--    current fitted value for the iteration. Please only use for 2D         --!
!--    parameter spaces.                                                      --!
!===============================================================================!

        subroutine SXYZPRINT(SSol, iteration)
            Type (sol), intent(in), dimension(:) :: SSol
            integer, intent(in)                  :: iteration
            character*80                         :: fileName, dummy
            integer                              :: i

            write(dummy,*) iteration

            fileName = "/home/jake/PSOFortran/xyz/swarm" // trim(adjustl(dummy)) // ".dat"

            open(unit=21,file=trim(fileName), status='unknown')

            do i=1, (size(SSol) - 1)

                write(21,48) SSol(i)%vars( 1), SSol(i)%vars( 2), SSol(i)%fValue

            end do

            48 format(f8.2, 8x, f8.2, 8x, f8.2)

            close(21)

        end subroutine

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME    SOPTPRINT                                                      --!
!--    Prints the optimal parameter space position of each swarm and the      --!
!--    current fitted value for the iteration. Please only use for 2D         --!
!--    parameter spaces.                                                      --!
!===============================================================================!

        subroutine SOPTPRINT(SSol)
            Type (sol), intent(in), dimension(:) :: SSol
            character*80                         :: fileName
            integer                              :: i


            fileName = "/home/jake/PSOFortran/PSOData2.dat"

            open(unit=21,file=trim(fileName), status='unknown')

            do i=1, (size(SSol) - 1)

                write(21,49) SSol(i)%vars( 1), SSol(i)%vars( 2), SSol(i)%fValue

            end do

            49 format(f16.8, 8x, f16.8, 8x, f16.8)

            close(21)

        end subroutine

!===============================================================================!
!==============================-----END-----====================================!
!===============================================================================!
end module
