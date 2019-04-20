program PSOTest
!===============================================================================!
!==============================----START----====================================!
!===============================================================================!

!===============================================================================!
!            Particle Swarm Optimiser                                           !
!===============================================================================!
!    Code performs particle swarm optimisation based on three (four) inputs:    !
!    1. Function with parameters to be optomised (test function).               !
!    2. Fitting function to test accuracy of chosen parameter.                  !
!    3. Data to test against.                                                   !
!    (4.) Number of solutions and number of iterations to try.                  !
!===============================================================================!

!<<----------------------------------- Module imports
    use    Soln      ! Contains the solution class.
    use    Startup   ! Initialises the solutions.
    use    functions ! Contains the fitting and test function.
    use    velocity  ! Performs the loop of solution velocity and position updates.
!<<----------------------------------- End imports

    implicit none

    integer                               :: Num, NumS, VarNo, ItNum, i, NoLoops, j, n
    type (sol), dimension(:), allocatable :: solArray, solBest, solBoftheB, solMin, solMax
    character*80                          :: fileNameD, fileNameT, fileOut
    real*8, dimension(:), allocatable     :: Darray, Tarray
    integer                               :: time_array(8), inint
    real*8                                :: ScalingOMax, ScalingOMin

    !Time.
    call date_and_time(values=time_array)
    inint = (60000 * time_array(6)) + (1000 * time_array(7)) + time_array(8)

    !General
    NumS    = 10  !Number of swarms
    Num     = 10  !Number of solutions to try per swarm
    VarNo   = 3   !Number of parameters to be optimised
    ItNum   = 150 !Number of iterations
    NoLoops = 100 !Number of attempts of the algorithm

    !Omega
    ScalingOMax = 0.9 !Max value for omega update scaling
    ScalingOMin = 0.4 !Min value for omega update scaling

    allocate(solArray( Num * NumS ))    !Allocates the number of solutions.
    allocate(solBest( NumS + 1 ))       !Allocates the array with the number of swarms + 1 for the best solution found.
    allocate(solBoftheB( NoLoops + 1 )) !Allocates the best of each attempt
    allocate(solMin( NumS ))            !Allocates the swarm min value for initiation
    allocate(solMax( NumS ))            !Allocates the swarm max value for initiation


    call INITIATEA( solMin     , VarNo ) !-----------------------------------!
    call INITIATEA( solMax     , VarNo ) !Sets up the arrays for the sol type.
    call INITIATEA( solArray   , VarNo ) !-----------------------------------!
    call INITIATEA( solBest    , VarNo ) !-----------------------------------!
    call INITIATEA( solBoftheB , VarNo ) !-----------------------------------!

    !Presets the best value so it is updated automatically in start up routine.
    solBest%fValue = -1

    !Presets the Best of the best value
    solBoftheB( NoLoops + 1 )%vars   = 0
    solBoftheB( NoLoops + 1 )%fValue = -1

    !Parameter minimum and maximum values for search space
    do i = 1, size( solMin )
        solMin(i)%vars( 1 ) = -200 ; solMin(i)%vars( 2 ) = -200  ; solMin(i)%vars( 3 ) = -200
        solMax(i)%vars( 1 ) =  200 ; solMax(i)%vars( 2 ) =  200  ; solMax(i)%vars( 3 ) = 200
        solMax(i)%vmax( 1 ) =  10  ; solMax(i)%vmax( 2 ) =  10   ; solMax(i)%vmax( 3 ) = 10
        solMax(i)%vmin( 1 ) =  1   ; solMax(i)%vmin( 2 ) =  1    ; solMax(i)%vmin( 3 ) = 1
    end do

    !Parameter minimum and maximum values for local search space
    solMin(1)%vars( 1 ) = -2     ; solMin(1)%vars( 2 ) = -2    ; solMin(1)%vars( 3 ) = -2
    solMax(1)%vars( 1 ) =  2     ; solMax(1)%vars( 2 ) =  2    ; solMax(1)%vars( 3 ) = 2
    solMax(1)%vmax( 1 ) =  0.1   ; solMax(1)%vmax( 2 ) =  0.1  ; solMax(1)%vmax( 3 ) = 0.1
    solMax(1)%vmin( 1 ) =  0.01  ; solMax(1)%vmin( 2 ) =  0.01 ; solMax(1)%vmin( 3 ) = 0.01

    !Files names
    fileNameD = "xdata.dat"        !Contains the dependent variable.
    fileNameT = "ydata.dat"        !Contains the cross reference data used to test the fits.
    fileOut   = "PSOData.dat"    !Output file name.

    !Counts the best of the best updates
    n = 0

    call STARTREAD(fileNameD, Darray, fileNameT, Tarray)

    do i=1, NoLoops

!===============================================================================!
!=======    START initialises the solutions (Startup module)    ========!
!===============================================================================!
        call STARTARRAY(SolArray, solMin, solMax, VarNo, Num)

!===============================================================================!
!=======    LOOP performs the solution iterations (velocity module)    ========!
!===============================================================================!
        call LOOP(SolArray,    &
                  Darray,      &
                  Tarray,      &
                  solBest,     &
                  solMin,      &
                  solMax,      &
                  ItNum,       &
                  fileOut,     &
                  ScalingOMax, &
                  ScalingOMin, &
                  solBoftheB( NoLoops + 1 ))
!===============================================================================!
!==========================---break---==========================================!
!===============================================================================!


        !Output commands.
        !print *, "Attempt ", i
        !print *, "fBest = ", solBest(size(solBest))%fValue
        !print *, "a = ", solBest(size(solBest))%vars(1), "b = ", solBest(size(solBest))%vars(2)

        solBoftheB(i) = solBest(size(solBest))

        if ( solBoftheB(i)%fValue .lt. solBoftheB( NoLoops + 1 )%fValue  .or. solBoftheB( NoLoops + 1 )%fValue .eq. -1  ) then
                solBoftheB( NoLoops + 1 ) = solBoftheB(i)
                if(i .gt. 1) write(*,fmt="(4(a,f9.4),a,i4)") &
                                " | a:", solBoftheB( NoLoops + 1 )%vars(1) &
                               , " b:", solBoftheB( NoLoops + 1 )%vars(2) &
                               , " c:", solBoftheB( NoLoops + 1 )%vars(3) &
                               , " f:", solBoftheB( NoLoops + 1 )%fValue, &
                               " Update ", n

                n = n + 1

                do j= 1, size(solMin(1)%vars)
                    !Update for parameter minimum and maximum values for local search space
                    solMin(1)%vars(j) = solBoftheB(i)%vars(j) - 0.01
                    solMax(1)%vars(j) = solBoftheB(i)%vars(j) + 0.01
                end do
        end if


        if (  NoLoops > 1 ) then
            deallocate(solBest); deallocate(solArray)

            allocate( solArray(Num * NumS)); allocate( solBest( NumS + 1  ) )

            call INITIATEA( solArray,  VarNo )    !-----------------------------------!
            call INITIATEA( solBest ,  VarNo )    !-----------------------------------!
            solBest%fValue = -1    !Presets the best value so it is updated automatically in start up routine.

            !allow one swarm to remember the best position from all the data
            solBest(1) = solBoftheB( NoLoops + 1 )
        end if

!        write(*,fmt="(a1,a,f6.2,a,a,f8.3,a,f8.3,a,f8.3)",advance="no") achar(13), "Percentage complete:"   , &
!                                            (real(i)/real(NoLoops)) * 100, &
!                            "%", " a:", solBoftheB( NoLoops + 1 )%vars(1) &
!                               , " b:", solBoftheB( NoLoops + 1 )%vars(2) &
!                               , " f:", solBoftheB( NoLoops + 1 )%fValue

        write(*,fmt="(a1,a,f6.2,a)",advance="no") achar(13), "Percentage complete:",(real(i)/real(NoLoops)) * 100,"%"


    end do


    !Time.
    call date_and_time(values=time_array); inint = (60000 * time_array(6)) + (1000 * time_array(7)) + time_array(8) - inint

    write(*,fmt="(a)",advance="yes") ""
    write(*,fmt="(a,i6)",advance="yes") "Final time: ", inint
    write(*,fmt="(a)",advance="yes") "Finished"

    if ( NoLoops > 1 ) then
        call SOPTPRINT( solBoftheB )    !Prints swarm centres
    else
        call SOPTPRINT( solBest )    !Prints swarm centres
    end if



    open(unit=10, file="stats")

    write(10,50)
    write(10,51)
    write(10,50)
    write(10,40)    NumS
    write(10,50)
    write(10,41)    Num
    write(10,50)
    write(10,42)    ItNum
    write(10,50)
    write(10,43)    NoLoops
    write(10,50)
    write(10,48)    ScalingOMax
    write(10,50)
    write(10,49)    ScalingOMin
    write(10,50)
    write(10,46)    inint
    write(10,50)
    write(10,50)
    do j=1, NumS
        write(10,44)    j
        write(10,50)
        do i=1, size(solMax(size(solMax))%vars)
            write(10,47)    i, solMax(j)%vmax(i)
            write(10,45)    i, solMax(j)%vmin(i)
            write(10,50)
        end do
        write(10,50)
    end do
    write(10,50)

    close(10)

    40 format("Number of Swarms                  |", i15)
    41 format("particles per Swarm               |", i15)
    42 format("Number of iterations              |", i15)
    43 format("Number of attempts                |", i15)
    48 format("Maximum omega scaling             |", 10x,f5.2)
    49 format("Minimum omega scaling             |", 10x,f5.2)
    46 format("Time taken                        |", i15)
    47 format("Maximum velocity for parameter ",i2," |", f15.2)
    45 format("Minimum velocity for parameter ",i2," |", f15.2)
    44 format("Swarm Number                      |", i15)
    50 format(50("-"))
    51 format(20("-"),2x,"STATS",3x,20("-"))


!===============================================================================!
!==============================-----END-----====================================!
!===============================================================================!
end program
