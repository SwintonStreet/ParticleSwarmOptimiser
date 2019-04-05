module Startup
!===============================================================================!
!==============================----START----====================================!
!===============================================================================!

!<<----------------------------------- Module imports
    use    Soln
!<<----------------------------------- End imports
    implicit none

    contains

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME STARTARRAY                                                        --!
!--    Initialises position, velocity of the solutions. In addition it reads  --!
!===============================================================================!

        subroutine STARTARRAY( SolArray, solMin, solMax, varNo, NumS )
            type (sol), dimension(:), intent(inout) :: SolArray, solMin, solMax
            type (sol)                              :: curMin, curMax
            real*8                                  :: input, Curvscale
            integer                                 :: time_array(8), inint, i, j, k, kk
            integer, intent(in)                     :: varNo, NumS

!===============================================================================!
!==========================---break---==========================================! << VALUES INITIALISATION
!===============================================================================!

            call date_and_time(values=time_array)
            inint = (1000 * time_array(7)) + time_array(8)
            input = rand(inint)

            k=1; kk=0


            do i=1, size(SolArray)

                SolArray(i)%bValue = -1

                SolArray(i)%SwarmNo = k

                ! kk is used to ensure we iterate over all the swarms
                if (kk .eq. NumS) then
                    kk = 0; k = k + 1
                end if

                kk = kk + 1

                do j=1, varNo

                    input = (rand() * (solMax(k)%vars(j) - solMin(k)%vars(j)) ) + solMin(k)%vars(j)

                    SolArray(i)%vars(j)  = input
                    SolArray(i)%bvars(j) = input

                    input = solMax(k)%vmax(j) * ((rand() * (solMax(k)%vars(j) - solMin(k)%vars(j)) ) + solMin(k)%vars(j))
                    SolArray(i)%velo(j)  = input
                end do

            end do

        end subroutine

!===============================================================================!
!=======================----Subroutine----======================================!
!--    NAME STARTREAD                                                         --!
!--    Reads in data for the dependent variable and the fitting data.         --!
!===============================================================================!

        subroutine STARTREAD( fileNameD, Darray, fileNameT, Tarray )
            integer                           :: n, stat, stat2
            character*80, intent(in)          :: fileNameD, fileNameT
            real*8, dimension(:), allocatable :: Darray, Tarray
            real*8                            :: dummy

            open(1,file=trim(fileNameD),iostat=stat ,status='old',access='sequential') !Data for fitting
            open(2,file=trim(fileNameT),iostat=stat2,status='old',access='sequential') !Dependent variable

            if(stat .ne. 0 .or. stat2 .ne. 0) then
                print *, "File didn't open!";
                go to 98
            end if

            n=0

            do
                read (1, *,end=97,err=98), dummy    ; n=n+1
            end do

            97 continue

            rewind(1); allocate(Darray(n)); allocate(Tarray(n))

            read (1, 50,err=98), Darray; read (2, 51,err=98), Tarray

            50 format(f3.1); 51 format(f5.2)

            98 close(1); close(2)


            close(1)
            close(2)


        end subroutine


!===============================================================================!
!==============================-----END-----====================================!
!===============================================================================!
end module
