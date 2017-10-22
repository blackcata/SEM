!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_main.f90                                                     !
!                                                                              !
!   PURPOSE : To make inflow generation by using SEM (Synthetic Eddy Method)   !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        PROGRAM main
          USE SEM_module
          USE flow_module

          IMPLICIT NONE
          INTEGER :: it
          REAL(KIND=8) :: time_sta, time_end

          CALL FLOW_SETUP
          
          DO it = 1,Nt
            time = it * dt
            CALL CPU_TIME(time_sta)

            CALL SEM_main(it)

            CALL CPU_TIME(time_end)

            WRITE(*,"(A,1X,I7,2X,A,F10.6,A)")                                   &
                       'SEM for',it,'iteration time : ',time_end - time_sta,' s'
            IF( mod(it,OUT_NUM) == 0 ) CALL OUTPUT
          END DO

        END PROGRAM main
