!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_main.f90                                                     !
!                                                                              !
!   PURPOSE : To make inflow generation by using SEM (Synthetic Eddy Method)   !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        PROGRAM SEM_main
          USE SEM_module,                                                     &
            ONLY : Nt

          IMPLICIT NONE
          INTEGER :: it
          REAL(KIND=8) :: time_sta, time_end

          CALL SETUP
          CALL READ_DATA
          CALL EDDY_SETTING

          DO it = 1,Nt
            CALL CPU_TIME(time_sta)
            CALL FLUCT_GEN
            CALL COMB_SLICE
            ! CALL CONVECT_EDDY
            CALL CPU_TIME(time_end)

            WRITE(*,"(A,I5,2X,A,F10.6,A)") 'SEM for',it,'iteration time : ',time_end - time_sta,' s'
          END DO

          CALL OUTPUT

        END PROGRAM SEM_main
