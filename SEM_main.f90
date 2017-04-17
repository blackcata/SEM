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
          USE SEM_module,                                                       &
            ONLY : Nt, time, dt, OUT_NUM

          IMPLICIT NONE
          INTEGER :: it
          REAL(KIND=8) :: time_sta, time_end

          CALL SETUP
          CALL READ_DATA
          CALL EDDY_SETTING

          DO it = 1,Nt
            time = it * dt
            CALL CPU_TIME(time_sta)
            CALL FLUCT_GEN
            CALL COMB_SLICE
            CALL CONVECT_EDDY
            CALL SEM_STAT
            CALL CPU_TIME(time_end)

            WRITE(*,"(A,I5,2X,A,F10.6,A)")                                      &
                       'SEM for',it,'iteration time : ',time_end - time_sta,' s'
            IF( mod(it,OUT_NUM) == 0 ) CALL OUTPUT
          END DO

        END PROGRAM SEM_main
