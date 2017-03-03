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

          IMPLICIT NONE

          CALL SETUP
          CALL READ_DATA

          CALL EDDY_SETTING
          CALL FLUCT_GEN
          CALL COMB_SLICE
          CALL CONVECT_EDDY

          CALL OUTPUT

        END PROGRAM SEM_main
