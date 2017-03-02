!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_main.f90                                                     !
!                                                                              !
!   PURPOSE : To make inflow generation by using SEM (Synthetic Eddy Method)   !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        PROGRAM LES_FILTERING

          IMPLICIT NONE

          CALL SETUP
          CALL READ_DATA

          CALL OUTPUT

        END PROGRAM LES_FILTERING
