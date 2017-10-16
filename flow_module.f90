!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : flow_module.f90                                                   !
!                                                                              !
!   PURPOSE : Module for SEM inflow generator                                  !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        MODULE flow_module

          IMPLICIT NONE

          INTEGER :: Nt
          REAL(KIND=8) :: dt

          !--------------------------------------------------------------------!
          !                  Interfaces of flow Subroutines                    !
          !--------------------------------------------------------------------!
            !--------Basic Setup for flow solver's properties
            INTERFACE FLOW_SETUP
              MODULE PROCEDURE FLOW_SETUP
            END INTERFACE FLOW_SETUP

        CONTAINS
          SUBROUTINE FLOW_SETUP
            IMPLICIT NONE

            Nt    = 30000
            dt    = 3e-1

          END SUBROUTINE FLOW_SETUP
        END MODULE flow_module
