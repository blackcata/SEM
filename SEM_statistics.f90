!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_statistics.f90                                               !
!                                                                              !
!   PURPOSE : Make a statistics about flows including mean and rms data        !
!                                                                              !
!                                                             2017.03.06 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE SEM_STAT

            USE SEM_module,                                                     &
              ONLY : N, Ny, Nz, Nt, dt, SIGMA, V_b

            USE SEM_module,                                                     &
              ONLY : Y, Z, U, V, W, RS, SEM_EDDY, U_INLET, V_INLET, W_INLET,    &
                     U_COMB, V_COMB, W_COMB

            IMPLICIT NONE
            INTEGER :: j,k

            
        END SUBROUTINE SEM_STAT
