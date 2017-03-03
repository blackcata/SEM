!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_write.f90                                                    !
!                                                                              !
!   PURPOSE : Write each variables in the RESULT folder.                       !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!
          SUBROUTINE OUTPUT

            USE SEM_module,                                                     &
              ONLY : Y, Z, U, V, W, RS, U_INLET, V_INLET, W_INLET, SEM_EDDY,    &
                     U_COMB, V_COMB, W_COMB

              IMPLICIT NONE
              REAL(KIND=8) :: time_sta, time_end

              WRITE(*,*) '----------------------------------------------------'
              WRITE(*,*) '              WRITING PROCESS STARTED               '
              CALL CPU_TIME(time_sta)

              CALL CPU_TIME(time_end)

              WRITE(*,*) '           WRITING PROCESS IS COMPLETED            '
              WRITE(*,*) '  Total Writing time : ',time_end - time_sta,' s'
              WRITE(*,*) '----------------------------------------------------'
              WRITE(*,*) ''

              DEALLOCATE(Y,Z,U,V,W,RS,U_INLET,V_INLET,W_INLET,SEM_EDDY)
              DEALLOCATE(U_COMB,V_COMB,W_COMB)

          END SUBROUTINE OUTPUT
