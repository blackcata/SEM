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
              ONLY : N, Ny, Nz, Nt, dt, time

            USE SEM_module,                                                     &
              ONLY : Y, Z, U, V, W, RS, SEM_EDDY, U_INLET, V_INLET, W_INLET,    &
                     U_COMB, V_COMB, W_COMB, U_c, U_pr, rms_pr

            IMPLICIT NONE
            INTEGER :: it,j,k,tt
            REAL(KIND=8) :: U_tmp(3,1:Ny), rms_tmp(4,1:Ny)
            tt = INT(time / dt)

            U_tmp(1:3,1:Ny)   = 0.0
            rms_tmp(1:4,1:Ny) = 0.0

            DO j = 1,Ny
              DO k = 1,Nz
                  U_tmp(1,j) = U_tmp(1,j) + U_COMB(j,k)
                  U_tmp(2,j) = U_tmp(2,j) + V_COMB(j,k)
                  U_tmp(3,j) = U_tmp(3,j) + W_COMB(j,k)

                  rms_tmp(1,j) = rms_tmp(1,j) + (U_COMB(j,k) - U(j,k))**2
                  rms_tmp(2,j) = rms_tmp(2,j) + (V_COMB(j,k) - V(j,k))**2
                  rms_tmp(3,j) = rms_tmp(3,j) + (W_COMB(j,k) - W(j,k))**2
                  rms_tmp(4,j) = rms_tmp(4,j) +                                 &
                                 (U_COMB(j,k) - U(j,k))*(V_COMB(j,k) - V(j,k))
              END DO

              U_tmp(1:3,j) = U_tmp(1:3,j)/Nz
              U_pr(1:3,j)  = ( U_pr(1:3,j) * (tt - 1) + U_tmp(1:3,j) )/tt

              rms_tmp(1:4,j) = rms_tmp(1:4,j)/Nz
              rms_pr(1:4,j) = ( rms_pr(1:4,j) * (tt - 1) + rms_tmp(1:4,j) )/tt

            END DO

        END SUBROUTINE SEM_STAT
