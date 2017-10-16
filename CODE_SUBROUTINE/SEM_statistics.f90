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
              ONLY : Y, Z, U_READ, V_READ, W_READ, T_READ, RS, THS,             &
                     SEM_EDDY, U_INLET, V_INLET, W_INLET,                       &
                     U_COMB, V_COMB, W_COMB, T_COMB, T_INLET, U_c, U_pr, rms_pr

            IMPLICIT NONE
            INTEGER :: it,j,k,tt
            REAL(KIND=8) :: U_tmp(4,1:Ny), rms_tmp(8,1:Ny)
            tt = INT(time / dt)

            U_tmp(1:4,1:Ny)   = 0.0
            rms_tmp(1:8,1:Ny) = 0.0

            !$OMP PARALLEL DO private(k,j)
            DO j = 1,Ny
              DO k = 1,Nz
                  U_tmp(1,j) = U_tmp(1,j) + U_COMB(j,k)
                  U_tmp(2,j) = U_tmp(2,j) + V_COMB(j,k)
                  U_tmp(3,j) = U_tmp(3,j) + W_COMB(j,k)
                  U_tmp(4,j) = U_tmp(4,j) + T_COMB(j,k)

                  rms_tmp(1,j) = rms_tmp(1,j) + (U_COMB(j,k) - U_READ(j,k))**2
                  rms_tmp(2,j) = rms_tmp(2,j) + (V_COMB(j,k) - V_READ(j,k))**2
                  rms_tmp(3,j) = rms_tmp(3,j) + (W_COMB(j,k) - W_READ(j,k))**2
                  rms_tmp(4,j) = rms_tmp(4,j) + (T_COMB(j,k) - T_READ(j,k))**2
                  rms_tmp(5,j) = rms_tmp(5,j) +                                 &
                                 (U_COMB(j,k) - U_READ(j,k))                    &
                                *(V_COMB(j,k) - V_READ(j,k))
                  rms_tmp(6,j) = rms_tmp(6,j) +                                 &
                                 (U_COMB(j,k) - U_READ(j,k))                    &
                                *(T_COMB(j,k) - T_READ(j,k))
                  rms_tmp(7,j) = rms_tmp(7,j) +                                 &
                                 (V_COMB(j,k) - V_READ(j,k))                    &
                                *(T_COMB(j,k) - T_READ(j,k))
                  rms_tmp(8,j) = rms_tmp(8,j) +                                 &
                                 (W_COMB(j,k) - W_READ(j,k))                    &
                                *(T_COMB(j,k) - T_READ(j,k))


              END DO

              U_tmp(1:4,j) = U_tmp(1:4,j)/Nz
              U_pr(1:4,j)  = ( U_pr(1:4,j) * (tt - 1) + U_tmp(1:4,j) )/tt

              rms_tmp(1:8,j) = rms_tmp(1:8,j)/Nz
              rms_pr(1:8,j) = ( rms_pr(1:8,j) * (tt - 1) + rms_tmp(1:8,j) )/tt

            END DO
            !OMP END PARALLEL

        END SUBROUTINE SEM_STAT
