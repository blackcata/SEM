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
              ONLY : Y, Z, U, V, W, T, RS, THS,                                 &
                     SEM_EDDY, U_INLET, V_INLET, W_INLET,                       &
                     U_COMB, V_COMB, W_COMB, T_COMB, T_INLET, U_c, U_pr, rms_pr

            IMPLICIT NONE
            INTEGER :: it,j,k,tt
            REAL(KIND=8) :: U_tmp(4,1:Nz), rms_tmp(10,1:Nz)
            tt = INT(time / dt)

            U_tmp(1:4,1:Nz)   = 0.0
            rms_tmp(1:10,1:Nz) = 0.0

            !$OMP PARALLEL DO private(k,j)
            DO k = 1,Nz
              DO j = 1,Ny
                  U_tmp(1,k) = U_tmp(1,k) + U_COMB(j,k)
                  U_tmp(2,k) = U_tmp(2,k) + V_COMB(j,k)
                  U_tmp(3,k) = U_tmp(3,k) + W_COMB(j,k)
                  U_tmp(4,k) = U_tmp(4,k) + T_COMB(j,k)

                  rms_tmp(1,k) = rms_tmp(1,k) + (U_COMB(j,k) - U(j,k))**2
                  rms_tmp(2,k) = rms_tmp(2,k) + (V_COMB(j,k) - V(j,k))**2
                  rms_tmp(3,k) = rms_tmp(3,k) + (W_COMB(j,k) - W(j,k))**2
                  rms_tmp(4,k) = rms_tmp(4,k) + (T_COMB(j,k) - T(j,k))**2
                  rms_tmp(5,k) = rms_tmp(5,k) +                                 &
                                 (U_COMB(j,k) - U(j,k))*(V_COMB(j,k) - V(j,k))
                  rms_tmp(6,k) = rms_tmp(6,k) +                                 &
                                 (U_COMB(j,k) - U(j,k))*(W_COMB(j,k) - W(j,k))
                  rms_tmp(7,k) = rms_tmp(7,k) +                                 &
                                 (V_COMB(j,k) - V(j,k))*(W_COMB(j,k) - W(j,k))
                  rms_tmp(8,k) = rms_tmp(8,k) +                                 &
                                 (U_COMB(j,k) - U(j,k))*(T_COMB(j,k) - T(j,k))
                  rms_tmp(9,k) = rms_tmp(9,k) +                                 &
                                 (V_COMB(j,k) - V(j,k))*(T_COMB(j,k) - T(j,k))
                  rms_tmp(10,k) = rms_tmp(10,k) +                                 &
                                 (W_COMB(j,k) - W(j,k))*(T_COMB(j,k) - T(j,k))


              END DO

              U_tmp(1:4,k) = U_tmp(1:4,k)/Ny
              U_pr(1:4,k)  = ( U_pr(1:4,k) * (tt - 1) + U_tmp(1:4,k) )/tt

              rms_tmp(1:10,k) = rms_tmp(1:10,k)/Ny
              rms_pr(1:10,k) = ( rms_pr(1:10,k) * (tt - 1) + rms_tmp(1:10,k) )/tt

            END DO
            !OMP END PARALLEL

        END SUBROUTINE SEM_STAT
