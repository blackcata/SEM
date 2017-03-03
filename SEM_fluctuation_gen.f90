!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_fluctuation_gen.f90                                          !
!                                                                              !
!   PURPOSE : Generate fluctuations without combining mean and rms data        !
!                                                                              !
!                                                             2017.03.03 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE FLUCT_GEN

            USE SEM_module,                                                     &
              ONLY : N, Ny, Nz, SIGMA, V_b

            USE SEM_module,                                                     &
              ONLY : Y, Z, U_INLET, V_INLET, W_INLET, SEM_EDDY

            IMPLICIT NONE
            INTRINSIC :: sqrt, abs

            INTEGER :: it,j,k
            REAL(KIND=8) :: time_sta, time_end, x0, y0, z0, f

            WRITE(*,*) '----------------------------------------------------'
            WRITE(*,*) '             FLUCTUATION PROCESS STARTED            '
            CALL CPU_TIME(time_sta)

            DO k = 1,Nz
              DO j = 1,Ny

                DO it = 1,N
                  x0 = (0    - SEM_EDDY(it)%X_pos)/SEM_EDDY(it)%eddy_len
                  y0 = (Y(j) - SEM_EDDY(it)%Y_pos)/SEM_EDDY(it)%eddy_len
                  z0 = (Z(k) - SEM_EDDY(it)%Z_pos)/SEM_EDDY(it)%eddy_len

                  !------------------------------------------------------------!
                  !                        Shape function                      !
                  !------------------------------------------------------------!
                  IF ( abs(x0) <=1 .AND. abs(y0) <=1 .AND. abs(z0) <=1) THEN
                    f = sqrt(1.5) * (1- abs(x0)) *                              &
                        sqrt(1.5) * (1- abs(y0)) *                              &
                        sqrt(1.5) * (1- abs(z0))

                    U_INLET(j,k) = U_INLET(j,k) +                               &
                                   sqrt(V_b/SEM_EDDY(it)%eddy_len**3) *         &
                                   SEM_EDDY(it)%X_int*f

                    V_INLET(j,k) = V_INLET(j,k) +                               &
                                  sqrt(V_b/SEM_EDDY(it)%eddy_len**3) *          &
                                  SEM_EDDY(it)%Y_int*f

                    W_INLET(j,k) = W_INLET(j,k) +                               &
                                  sqrt(V_b/SEM_EDDY(it)%eddy_len**3) *          &
                                  SEM_EDDY(it)%Z_int*f
                  END IF
                END DO

              END DO
            END DO

            U_INLET(1:Ny,1:Nz) = U_INLET(1:Ny,1:Nz) / sqrt(REAL(N,8))
            V_INLET(1:Ny,1:Nz) = V_INLET(1:Ny,1:Nz) / sqrt(REAL(N,8))
            W_INLET(1:Ny,1:Nz) = W_INLET(1:Ny,1:Nz) / sqrt(REAL(N,8))

            CALL CPU_TIME(time_end)
            WRITE(*,*) '         FLUCTUATION PROCESS IS COMPLETED           '
            WRITE(*,*) '  Total Reading time : ',time_end - time_sta,' s'
            WRITE(*,*) '----------------------------------------------------'
            WRITE(*,*) ''

        END SUBROUTINE FLUCT_GEN
