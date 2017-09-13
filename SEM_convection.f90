!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_convection.f90                                               !
!                                                                              !
!   PURPOSE : Convect each eddies by convective velocity                       !
!                                                                              !
!                                                             2017.03.03 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE CONVECT_EDDY

            USE SEM_module,                                                     &
              ONLY : N, Ny, Nz, SIGMA, dt

            USE SEM_module,                                                     &
              ONLY : Y, Z, SEM_EDDY, U_COMB, V_COMB, W_COMB, T_COMB

            USE SEM_module,                                                     &
              ONLY : INTENSITY_det

            IMPLICIT NONE
            INTRINSIC :: sqrt, abs

            INTEGER :: it,j,k
            REAL(KIND=8) :: time_sta, time_end, Y_start, Y_end, Z_start, Z_end, &
                            U_conv, V_conv, W_conv, tmp_z
            REAL(KIND=8) :: tmp(1:6)

            ! WRITE(*,*) '----------------------------------------------------'
            ! WRITE(*,*) '             CONVECTION PROCESS STARTED             '
            ! CALL CPU_TIME(time_sta)

            Y_start = Y(1) - SIGMA
            Y_end   = Y(Ny) + SIGMA

            Z_start = Z(1) - SIGMA
            Z_end   = Z(Nz) + SIGMA

            CALL RANDOM_SEED

            U_conv = SUM(U_COMB)/(Ny*Nz)
            V_conv = SUM(V_COMB)/(Ny*Nz)
            W_conv = SUM(W_COMB)/(Ny*Nz)

            DO it = 1,N
              SEM_EDDY(it)%X_pos = SEM_EDDY(it)%X_pos + U_conv*dt
              SEM_EDDY(it)%Y_pos = SEM_EDDY(it)%Y_pos + V_conv*dt
              SEM_EDDY(it)%Z_pos = SEM_EDDY(it)%Z_pos + W_conv*dt

              IF ( (SEM_EDDY(it)%X_pos-(-SIGMA))*                               &
                   (SEM_EDDY(it)%X_pos-(SIGMA)) > 0 .OR.                        &
                   (SEM_EDDY(it)%Y_pos - Y_start)*                              &
                   (SEM_EDDY(it)%Y_pos - Y_end) > 0 ) THEN

                   SEM_EDDY(it)%eddy_len = SIGMA
                   SEM_EDDY(it)%X_pos = - SIGMA

                   CALL RANDOM_NUMBER(tmp)
                   SEM_EDDY(it)%Y_pos = Y_start + (Y_end-Y_start)*tmp(1)
                   SEM_EDDY(it)%Z_pos = Z_start + (Z_end-Z_start)*tmp(2)

                   SEM_EDDY(it)%X_int = INTENSITY_det(tmp(3)-0.5)
                   SEM_EDDY(it)%Y_int = INTENSITY_det(tmp(4)-0.5)
                   SEM_EDDY(it)%Z_int = INTENSITY_det(tmp(5)-0.5)
                   SEM_EDDY(it)%T_int = INTENSITY_det(tmp(6)-0.5)

              END IF
              !----------------------------------------------------------------!
              !                 Periodic boundary conditions                   !
              !----------------------------------------------------------------!
              IF ( SEM_EDDY(it)%Z_pos < Z_start ) THEN
                tmp_z = Z_start - SEM_EDDY(it)%Z_pos
                SEM_EDDY(it)%Z_pos = Z_end - tmp_z
              END IF

              IF ( SEM_EDDY(it)%Z_pos > Z_end ) THEN
                tmp_z = SEM_EDDY(it)%Z_pos - Z_end
                SEM_EDDY(it)%Z_pos = Z_start + tmp_z
              END IF

            END DO

            ! CALL CPU_TIME(time_end)
            ! WRITE(*,*) '           CONVECTION PROCESS IS COMPLETED          '
            ! WRITE(*,*) '  Total Reading time : ',time_end - time_sta,' s'
            ! WRITE(*,*) '----------------------------------------------------'
            ! WRITE(*,*) ''

        END SUBROUTINE CONVECT_EDDY
