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

            INTEGER :: it,j,k, tt
            REAL(KIND=8) :: Y_start, Y_end, Z_start, Z_end, &
                            U_conv, V_conv, W_conv, tmp_y
            REAL(KIND=8) :: tmp(1:6)

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
                   (SEM_EDDY(it)%Z_pos - Z_start)*                              &
                   (SEM_EDDY(it)%Z_pos - Z_end) > 0 ) THEN

                   SEM_EDDY(it)%eddy_len = SIGMA
                   SEM_EDDY(it)%X_pos = - SIGMA

                   CALL RANDOM_NUMBER(tmp)
                   SEM_EDDY(it)%X_pos = 0 - SIGMA
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
              IF ( SEM_EDDY(it)%Y_pos < Y_start ) THEN
                tmp_y = Y_start - SEM_EDDY(it)%Z_pos
                SEM_EDDY(it)%Y_pos = Y_end - tmp_y
              END IF

              IF ( SEM_EDDY(it)%Y_pos > Y_end ) THEN
                tmp_y = SEM_EDDY(it)%Z_pos - Y_end
                SEM_EDDY(it)%Y_pos = Y_start + tmp_y
              END IF

            END DO

        END SUBROUTINE CONVECT_EDDY
