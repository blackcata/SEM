!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_eddy_setting.f90                                             !
!                                                                              !
!   PURPOSE : Setup each eddys characteristic including eddy length,           !
!             positions and intensities of each diretions.                     !
!                                                                              !
!                                                             2017.03.03 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE EDDY_SETTING

            USE SEM_module,                                                     &
              ONLY : N, Ny, Nz, SIGMA

            USE SEM_module,                                                     &
              ONLY : Y, Z, SEM_EDDY

            USE SEM_module,                                                     &
              ONLY : INTENSITY_det

            IMPLICIT NONE

            INTEGER :: it
            REAL(KIND=8) :: Y_start, Y_end, Z_start, Z_end, tmp ,               &
                            time_sta, time_end
            REAL(KIND=8) :: INT_X(1:N), INT_Y(1:N), INT_Z(1:N)

            WRITE(*,*) '----------------------------------------------------'
            WRITE(*,*) '            EDDY SETTING PROCESS STARTED            '
            CALL CPU_TIME(time_sta)


            !------------------------------------------------------------------!
            !                  Make & Initialize Result folder                 !
            !------------------------------------------------------------------!
            Y_start = Y(1) - SIGMA
            Y_end   = Y(Ny) + SIGMA

            Z_start = Z(1) - SIGMA
            Z_end   = Z(Nz) + SIGMA

            CALL RANDOM_SEED
            CALL RANDOM_NUMBER(INT_X)
            CALL RANDOM_NUMBER(INT_Y)
            CALL RANDOM_NUMBER(INT_Z)

            DO it = 1,N
              SEM_EDDY(it)%eddy_num = it
              SEM_EDDY(it)%eddy_len = SIGMA

              CALL RANDOM_NUMBER(tmp)
              SEM_EDDY(it)%X_pos = -SIGMA + 2*SIGMA*tmp

              CALL RANDOM_NUMBER(tmp)
              SEM_EDDY(it)%Y_pos = Y_start + (Y_end-Y_start)*tmp

              CALL RANDOM_NUMBER(tmp)
              SEM_EDDY(it)%Z_pos = Z_start + (Z_end-Z_start)*tmp

              SEM_EDDY(it)%X_int = INTENSITY_det(INT_X(it)-0.5)
              SEM_EDDY(it)%Y_int = INTENSITY_det(INT_Y(it)-0.5)
              SEM_EDDY(it)%Z_int = INTENSITY_det(INT_Z(it)-0.5)
            END DO

            CALL CPU_TIME(time_end)
            WRITE(*,*) '        EDDY SETTING PROCESS IS COMPLETED           '
            WRITE(*,*) '  Total Reading time : ',time_end - time_sta,' s'
            WRITE(*,*) '----------------------------------------------------'
            WRITE(*,*) ''
        END SUBROUTINE EDDY_SETTING
