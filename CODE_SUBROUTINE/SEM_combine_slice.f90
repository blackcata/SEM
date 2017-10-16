!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_combine_slice.f90                                            !
!                                                                              !
!   PURPOSE : Combine slice mean,rms data with generated fluctuation variables !
!                                                                              !
!                                                             2017.03.03 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE COMB_SLICE

            USE SEM_module,                                                     &
              ONLY : Ny, Nz

            USE SEM_module,                                                     &
              ONLY : U_INLET, V_INLET, W_INLET, SEM_EDDY,                       &
                     U_READ, V_READ, W_READ, T_READ, RS, THS,                   &
                     U_COMB, V_COMB, W_COMB, T_INLET, T_COMB

            USE SEM_module,                                                     &
             ONLY : CHOL,MAT_MUL

            IMPLICIT NONE
            INTRINSIC :: sqrt, abs

            INTEGER :: it,j,k
            REAL(KIND=8) :: time_sta, time_end, R_loc(4,4), A(4,4),             &
                            u_ins(4,1), u_mean(4,1), u_fluc(4,1), u_tmp(4,1)

            ! WRITE(*,*) '----------------------------------------------------'
            ! WRITE(*,*) '               COMBINE PROCESS STARTED              '
            ! CALL CPU_TIME(time_sta)

            DO k = 1,Nz
              DO j = 1,Ny
                R_loc(1:4,1:4) = 0.0
                u_ins(1:4,1)   = 0.0
                u_fluc(1:4,1)  = 0.0

                U_mean(1:4,1) = (/U_READ(j,k),V_READ(j,k),W_READ(j,k),T_READ(j,k)/)
                u_tmp(1:4,1)  = (/U_INLET(j,k),V_INLET(j,k),                    &
                                  W_INLET(j,k),T_INLET(j,k)/)
                R_loc(1,1:4)  = (/RS(1,j,k),RS(4,j,k),RS(5,j,k),THS(2,j,k)/)
                R_loc(2,1:4)  = (/RS(4,j,k),RS(2,j,k),RS(6,j,k),THS(3,j,k)/)
                R_loc(3,1:4)  = (/RS(5,j,k),RS(6,j,k),RS(3,j,k),THS(4,j,k)/)
                R_loc(4,1:4)  = (/THS(2,j,k),THS(3,j,k),THS(4,j,k),THS(1,j,k)/)

                CALL CHOL(A,R_loc,4)
                CALL MAT_MUL(A,u_tmp,u_fluc,4,1,4)

                u_ins(1:4,1) = U_mean(1:4,1) + u_fluc(1:4,1)
WRITE(*,"(2I,3F15.10)"),k,j,u_ins(1,1),u_mean(1,1),u_fluc(1,1)
                U_COMB(j,k) = u_ins(1,1)
                V_COMB(j,k) = u_ins(2,1)
                W_COMB(j,k) = u_ins(3,1)
                T_COMB(j,k) = u_ins(4,1)

              END DO
            END DO

            ! CALL CPU_TIME(time_end)
            ! WRITE(*,*) '           COMBINE PROCESS IS COMPLETED             '
            ! WRITE(*,*) '  Total Reading time : ',time_end - time_sta,' s'
            ! WRITE(*,*) '----------------------------------------------------'
            ! WRITE(*,*) ''

        END SUBROUTINE COMB_SLICE
