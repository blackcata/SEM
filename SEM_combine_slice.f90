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
              ONLY : U_INLET, V_INLET, W_INLET, SEM_EDDY, U, V, W, T, RS, THS,  &
                     U_COMB, V_COMB, W_COMB, T_INLET, T_COMB

            USE SEM_module,                                                     &
             ONLY : CHOL,MAT_MUL

            IMPLICIT NONE
            INTRINSIC :: sqrt, abs

            INTEGER :: j,k
            REAL(KIND=8) :: R_loc(4,4), A(4,4),             &
                            u_ins(4,1), u_mean(4,1), u_fluc(4,1), u_tmp(4,1)

            DO k = 1,Nz
              DO j = 1,Ny
                A(1:4,1:4)     = 0.0
                R_loc(1:4,1:4) = 0.0
                u_ins(1:4,1)   = 0.0
                u_fluc(1:4,1)  = 0.0
                u_mean(1:4,1)  = 0.0
                u_tmp(1:4,1)   = 0.0

                U_mean(1:4,1) = (/U(j,k),V(j,k),W(j,k),T(j,k)/)
                u_tmp(1:4,1)  = (/U_INLET(j,k),V_INLET(j,k),                    &
                                  W_INLET(j,k),T_INLET(j,k)/)
                R_loc(1,1:4)  = (/RS(1,j,k),RS(4,j,k),RS(5,j,k),THS(2,j,k)/)
                R_loc(2,1:4)  = (/RS(4,j,k),RS(2,j,k),RS(6,j,k),THS(3,j,k)/)
                R_loc(3,1:4)  = (/RS(5,j,k),RS(6,j,k),RS(3,j,k),THS(4,j,k)/)
                R_loc(4,1:4)  = (/THS(2,j,k),THS(3,j,k),THS(4,j,k),THS(1,j,k)/)

                CALL CHOL(A,R_loc,4)
                CALL MAT_MUL(A,u_tmp,u_fluc,4,1,4)

                u_ins(1:4,1) = U_mean(1:4,1) + u_fluc(1:4,1)

                U_COMB(j,k) = u_ins(1,1)
                V_COMB(j,k) = u_ins(2,1)
                W_COMB(j,k) = u_ins(3,1)
                T_COMB(j,k) = u_ins(4,1)

              END DO
            END DO

        END SUBROUTINE COMB_SLICE
