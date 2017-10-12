!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_module.f90                                                   !
!                                                                              !
!   PURPOSE : Module for SEM inflow generator                                  !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!   VARIABLES : dt    : Time step                                              !
!               N     : The number of eddies                                   !
!               SIGMA : Eddy length scale                                      !
!               V_b   : Volume of box including eddies                         !
!               Nt    : The number of iterations                               !
!                                                                              !
!               Y,Z      : Y,Z coordinates                                     !
!               U,V,W    : Mean velocity arrays                                !
!               T        : Mean temperature arrays                             !
!               RS       : Reynolds stress                                     !
!               SEM_EDDY : Each eddies properties including positions,         !
!                          intensities, length scales.                         !
!               U,V,W,T_INLET : Stochastic components of inflow surface        !
!               U,V,W,T_COMB  : Reconstructed compoents of inflow              !
!                                                                              !
!               U_c    : Local convection velocities                           !
!               U_pr   : Mean profiles (U,V,W,T)                               !
!               rms_pr : Reynolds stress profiles (uu,vv,ww,tt,uv,ut,vt,wt)    !
!                                                                              !
!------------------------------------------------------------------------------!

        MODULE SEM_module

          IMPLICIT NONE

          TYPE EDDY_CHAR
            INTEGER :: eddy_num       ! Eddy specification number
            REAL(KIND=8) :: eddy_len  ! Eddy length scale
            REAL(KIND=8) :: X_pos     ! Eddy's X position
            REAL(KIND=8) :: Y_pos     ! Eddy's Y position
            REAL(KIND=8) :: Z_pos     ! Eddy's Z position
            REAL(KIND=8) :: X_int     ! Eddy's X intensity
            REAL(KIND=8) :: Y_int     ! Eddy's Y intensity
            REAL(KIND=8) :: Z_int     ! Eddy's Z intensity
            REAL(KIND=8) :: T_int     ! Eddy's Z intensity
          END TYPE EDDY_CHAR

          INTEGER :: N, Ny, Nz, Nt, OUT_NUM
          REAL(KIND=8) :: dt, SIGMA, V_b, time, eps
          CHARACTER(LEN=65) :: file_name, dir_name, path_name

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Y,Z
          REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: U,V,W,T,                   &
                                                     U_INLET,V_INLET,W_INLET,   &
                                                     U_COMB,V_COMB,W_COMB,      &
                                                     T_INLET, T_COMB,           &
                                                     U_pr, rms_pr, U_c
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: RS, THS
          TYPE(EDDY_CHAR),DIMENSION(:),ALLOCATABLE  :: SEM_EDDY

        CONTAINS
          !--------------------------------------------------------------------!
          !                  Intensity determination Function                  !
          !--------------------------------------------------------------------!
          FUNCTION INTENSITY_det(x_int)
            REAL(KIND=8) :: INTENSITY_det
            REAL(KIND=8),INTENT(IN) :: x_int

            IF ( x_int > 0 ) THEN
              INTENSITY_det = 1
            ELSE
              INTENSITY_det = -1
            END IF
          END FUNCTION INTENSITY_det

          !--------------------------------------------------------------------!
          !                   Cholesky Decomposition Function                  !
          !--------------------------------------------------------------------!
          SUBROUTINE CHOL(A,R,N)
            IMPLICIT NONE

            INTEGER,INTENT(IN) :: N
            REAL(KIND=8) :: A(N,N)
            REAL(KIND=8),INTENT(IN) :: R(N,N)

            INTEGER :: i,j,k
            A(1:N,1:N) = 0.0

            DO i = 1,N
              DO j = 1,i
                A(i,j) = R(i,j)

                IF (i==j) THEN
                  DO k = 1,j-1
                    A(i,j) = A(i,j) - A(j,k)**2
                  END DO
                  A(i,j) = sqrt(abs(A(i,j)))

                ELSE
                  DO k = 1,j-1
                    A(i,j) = A(i,j) - A(i,k)*A(j,k)
                  END DO
                  A(i,j) = A(i,j)/(A(j,j) + eps)

                END IF

              END DO
            END DO

          END SUBROUTINE

          !--------------------------------------------------------------------!
          !                   Matrix multiplication function                   !
          !--------------------------------------------------------------------!
          SUBROUTINE MAT_MUL(A,B,AB,Ni,Nj,Nk)
            IMPLICIT NONE

            INTEGER,INTENT(IN) :: Ni,Nj,Nk
            REAL(KIND=8) :: AB(Ni,Nj)
            REAL(KIND=8),INTENT(IN) :: A(Ni,Nk),B(Nk,Nj)

            INTEGER :: i,j,k
            AB(1:Ni,1:Nj) = 0.0

            DO i = 1,Ni
              DO j = 1,Nj
                DO k = 1,Nk
                  AB(i,j) = AB(i,j) + A(i,k)*B(k,j)
                END DO
              END DO
            END DO

          END SUBROUTINE
        END MODULE
