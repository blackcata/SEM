!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_module.f90                                                   !
!                                                                              !
!   PURPOSE : Module for SEM inflow generator                                  !
!                                                                              !
!                                                             2016.03.02 K.Noh !
!                                                                              !
!   VARIABLES : dt    : Time step                                              !
!               N     : The number of eddies                                   !
!               SIGMA : Eddy length scale                                      !
!               V_b   : Volume of box including eddies                         !
!               Nt    : The number of iterations
!                                                                              !
!               Y,Z   : Y,Z coordinates                                        !
!               U,V,W : Mean velocity arrays                                   !
!               RS    : Reynolds stress                                        !
!                                                                              !
!------------------------------------------------------------------------------!

        MODULE SEM_module

          IMPLICIT NONE

          TYPE EDDY_CHAR
            INTEGER :: eddy_num
            REAL(KIND=8) :: eddy_len
            REAL(KIND=8) :: X_pos
            REAL(KIND=8) :: Y_pos
            REAL(KIND=8) :: Z_pos
            REAL(KIND=8) :: X_int
            REAL(KIND=8) :: Y_int
            REAL(KIND=8) :: Z_int
          END TYPE EDDY_CHAR

          INTEGER :: N, Ny, Nz, Nt
          REAL(KIND=8) :: dt, SIGMA, V_b
          CHARACTER(LEN=65) :: file_name, dir_name, path_name

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Y,Z
          REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: U,V,W,                     &
                                                     U_INLET,V_INLET,W_INLET,   &
                                                     U_COMB,V_COMB,W_COMB
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: RS
          TYPE(EDDY_CHAR),DIMENSION(:),ALLOCATABLE :: SEM_EDDY

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
            INTEGER :: N
            REAL(KIND=8) :: A(N,N)
            REAL(KIND=8),INTENT(IN) :: R(N,N)

            INTEGER :: i,j,k
            A(1:N,1:N) = 0.0

            DO i = 1,N
              DO j = 1,i
                IF (i==j) THEN

                  A(i,j) = R(i,j)
                  DO k = 1,j-1
                    A(i,j) = A(i,j) - A(j,k)**2
                  END DO
                  A(i,j) = sqrt(A(i,j))
                  print*,i,j,A(i,j)
                ELSE

                  A(i,j) = R(i,j)
                  DO k = 1,j-1
                    A(i,j) = A(i,j) - A(i,k)*A(j,k)
                  END DO
                  A(i,j) = A(i,j)/A(j,j)
                  print*,i,j,A(i,j)
                END IF
              END DO
            END DO

          END SUBROUTINE

        END MODULE
