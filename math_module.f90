!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : math_module.f90                                                  !
!                                                                              !
!   PURPOSE : math mathical function & subroutines                             !
!                                                                              !
!                                                             2017.10.22 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!
        MODULE math_module
          IMPLICIT NONE
          REAL(KIND=8) :: eps = 1e-8

        !--------------------------------------------------------------------!
        !                  Interfaces of math Subroutines                     !
        !--------------------------------------------------------------------!

          !--------Cholesky Decomposition Subrouitne
          INTERFACE CHOL
            MODULE PROCEDURE CHOL
          END INTERFACE CHOL

          !--------Matrix Multiplication Subrouine
          INTERFACE MAT_MUL
            MODULE PROCEDURE MAT_MUL
          END INTERFACE MAT_MUL

          CONTAINS

        !--------------------------------------------------------------------!
        !                   Cholesky Decomposition Function                  !
        !--------------------------------------------------------------------!
        SUBROUTINE CHOL(A,R,N)
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: N
          REAL(KIND=8),INTENT(OUT) :: A(N,N)
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
                A(i,j) = sqrt(abs(A(i,j)) + eps)

              ELSE

                A(i,j) = R(i,j)
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

      END MODULE math_module
