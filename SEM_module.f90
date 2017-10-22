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
!               U,V,W,T_COMB  : Reconstructed components of inflow              !
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
            REAL(KIND=8) :: T_int     ! Eddy's T intensity
          END TYPE EDDY_CHAR

          INTEGER :: N, Ny, Nz, OUT_NUM
          REAL(KIND=8) :: SIGMA, V_b, time, eps

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Y,Z
          REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: U,V,W,T,                   &
                                                     U_INLET,V_INLET,W_INLET,   &
                                                     U_COMB,V_COMB,W_COMB,      &
                                                     T_INLET, T_COMB,           &
                                                     U_pr, rms_pr, U_c
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: RS, THS
          TYPE(EDDY_CHAR),DIMENSION(:),ALLOCATABLE  :: SEM_EDDY

          !--------------------------------------------------------------------!
          !                  Interfaces of SEM Subroutines                     !
          !--------------------------------------------------------------------!

            !--------Intensity Determination function
            INTERFACE INTENSITY_det
              MODULE PROCEDURE INTENSITY_det
            END INTERFACE INTENSITY_det

            !--------Cholesky Decomposition Subrouitne
            INTERFACE CHOL
              MODULE PROCEDURE CHOL
            END INTERFACE CHOL

            !--------Matrix Multiplication Subrouine
            INTERFACE MAT_MUL
              MODULE PROCEDURE MAT_MUL
            END INTERFACE MAT_MUL

            !--------SEM Main
            INTERFACE SEM_main
              MODULE PROCEDURE SEM_main
            END INTERFACE SEM_main

            !--------Setup variables and constants
            INTERFACE SETUP
              MODULE PROCEDURE SETUP
            END INTERFACE SETUP

            !--------Reading datas from external data
            INTERFACE READ_DATA
              MODULE PROCEDURE READ_DATA
            END INTERFACE READ_DATA
            !
            ! !--------Initializing eddies position and intensities
            ! INTERFACE EDDY_SETTING
            !   MODULE PROCEDURE EDDY_SETTING
            ! END INTERFACE EDDY_SETTING
            !
            ! !--------Compose flucutuation components by combining each eddies
            ! INTERFACE FLUCT_GEN
            !   MODULE PROCEDURE FLUCT_GEN
            ! END INTERFACE FLUCT_GEN
            !
            ! !--------Combining given mean & Reynolds stress with flucutuations
            ! INTERFACE COMB_SLICE
            !   MODULE PROCEDURE COMB_SLICE
            ! END INTERFACE COMB_SLICE
            !
            ! !--------Convecting each eddies with periodic boundaery conditions
            ! INTERFACE CONVECT_EDDY
            !   MODULE PROCEDURE CONVECT_EDDY
            ! END INTERFACE CONVECT_EDDY
            !
            ! !--------Making statistics of SEM results
            ! INTERFACE SEM_STAT
            !   MODULE PROCEDURE SEM_STAT
            ! END INTERFACE SEM_STAT

            SAVE

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

!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_main.f90                                                     !
!                                                                              !
!   PURPOSE : To make inflow generation by using SEM (Synthetic Eddy Method)   !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE SEM_main(it)
          IMPLICIT NONE
          INTEGER :: it

          IF ( it == 1 ) THEN
            CALL SETUP
            CALL READ_DATA
            CALL EDDY_SETTING
          END IF

          CALL FLUCT_GEN
          CALL COMB_SLICE
          CALL CONVECT_EDDY
          CALL SEM_STAT

        END SUBROUTINE SEM_main

!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_setup.f90                                                    !
!                                                                              !
!   PURPOSE : Setup for SEM inflow generator                                   !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE SETUP

            USE flow_module

            IMPLICIT NONE

            !------------------------------------------------------------------!
            !                         Constants for SEM                        !
            !------------------------------------------------------------------!
            N  = 1000
            Ny = 65
            Nz = 66

            SIGMA = 13.00

            OUT_NUM = 1

            eps = 1e-8

            !------------------------------------------------------------------!
            !                         Allocate variables                       !
            !------------------------------------------------------------------!
            ALLOCATE( Y(1:Ny),Z(1:Nz) )
            ALLOCATE( U(1:Ny,1:Nz), V(1:Ny,1:Nz), W(1:Ny,1:Nz), T(1:Ny,1:Nz) )
            ALLOCATE( U_INLET(1:Ny,1:Nz),V_INLET(1:Ny,1:Nz),W_INLET(1:Ny,1:Nz) )
            ALLOCATE( U_COMB(1:Ny,1:Nz),V_COMB(1:Ny,1:Nz),W_COMB(1:Ny,1:Nz) )
            ALLOCATE( T_INLET(1:Ny,1:Nz), T_COMB(1:Ny,1:Nz)  )
            ALLOCATE( RS(6,1:Ny,1:Nz), THS(4,1:Ny,1:Nz), U_c(1:Ny,1:Nz) )
            ALLOCATE( SEM_EDDY(1:N), U_pr(4,1:Nz), rms_pr(10,1:Nz) )

            !------------------------------------------------------------------!
            !                         Initial Conditions                       !
            !------------------------------------------------------------------!
            Y(1:Ny) = 0.0
            Z(1:Nz) = 0.0

            U(1:Ny,1:Nz) = 0.0
            V(1:Ny,1:Nz) = 0.0
            W(1:Ny,1:Nz) = 0.0
            T(1:Ny,1:Nz) = 0.0

            RS(1:6,1:Ny,1:Nz) = 0.0
            THS(1:4,1:Ny,1:Nz) = 0.0

            U_INLET(1:Ny,1:Nz) = 0.0
            V_INLET(1:Ny,1:Nz) = 0.0
            W_INLET(1:Ny,1:Nz) = 0.0
            T_INLET(1:Ny,1:Nz) = 0.0

            U_COMB(1:Ny,1:Nz) = 0.0
            V_COMB(1:Ny,1:Nz) = 0.0
            W_COMB(1:Ny,1:Nz) = 0.0
            T_COMB(1:Ny,1:Nz) = 0.0

            U_c(1:Ny,1:Nz)   = 0.0
            U_pr(1:4,1:Ny)   = 0.0
            rms_pr(1:10,1:Ny) = 0.0

            SEM_EDDY(1:N)%eddy_num = 0
            SEM_EDDY(1:N)%eddy_len = 0.0
            SEM_EDDY(1:N)%X_pos    = 0.0
            SEM_EDDY(1:N)%Y_pos    = 0.0
            SEM_EDDY(1:N)%Z_pos    = 0.0
            SEM_EDDY(1:N)%X_int    = 0.0
            SEM_EDDY(1:N)%Y_int    = 0.0
            SEM_EDDY(1:N)%Z_int    = 0.0
            SEM_EDDY(1:N)%T_int    = 0.0

        END SUBROUTINE SETUP
!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_read.f90                                                     !
!                                                                              !
!   PURPOSE : Reading datas from another simulations including mean flow and   !
!             Reynolds stress.                                                 !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE READ_DATA

          USE IO_module

          IMPLICIT NONE

          INTEGER :: j,k
          REAL(KIND=8) :: time_sta, time_end, tmp_y, tmp_z, S
          CHARACTER(20) :: header

          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) '              READING PROCESS STARTED               '
          CALL CPU_TIME(time_sta)

          OPEN(100,FILE=path_name,FORM='FORMATTED',STATUS='OLD')
          READ(100,*) header
          READ(100,*) header

          !--------------------------------------------------------------------!
          !                  Main loop of reading SLICE datas                  !
          !--------------------------------------------------------------------!
          DO j = 1,Ny
            DO k = 1,Nz

                READ(100,*) tmp_z, tmp_y, U(j,k), V(j,k), W(j,k), T(j,k),       &
                            RS(1,j,k), RS(2,j,k), RS(3,j,k), RS(4,j,k),         &
                            RS(5,j,k), RS(6,j,k),                               &
                            THS(1,j,k), THS(2,j,k), THS(3,j,k), THS(4,j,k)

                IF (k==1) Y(j)      = tmp_y
                IF (j==1) Z(k)      = tmp_z

            END DO
          END DO

          CLOSE(100)

          S = ( Y(Ny) - Y(1) + 2*SIGMA ) * ( Z(Nz) - Z(1) + 2*SIGMA )
          V_b = S*2*SIGMA

          CALL CPU_TIME(time_end)
          WRITE(*,*) '            READING PROCESS IS COMPLETED            '
          WRITE(*,*) '  Total Reading time : ',time_end - time_sta,' s'
          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) ''

        END SUBROUTINE READ_DATA

        END MODULE
