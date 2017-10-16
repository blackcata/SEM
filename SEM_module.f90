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
!               Y,Z                     : Y,Z coordinates                      !
!               U_READ,V_READ,W_READ    : Mean velocity arrays                 !
!               T_READ                  : Mean temperature arrays              !
!               RS                      : Reynolds stress                      !
!               THS                     : Covarainces velocity & temperature   !
!                                                                              !
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

          INTEGER :: N, Ny, Nz, OUT_NUM
          REAL(KIND=8) :: SIGMA, V_b, time, eps

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Y,Z
          REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: U_READ,V_READ,W_READ,T_READ,&
                                                     U_INLET,V_INLET,W_INLET,    &
                                                     U_COMB,V_COMB,W_COMB,       &
                                                     T_INLET, T_COMB,            &
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
          SAVE

          !--------Setup variables and constants
          INTERFACE SETUP
            MODULE PROCEDURE SETUP
          END INTERFACE SETUP

          !--------Reading datas from external data
          INTERFACE READ_DATA
            MODULE PROCEDURE READ_DATA
          END INTERFACE READ_DATA

          !--------Initializing eddies position and intensities
          INTERFACE EDDY_SETTING
            MODULE PROCEDURE EDDY_SETTING
          END INTERFACE EDDY_SETTING

          !--------Compose flucutuation components by combining each eddies
          INTERFACE FLUCT_GEN
            MODULE PROCEDURE FLUCT_GEN
          END INTERFACE FLUCT_GEN

          !--------Combining given mean & Reynolds stress with flucutuations
          INTERFACE COMB_SLICE
            MODULE PROCEDURE COMB_SLICE
          END INTERFACE COMB_SLICE

          !--------Convecting each eddies with periodic boundaery conditions
          INTERFACE CONVECT_EDDY
            MODULE PROCEDURE CONVECT_EDDY
          END INTERFACE CONVECT_EDDY

          !--------Making statistics of SEM results
          INTERFACE SEM_STAT
            MODULE PROCEDURE SEM_STAT
          END INTERFACE SEM_STAT

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

          END SUBROUTINE CHOL

          !--------------------------------------------------------------------!
          !                   Matrix multiplication function                   !
          !--------------------------------------------------------------------!
          SUBROUTINE MAT_MUL(A,B,AB,Ni,Nj,Nk)
            IMPLICIT NONE

            INTEGER,INTENT(IN) :: Ni,Nj,Nk
            REAL(KIND=8),INTENT(OUT) :: AB(Ni,Nj)
            REAL(KIND=8),INTENT(IN) :: A(Ni,Nk),B(Nk,Nj)

            INTEGER :: i,j,k
            AB(1:Ni,1:Nj) = 0.0

            DO i = 1,Ni
              DO j = 1,Nj
                DO k = 1,Nk
                    AB(i,j) = A(i,k)*B(k,j)
                END DO
              END DO
            END DO

          END SUBROUTINE MAT_MUL

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

            IMPLICIT NONE

            !------------------------------------------------------------------!
            !                         Constants for SEM                        !
            !------------------------------------------------------------------!
            N  = 1000
            Ny = 64
            Nz = 64

            SIGMA = 2.00

            OUT_NUM = 1

            eps = 1e-8

            !------------------------------------------------------------------!
            !                         Allocate variables                       !
            !------------------------------------------------------------------!
            ALLOCATE( Y(1:Ny),Z(1:Nz) )
            ALLOCATE( U_READ(1:Ny,1:Nz), V_READ(1:Ny,1:Nz),                     &
                      W_READ(1:Ny,1:Nz), T_READ(1:Ny,1:Nz) )
            ALLOCATE( U_INLET(1:Ny,1:Nz),V_INLET(1:Ny,1:Nz),W_INLET(1:Ny,1:Nz) )
            ALLOCATE( U_COMB(1:Ny,1:Nz),V_COMB(1:Ny,1:Nz),W_COMB(1:Ny,1:Nz) )
            ALLOCATE( T_INLET(1:Ny,1:Nz), T_COMB(1:Ny,1:Nz)  )
            ALLOCATE( RS(6,1:Ny,1:Nz), THS(4,1:Ny,1:Nz), U_c(1:Ny,1:Nz) )
            ALLOCATE( SEM_EDDY(1:N), U_pr(4,1:Ny), rms_pr(10,1:Ny) )

            !------------------------------------------------------------------!
            !                         Initial Conditions                       !
            !------------------------------------------------------------------!
            Y(1:Ny) = 0.0
            Z(1:Nz) = 0.0

            U_READ(1:Ny,1:Nz) = 0.0
            V_READ(1:Ny,1:Nz) = 0.0
            W_READ(1:Ny,1:Nz) = 0.0
            T_READ(1:Ny,1:Nz) = 0.0

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

          dir_name = 'DATA'
          file_name = 'Slice_TBL.plt'
          path_name = TRIM(dir_name)//'/'//TRIM(file_name)

          OPEN(100,FILE=path_name,FORM='FORMATTED',STATUS='OLD')
          READ(100,*) header
          READ(100,*) header

          !--------------------------------------------------------------------!
          !                  Main loop of reading SLICE datas                  !
          !--------------------------------------------------------------------!
          DO j = 1,Ny
            DO k = 1,Nz

                READ(100,*) tmp_z, tmp_y, U_READ(j,k), V_READ(j,k),             &
                            W_READ(j,k), T_READ(j,k), RS(1,j,k), RS(2,j,k),     &
                            RS(3,j,k), RS(4,j,k), RS(5,j,k), RS(6,j,k),         &
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

            IMPLICIT NONE

            INTEGER :: it
            REAL(KIND=8) :: Y_start, Y_end, Z_start, Z_end,                     &
                            time_sta, time_end
            REAL(KIND=8) :: INT_X(1:N), INT_Y(1:N), INT_Z(1:N), tmp(1:3)

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
              SEM_EDDY(it)%X_pos = -SIGMA + 2*SIGMA*tmp(1)
              SEM_EDDY(it)%Y_pos = Y_start + (Y_end-Y_start)*tmp(2)
              SEM_EDDY(it)%Z_pos = Z_start + (Z_end-Z_start)*tmp(3)

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

            IMPLICIT NONE
            INTRINSIC :: sqrt, abs

            INTEGER :: it,j,k
            REAL(KIND=8) :: x0, y0, z0, f

            !$OMP PARALLEL DO private(k,j,it,x0,y0,z0,f)
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

                    T_INLET(j,k) = T_INLET(j,k) +                               &
                                  sqrt(V_b/SEM_EDDY(it)%eddy_len**3) *          &
                                  SEM_EDDY(it)%Z_int*f
                  END IF
                END DO

              END DO
            END DO
            !OMP END PARALLEL

            U_INLET(1:Ny,1:Nz) = U_INLET(1:Ny,1:Nz) / sqrt(REAL(N,8))
            V_INLET(1:Ny,1:Nz) = V_INLET(1:Ny,1:Nz) / sqrt(REAL(N,8))
            W_INLET(1:Ny,1:Nz) = W_INLET(1:Ny,1:Nz) / sqrt(REAL(N,8))
            T_INLET(1:Ny,1:Nz) = T_INLET(1:Ny,1:Nz) / sqrt(REAL(N,8))

        END SUBROUTINE FLUCT_GEN

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

            IMPLICIT NONE
            INTRINSIC :: sqrt, abs

            INTEGER :: j,k,ii,jj,kk
            REAL(KIND=8) :: R_loc(4,4), A(4,4),                                 &
                            u_ins(4,1), u_mean(4,1), u_fluc(4,1), u_tmp(4,1)

            DO k = 1,Nz
              DO j = 1,Ny
                A(1:4,1:4)     = 0.0
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
                ! CALL MAT_MUL(A,u_tmp,u_fluc,4,1,4)
                DO ii = 1,4
                  DO jj = 1,1
                    DO kk = 1,4
                        u_fluc(ii,jj) = u_fluc(ii,jj) + A(ii,kk) * u_tmp(kk,jj)
                    END DO
                  END DO
                END DO

                u_ins(1:4,1) = U_mean(1:4,1) + u_fluc(1:4,1)

                U_COMB(j,k) = u_ins(1,1)
                V_COMB(j,k) = u_ins(2,1)
                W_COMB(j,k) = u_ins(3,1)
                T_COMB(j,k) = u_ins(4,1)

              END DO
            END DO

        END SUBROUTINE COMB_SLICE

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

            USE flow_module
            IMPLICIT NONE

            INTRINSIC :: sqrt, abs

            INTEGER :: it
            REAL(KIND=8) :: Y_start, Y_end, Z_start, Z_end, &
                            U_conv, V_conv, W_conv, tmp_z
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
                   (SEM_EDDY(it)%Y_pos - Y_start)*                              &
                   (SEM_EDDY(it)%Y_pos - Y_end) > 0 ) THEN

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
              IF ( SEM_EDDY(it)%Z_pos < Z_start ) THEN
                tmp_z = Z_start - SEM_EDDY(it)%Z_pos
                SEM_EDDY(it)%Z_pos = Z_end - tmp_z
              END IF

              IF ( SEM_EDDY(it)%Z_pos > Z_end ) THEN
                tmp_z = SEM_EDDY(it)%Z_pos - Z_end
                SEM_EDDY(it)%Z_pos = Z_start + tmp_z
              END IF

            END DO

        END SUBROUTINE CONVECT_EDDY

!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_statistics.f90                                               !
!                                                                              !
!   PURPOSE : Make a statistics about flows including mean and rms data        !
!                                                                              !
!                                                             2017.03.06 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE SEM_STAT
            USE flow_module

            IMPLICIT NONE
            INTEGER :: j,k,tt
            REAL(KIND=8) :: U_tmp(4,1:Ny), rms_tmp(10,1:Ny)
            tt = INT(time / dt)

            U_tmp(1:4,1:Ny)   = 0.0
            rms_tmp(1:10,1:Ny) = 0.0

            !$OMP PARALLEL DO private(k,j)
            DO j = 1,Ny
              DO k = 1,Nz
                  U_tmp(1,j) = U_tmp(1,j) + U_COMB(j,k)
                  U_tmp(2,j) = U_tmp(2,j) + V_COMB(j,k)
                  U_tmp(3,j) = U_tmp(3,j) + W_COMB(j,k)
                  U_tmp(4,j) = U_tmp(4,j) + T_COMB(j,k)

                  rms_tmp(1,j) = rms_tmp(1,j) + (U_COMB(j,k) - U_READ(j,k))**2
                  rms_tmp(2,j) = rms_tmp(2,j) + (V_COMB(j,k) - V_READ(j,k))**2
                  rms_tmp(3,j) = rms_tmp(3,j) + (W_COMB(j,k) - W_READ(j,k))**2
                  rms_tmp(4,j) = rms_tmp(4,j) + (T_COMB(j,k) - T_READ(j,k))**2

                  rms_tmp(5,j) = rms_tmp(5,j) +                                 &
                                 (U_COMB(j,k) - U_READ(j,k))                    &
                                *(V_COMB(j,k) - V_READ(j,k))
                  rms_tmp(6,j) = rms_tmp(6,j) +                                 &
                                 (U_COMB(j,k) - U_READ(j,k))                    &
                                *(W_COMB(j,k) - W_READ(j,k))
                  rms_tmp(7,j) = rms_tmp(7,j) +                                 &
                                 (V_COMB(j,k) - V_READ(j,k))                    &
                                *(W_COMB(j,k) - W_READ(j,k))

                  rms_tmp(8,j) = rms_tmp(8,j) +                                 &
                                 (U_COMB(j,k) - U_READ(j,k))                    &
                                *(T_COMB(j,k) - T_READ(j,k))
                  rms_tmp(9,j) = rms_tmp(9,j) +                                 &
                                 (V_COMB(j,k) - V_READ(j,k))                    &
                                *(T_COMB(j,k) - T_READ(j,k))
                  rms_tmp(10,j) = rms_tmp(10,j) +                                 &
                                 (W_COMB(j,k) - W_READ(j,k))                    &
                                *(T_COMB(j,k) - T_READ(j,k))

              END DO

              U_tmp(1:4,j) = U_tmp(1:4,j)/Nz
              U_pr(1:4,j)  = ( U_pr(1:4,j) * (tt - 1) + U_tmp(1:4,j) )/tt

              rms_tmp(1:10,j) = rms_tmp(1:10,j)/Nz
              rms_pr(1:10,j) = ( rms_pr(1:10,j) * (tt - 1) + rms_tmp(1:10,j) )/tt

            END DO
            !OMP END PARALLEL

        END SUBROUTINE SEM_STAT


      END MODULE
