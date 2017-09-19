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

            USE SEM_module,                                                     &
              ONLY : N, Ny, Nz, Nt, dt, SIGMA, V_b, file_name, dir_name, OUT_NUM

            USE SEM_module,                                                     &
              ONLY : Y, Z, U, V, W, T, RS, THS, U_INLET, V_INLET, W_INLET,      &
                      SEM_EDDY, U_COMB, V_COMB, W_COMB, U_pr, rms_pr, U_c,      &
                      T_INLET, T_COMB

            IMPLICIT NONE
            INTEGER :: i,j,k

            !------------------------------------------------------------------!
            !                  Make & Initialize Result folder                 !
            !------------------------------------------------------------------!
            dir_name  = 'RESULT'
            CALL SYSTEM('mkdir '//TRIM(dir_name))
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/*.plt')

            !------------------------------------------------------------------!
            !                         Constants for SEM                        !
            !------------------------------------------------------------------!
            N  = 1000
            Ny = 96
            Nz = 256

            Nt    = 40000
            dt    = 5e-3
            SIGMA = 0.20

            OUT_NUM = 100

            !------------------------------------------------------------------!
            !                         Allocate variables                       !
            !------------------------------------------------------------------!
            ALLOCATE( Y(1:Ny),Z(1:Nz) )
            ALLOCATE( U(1:Ny,1:Nz), V(1:Ny,1:Nz), W(1:Ny,1:Nz), T(1:Ny,1:Nz) )
            ALLOCATE( U_INLET(1:Ny,1:Nz),V_INLET(1:Ny,1:Nz),W_INLET(1:Ny,1:Nz) )
            ALLOCATE( U_COMB(1:Ny,1:Nz),V_COMB(1:Ny,1:Nz),W_COMB(1:Ny,1:Nz) )
            ALLOCATE( T_INLET(1:Ny,1:Nz), T_COMB(1:Ny,1:Nz)  )
            ALLOCATE( RS(6,1:Ny,1:Nz), THS(4,1:Ny,1:Nz), U_c(1:Ny,1:Nz) )
            ALLOCATE( SEM_EDDY(1:N), U_pr(4,1:Ny), rms_pr(8,1:Ny) )

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
            rms_pr(1:8,1:Ny) = 0.0

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
