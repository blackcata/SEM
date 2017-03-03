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
          REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: U,V,W
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: RS
          TYPE(EDDY_CHAR),DIMENSION(:),ALLOCATABLE :: SEM_EDDY

        CONTAINS
          FUNCTION INTENSITY_det(x_int)
            REAL(KIND=8) :: INTENSITY_det
            REAL(KIND=8),INTENT(IN) :: x_int

            IF ( x_int > 0 ) THEN
              INTENSITY_det = 1
            ELSE
              INTENSITY_det = -1
            END IF
          END FUNCTION INTENSITY_det

        END MODULE
