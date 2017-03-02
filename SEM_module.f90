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
          INTEGER :: N, Ny, Nz, Nt
          REAL(KIND=8) :: dt, SIGMA, V_b
          CHARACTER(LEN=65) :: file_name, dir_name, path_name

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Y,Z
          REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: U,V,W
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: RS

        END MODULE
