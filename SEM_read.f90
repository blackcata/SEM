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

          USE SEM_module,                                                     &
            ONLY : N, Ny, Nz, SIGMA, V_b, file_name, dir_name, path_name

          USE SEM_module,                                                     &
            ONLY : Y, Z, U, V, W, RS

          IMPLICIT NONE

          INTEGER :: j,k
          REAL(KIND=8) :: time_sta, time_end, tmp_y, tmp_z, S
          CHARACTER(20) :: header

          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) '              READING PROCESS STARTED               '
          CALL CPU_TIME(time_sta)

          dir_name = 'DATA'
          file_name = 'SLICE.plt'
          path_name = TRIM(dir_name)//'/'//TRIM(file_name)

          OPEN(100,FILE=path_name,FORM='FORMATTED',STATUS='OLD')
          READ(100,*) header
          READ(100,*) header

          !--------------------------------------------------------------------!
          !                  Main loop of reading SLICE datas                  !
          !--------------------------------------------------------------------!
          DO k = 1,Nz
            DO j = 1,Ny

                READ(100,*) tmp_z, tmp_y, U(j,k), V(j,k), W(j,k),               &
                            RS(1,j,k), RS(2,j,k), RS(3,j,k), RS(4,j,k),         &
                            RS(5,j,k), RS(6,j,k)

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
