!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : SEM_write.f90                                                    !
!                                                                              !
!   PURPOSE : Write each variables in the RESULT folder.                       !
!                                                                              !
!                                                             2017.03.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!
          SUBROUTINE OUTPUT
            USE flow_module,                                                    &
              ONLY :  Nt, dt

            USE SEM_module,                                                     &
              ONLY : N, Ny, Nz, time

            USE SEM_module,                                                     &
              ONLY : Y, Z, U_READ, V_READ, W_READ, T_READ, RS, THS,             &
                     U_INLET, V_INLET, W_INLET, T_INLET, SEM_EDDY,    &
                     U_COMB, V_COMB, W_COMB, T_COMB, U_c, U_pr, rms_pr

            USE IO_module

              IMPLICIT NONE
              INTEGER :: it,j,k
              REAL(KIND=8) :: U_mean(4,1:Ny), rms_mean(10,1:Ny), UT, TT, NUT

              U_mean(1:4,1:Ny)   = 0.0
              rms_mean(1:10,1:Ny) = 0.0

              dir_name = 'RESULT'

              UT  = 5.5e-2
              TT  = 6.58e-2
              NUT = 6.06e-2

              !----------------------------------------------------------------!
              !                       Outputs for U Slice                      !
              !----------------------------------------------------------------!
            !  file_name = '/U_ins.plt'
            !  path_name = TRIM(dir_name)//TRIM(file_name)
            !
            !  OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
            !  WRITE(100,*) 'VARIABLES = Z,Y,U_ins,V_ins,W_ins,T_ins'
            !  WRITE(100,"(2(A,I3,2X))")' ZONE  I = ',Nz,' J = ',Ny
            !  WRITE(100,*) 'SOLUTIONTIME =',time
            !
            !  DO j = 1,Ny
            !    DO k = 1,Nz
            !
            !        WRITE(100,"(5F15.9)") Z(k),Y(j),                            &
            !                              U_COMB(j,k),V_COMB(j,k),W_COMB(j,k),  &
            !                              T_COMB(j,k)
            !
            !    END DO
            !  END DO
            !  WRITE(100,*)
            !  CLOSE(100)

              !----------------------------------------------------------------!
              !                   Outputs for U mean profiles                  !
              !----------------------------------------------------------------!
              file_name = '/Mean_profiles.plt'
              path_name = TRIM(dir_name)//TRIM(file_name)

              IF ( INT(time/dt) == Nt) THEN
                OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
                WRITE(100,*)'VARIABLES = Y,U_ins,U_exac,V_ins,V_exac,W_ins,W_exac,T_ins,T_exac'

                DO j = 1,Ny
                  DO k = 1,Nz
                    U_mean(1,j) = U_mean(1,j) + U_READ(j,k)
                    U_mean(2,j) = U_mean(2,j) + V_READ(j,k)
                    U_mean(3,j) = U_mean(3,j) + W_READ(j,k)
                    U_mean(4,j) = U_mean(4,j) + T_READ(j,k)
                  END DO
                    U_mean(1:4,j) = U_mean(1:4,j)/Nz
!                    WRITE(100,"(9F15.9)") Y(j), U_pr(1,j), U_mean(1,j),         &
!                                                U_pr(2,j), U_mean(2,j),         &
!                                                U_pr(3,j), U_mean(3,j),         &
!                                                U_pr(4,j), U_mean(4,j)

                    WRITE(100,"(9F15.9)") Y(j)/NUT, U_pr(1,j)/UT, U_mean(1,j)/UT,         &
                                                U_pr(2,j)/UT, U_mean(2,j)/UT,         &
                                                U_pr(3,j)/UT, U_mean(3,j)/UT,         &
                                                U_pr(4,j)/TT, U_mean(4,j)/TT
                END DO
                WRITE(100,*)
                CLOSE(100)
              END IF

              !----------------------------------------------------------------!
              !              Outputs for Reynolds stress profiles              !
              !----------------------------------------------------------------!
              file_name = '/RMS_profiles.plt'
              path_name = TRIM(dir_name)//TRIM(file_name)

              IF ( INT(time/dt) == Nt) THEN
                OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
                WRITE(100,*)'VARIABLES = Y,uu,uu_exac,vv,vv_exac,ww,ww_exac,tt,tt_exac,uv,uv_exac,uw,uw_exac,vw,vw_exac,ut,ut_exac,vt,vt_exac,wt,wt_exac'

                DO j = 1,Ny
                  DO k = 1,Nz
                    rms_mean(1,j) = rms_mean(1,j) + RS(1,j,k)
                    rms_mean(2,j) = rms_mean(2,j) + RS(2,j,k)
                    rms_mean(3,j) = rms_mean(3,j) + RS(3,j,k)
                    rms_mean(4,j) = rms_mean(4,j) + THS(1,j,k)

                    rms_mean(5,j) = rms_mean(5,j) + RS(4,j,k)
                    rms_mean(6,j) = rms_mean(6,j) + RS(5,j,k)
                    rms_mean(7,j) = rms_mean(7,j) + RS(6,j,k)

                    rms_mean(8,j) = rms_mean(8,j) + THS(2,j,k)
                    rms_mean(9,j) = rms_mean(9,j) + THS(3,j,k)
                    rms_mean(10,j) = rms_mean(10,j) + THS(4,j,k)
                  END DO
                    rms_mean(1:10,j) = rms_mean(1:10,j)/Nz
                  !  WRITE(100,"(21F15.9)") Y(j), rms_pr(1,j), rms_mean(1,j),    &
                  !                               rms_pr(2,j), rms_mean(2,j),    &
                  !                               rms_pr(3,j), rms_mean(3,j),    &
                  !                               rms_pr(4,j), rms_mean(4,j),    &
                  !                               rms_pr(5,j), rms_mean(5,j),    &
                  !                               rms_pr(6,j), rms_mean(6,j),    &
                  !                               rms_pr(7,j), rms_mean(7,j),    &
                  !                               rms_pr(8,j), rms_mean(8,j),    &
                  !                               rms_pr(9,j), rms_mean(9,j),    &
                  !                               rms_pr(10,j), rms_mean(10,j)
                    WRITE(100,"(21F15.9)") Y(j)/NUT, rms_pr(1,j)/UT**2, rms_mean(1,j)/UT**2,    &
                                                 rms_pr(2,j)/UT**2, rms_mean(2,j)/UT**2,    &
                                                 rms_pr(3,j)/UT**2, rms_mean(3,j)/UT**2,    &
                                                 rms_pr(4,j)/TT**2, rms_mean(4,j)/TT**2,    &
                                                 rms_pr(5,j)/UT**2, rms_mean(5,j)/UT**2,    &
                                                 rms_pr(6,j)/UT**2, rms_mean(6,j)/UT**2,    &
                                                 rms_pr(7,j)/UT**2, rms_mean(7,j)/UT**2,    &
                                                 rms_pr(8,j)/(UT*TT), rms_mean(8,j)/(UT*TT),    &
                                                 rms_pr(9,j)/(UT*TT), rms_mean(9,j)/(UT*TT),    &
                                                 rms_pr(10,j)/(UT*TT), rms_mean(10,j)/(UT*TT)
                END DO
                WRITE(100,*)
                CLOSE(100)
              END IF

              !----------------------------------------------------------------!
              !                   Outputs for Eddy posoitions                  !
              !----------------------------------------------------------------!
              ! file_name = '/EDDY_POS.plt'
              ! path_name = TRIM(dir_name)//TRIM(file_name)
              !
              ! OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
              ! WRITE(100,*) 'VARIABLES = X,Y,Z'
              ! WRITE(100,*) 'ZONE'
              ! WRITE(100,*) 'SOLUTIONTIME =',time
              !
              ! DO it = 1,N
              !   WRITE(100,*) SEM_EDDY(it)%X_pos*10,SEM_EDDY(it)%Y_pos,          &
              !                SEM_EDDY(it)%Z_pos
              ! END DO
              ! WRITE(100,*)
              ! CLOSE(100)

              IF ( INT(time/dt) == Nt) THEN
                DEALLOCATE(Y,Z,U_READ,V_READ,W_READ,RS,THS)
                DEALLOCATE(U_INLET,V_INLET,W_INLET,SEM_EDDY)
                DEALLOCATE(U_COMB,V_COMB,W_COMB,U_pr,rms_pr,U_c)
              END IF
          END SUBROUTINE OUTPUT
