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
            USE SEM_module,                                                     &
              ONLY : N, Ny, Nz, time, file_name, dir_name, path_name

            USE SEM_module,                                                     &
              ONLY : Y, Z, U, V, W, T, RS, THS,                                 &
                     U_INLET, V_INLET, W_INLET, T_INLET, SEM_EDDY,    &
                     U_COMB, V_COMB, W_COMB, T_COMB, U_c, U_pr, rms_pr

            USE flow_module

              IMPLICIT NONE
              INTEGER :: it,j,k
              REAL(KIND=8) :: time_sta, time_end, U_mean(4,1:Nz), rms_mean(10,1:Nz)

              U_mean(1:4,1:Nz)   = 0.0
              rms_mean(1:10,1:Nz) = 0.0

              ! WRITE(*,*) '----------------------------------------------------'
              ! WRITE(*,*) '              WRITING PROCESS STARTED               '
              CALL CPU_TIME(time_sta)
              dir_name = 'RESULT'
!               !----------------------------------------------------------------!
!               !                       Outputs for U Slice                      !
!               !----------------------------------------------------------------!
!               file_name = '/U_ins.plt'
!               path_name = TRIM(dir_name)//TRIM(file_name)
!
!               OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
!               WRITE(100,*) 'VARIABLES = Y,Z,U_ins,V_ins,W_ins,T_ins'
!               WRITE(100,"(2(A,I3,2X))")' ZONE  I = ',Nz,' J = ',Ny
!               WRITE(100,*) 'SOLUTIONTIME =',time
!
!               DO j = 1,Ny
!                 DO k = 1,Nz
!
!                     WRITE(100,"(6F15.9)") Y(j),Z(k),                            &
!                                           U_COMB(j,k),V_COMB(j,k),W_COMB(j,k),  &
!                                           T_COMB(j,k)
!
!                 END DO
!               END DO
!               WRITE(100,*)
!               CLOSE(100)
! !
              !----------------------------------------------------------------!
              !                   Outputs for U mean profiles                  !
              !----------------------------------------------------------------!
              file_name = '/Mean_profiles.plt'
              path_name = TRIM(dir_name)//TRIM(file_name)

              IF ( INT(time/dt) == Nt) THEN
                OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
                WRITE(100,*)'VARIABLES = Z,U_ins,U_exac,V_ins,V_exac,W_ins,W_exac,T_ins,T_exac'

                DO k = 1,Nz
                  DO j = 1,Ny
                    U_mean(1,k) = U_mean(1,k) + U(j,k)
                    U_mean(2,k) = U_mean(2,k) + V(j,k)
                    U_mean(3,k) = U_mean(3,k) + W(j,k)
                    U_mean(4,k) = U_mean(4,k) + T(j,k)
                  END DO
                    U_mean(1:4,k) = U_mean(1:4,k)/Ny

                    WRITE(100,"(9F15.9)") Z(k), U_pr(1,k), U_mean(1,k),         &
                                                U_pr(2,k), U_mean(2,k),         &
                                                U_pr(3,k), U_mean(3,k),         &
                                                U_pr(4,k), U_mean(4,k)
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
                WRITE(100,*)'VARIABLES = Z,uu,uu_exac,vv,vv_exac,ww,ww_exac,tt,tt_exac,uv,uv_exac,uw,uw_exac,vw,vw_exac,upt,upt_exac,vpt,vpt_exac,wpt,wpt_exac'

                DO k = 1,Nz
                  DO j = 1,Ny
                    rms_mean(1,k) = rms_mean(1,k) + RS(1,j,k)
                    rms_mean(2,k) = rms_mean(2,k) + RS(2,j,k)
                    rms_mean(3,k) = rms_mean(3,k) + RS(3,j,k)
                    rms_mean(4,k) = rms_mean(4,k) + THS(1,j,k)

                    rms_mean(5,k) = rms_mean(5,k) + RS(4,j,k)
                    rms_mean(6,k) = rms_mean(6,k) + RS(5,j,k)
                    rms_mean(7,k) = rms_mean(7,k) + RS(6,j,k)

                    rms_mean(8,k) = rms_mean(8,k) + THS(2,j,k)
                    rms_mean(9,k) = rms_mean(9,k) + THS(3,j,k)
                    rms_mean(10,k) = rms_mean(10,k) + THS(4,j,k)
                  END DO
                    rms_mean(1:10,k) = rms_mean(1:10,k)/Ny
                    WRITE(100,"(21F15.9)") Z(k), rms_pr(1,k), rms_mean(1,k),     &
                                                 rms_pr(2,k), rms_mean(2,k),     &
                                                 rms_pr(3,k), rms_mean(3,k),     &
                                                 rms_pr(4,k), rms_mean(4,k),     &
                                                 rms_pr(5,k), rms_mean(5,k),     &
                                                 rms_pr(6,k), rms_mean(6,k),     &
                                                 rms_pr(7,k), rms_mean(7,k),     &
                                                 rms_pr(8,k), rms_mean(8,k),     &
                                                 rms_pr(9,k), rms_mean(9,k),     &
                                                 rms_pr(10,k), rms_mean(10,k)
                END DO
                WRITE(100,*)
                CLOSE(100)
              END IF

              ! !----------------------------------------------------------------!
              ! !                   Outputs for Eddy posoitions                  !
              ! !----------------------------------------------------------------!
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

              CALL CPU_TIME(time_end)

              ! WRITE(*,*) '           WRITING PROCESS IS COMPLETED            '
              ! WRITE(*,*) '  Total Writing time : ',time_end - time_sta,' s'
              ! WRITE(*,*) '----------------------------------------------------'
              ! WRITE(*,*) ''

              IF ( INT(time/dt) == Nt) THEN
                DEALLOCATE(Y,Z,U,V,W,RS,U_INLET,V_INLET,W_INLET,SEM_EDDY)
                DEALLOCATE(U_COMB,V_COMB,W_COMB,U_pr,rms_pr,U_c)
              END IF
          END SUBROUTINE OUTPUT
