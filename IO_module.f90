!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : IO_module.f90                                                    !
!                                                                              !
!   PURPOSE : Write each variables in the RESULT folder.                       !
!                                                                              !
!                                                             2017.10.02 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

          MODULE IO_module

            CHARACTER(LEN=65) :: file_name, dir_name, path_name

          !--------------------------------------------------------------------!
          !                   Interfaces of IO Subroutines                     !
          !--------------------------------------------------------------------!

            !--------Folder Initializing Subroutine
            INTERFACE FOLDER_SETUP
              MODULE PROCEDURE FOLDER_SETUP
            END INTERFACE FOLDER_SETUP

          CONTAINS
            SUBROUTINE FOLDER_SETUP
            IMPLICIT NONE

            !------------------------------------------------------------------!
            !                  Make & Initialize Result folder                 !
            !------------------------------------------------------------------!
            dir_name  = 'RESULT'
            CALL SYSTEM('mkdir '//TRIM(dir_name))
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/*.plt')

            !------------------------------------------------------------------!
            !                  Reading Data Folder Designation                 !
            !------------------------------------------------------------------!
            dir_name = 'DATA'
            file_name = 'Slice_SBL.plt'
            path_name = TRIM(dir_name)//'/'//TRIM(file_name)
            
            END SUBROUTINE FOLDER_SETUP

          END MODULE IO_module
