SUBROUTINE INPUT_READIN
    USE CONTAIN
    IMPLICIT NONE
!-----------------------------------------------------------------------
!*******************READ PARAMETER FILE********************
!-----------------------------------------------------------------------
    CALL INPUT_PARAMETERS_READIN
!-----------------------------------------------------------------------
!*******************READ MESH INPUT FILE*******************
!-----------------------------------------------------------------------
    CALL MESH_READIN
END SUBROUTINE INPUT_READIN


SUBROUTINE INPUT_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE
    CHARACTER(80):: Header_name(9), line
!-----------------------------------------------------------------------
!*******************INITIALIZE CHARACTER*******************
!-----------------------------------------------------------------------
    Header_name = (/ '*STENCIL_WIDTH', '*RF', '*NPARTS  ', &
         '*PRINT', '*TIME_ORDER', '*EPSILON_ABS', &
         '*INIT_TIME   ', '*FINALTIME   ', '*DELTA_TIME '/)
!-----------------------------------------------------------------------
!*******************OPEN FILE******************************
!-----------------------------------------------------------------------
    OPEN(UNIT=100,FILE='Input/Input_Parameters')
!-----------------------------------------------------------------------
!*******************READ CONTENTS**************************
!-----------------------------------------------------------------------
    DO
       READ(100,*) line
       IF (line == Header_name(1)) THEN
          ! Read Stencil Width for Scheme
          READ(100,*) Stencil_Width
       ELSE IF (line == Header_name(2)) THEN
          ! RF Check
          READ(100,*) RF
       ELSE IF (line == Header_name(3)) THEN
          ! Number of Subdomains
          READ(100,*) Nparts
       ELSE IF (line == Header_name(4)) THEN
          ! Print Output Every (Print) Timesteps
          READ(100,*) Print
       ELSE IF (line == Header_name(5)) THEN
          ! Runge-Kutta Temporal Order
          READ(100,*) Runge_Kutta_Order
       ELSE IF (line == Header_name(6)) THEN
          ! Dormand-Prince Order Parameter
          READ(100,*) Time_Error_Threshold
       ELSE IF (line == Header_name(7)) THEN
          ! Initial Starting Time of Simulation
          READ(100,*) Initial_Time
       ELSE IF (line == Header_name(8)) THEN
          ! Maximum Number of Newton Subiterations
          READ(100,*) Final_Time
       ELSE IF (line == Header_name(9)) THEN
          ! Residual Tolerance Threshold
          READ(100,*) Delta_time
          EXIT
       ELSE
          ! Error Analyzing Input
          WRITE(*,*) 'Error Analyzing Input_Parameters by Processor' &
               , Myrank
          WRITE(*,*) line
          CALL MPI_FINALIZE(ierr)
          STOP
       END IF
    END DO
    CLOSE(100)
END SUBROUTINE INPUT_PARAMETERS_READIN


