MODULE CONTAIN
    IMPLICIT NONE
    SAVE
!-----------------------------------------------------------------------
!*******************PARALLEL VARIBLES**********************
!-----------------------------------------------------------------------
    ! Parallel Assembly/Solver Variables
    INTEGER:: Nparts, Neqn_Global
!-----------------------------------------------------------------------
!*******************DYNAMIC ARRAYS*************************
!-----------------------------------------------------------------------
    ! Relevant Allocatable Arrays
    REAL(8), ALLOCATABLE, DIMENSION(:):: x, y
    REAL(8), ALLOCATABLE, DIMENSION(:):: U, Material_Label, &
         Valboundary
    REAL(8), ALLOCATABLE, DIMENSION(:):: Valinitial
    INTEGER, ALLOCATABLE, DIMENSION(:):: nodeboundary_x, &
         nodeboundary_y, Dofboundary, Flux_Label
    INTEGER, ALLOCATABLE, DIMENSION(:):: Nodeinitial, Dofinitial
    INTEGER, ALLOCATABLE, DIMENSION(:, :):: Node_Type_Part_X, Node_Type_Part_Y
    INTEGER, ALLOCATABLE, DIMENSION(:, :):: Node_Type_Phi_X, Node_Type_Phi_Y
    INTEGER, ALLOCATABLE, DIMENSION(:, :):: Globaldof, Localdof
    REAL(8), ALLOCATABLE, DIMENSION(:, :):: n_i_plus, n_i_minus, n_i_minus_orig
    REAL(8), ALLOCATABLE, DIMENSION(:, :):: n_e_plus, n_e_minus, n_e_minus_orig
    REAL(8), ALLOCATABLE, DIMENSION(:, :):: Error_Offset_n_i, Error_Offset_n_e
    REAL(8), ALLOCATABLE, DIMENSION(:, :):: phi_plus
    REAL(8), ALLOCATABLE, DIMENSION(:, :):: Sigma_Charge
    REAL(8), ALLOCATABLE, DIMENSION(:, :, :):: K_i, K_e
!-----------------------------------------------------------------------
!*******************MISC. VARIBLES*************************
!-----------------------------------------------------------------------
    LOGICAL:: Petsc_Logical, Assembly_Check
    INTEGER:: Myrank, Mysize, Ierr, Iter, Fiter, Rf, Nonzeros_Per_Row
    REAL(8):: Clock_Start, Clock_Finish
    INTEGER:: Nnodes, Nnodes_x, Nnodes_y, Ndimensions, Nboundary, Ninitial
    INTEGER:: NFlux_Part, NFlux_Phi, Size
    INTEGER:: Ndof, Neqn, Print, Jacobian_Index, Max_Dof_of_Node, Runge_Kutta_Order
    REAL(8):: Initial_Time, Final_Time, Error_Offset, Scaling_Factor, Time_Error_Threshold
    REAL(8):: Delta_Time
    INTEGER:: Stencil_Width, Jacobian_Size
END MODULE CONTAIN


MODULE PROPERTIES
    IMPLICIT NONE
    SAVE

! non-dimensional parameters
    REAL(8), PARAMETER:: L_0 = 1.d-4
    REAL(8), PARAMETER:: Phi_0 = 1000.d0
    REAL(8), PARAMETER:: Tau_0 = 1.d-6
    REAL(8), PARAMETER:: n_0 = 1.d22

! fundamental constants
    REAL(8), PARAMETER:: Epsilon = 8.85418782d-12
    REAL(8), PARAMETER:: E_charge = 1.60217646d-19
    REAL(8), PARAMETER:: Pi = 4.d0*atan(1.d0)
    REAL(8), PARAMETER:: k_b = 1.3806503d-23

! case properties
    REAL(8), PARAMETER:: eps_d = 10.d0
    REAL(8), PARAMETER:: eps_p = 1.d0
    REAL(8), PARAMETER:: p = 760.d0

! electron properties
    REAL(8), PARAMETER:: m_e = 9.10938188d-31
    REAL(8), PARAMETER:: mu_e = 4.4d5/p*1.d-4*phi_0*tau_0/(L_0**2)
    REAL(8), PARAMETER:: T_e = 10000.d0
    REAL(8), PARAMETER:: v_e = Sqrt( (8.d0*k_b*T_e)/(Pi*m_e) )*Tau_0/L_0
    REAL(8), PARAMETER:: D_e = mu_e*(k_b*T_e/E_charge)/Phi_0

! positive ion properties
    REAL(8), PARAMETER:: m_i = 1.67262178d-27 * 28.d0
    REAL(8), PARAMETER:: mu_i = 1.45d3/p*1.d-4*phi_0*tau_0/(L_0**2)
    REAL(8), PARAMETER:: T_i = 350.d0
    REAL(8), PARAMETER:: v_i = Sqrt( (8.d0*k_b*T_i)/(Pi*m_i) )*Tau_0/L_0
    REAL(8), PARAMETER:: D_i = mu_i*(k_b*T_e/E_charge)/Phi_0

! reactions
    REAL(8), PARAMETER:: A = 12.d2*L_0
    REAL(8), PARAMETER:: B = 342.d2*L_0/Phi_0
    REAL(8), PARAMETER:: Beta = 1.d-13*Tau_0*n_0
    REAL(8), PARAMETER:: Gamma = 5.0d-2

END MODULE PROPERTIES


PROGRAM MAIN
    USE CONTAIN
    IMPLICIT NONE
    INTEGER:: Timestep
    REAL(8):: Tot_Time


!-----------------------------------------------------------------------
!*******************READ INPUT FILE************************
!-----------------------------------------------------------------------
    ! CALL INPUT_READIN
    CALL INPUT_PARAMETERS_READIN
!-----------------------------------------------------------------------
!*******************SET INITIAL CONDITIONS*****************
!-----------------------------------------------------------------------
    tot_time = initial_time
    Timestep = 0
    assembly_check = .true.
    ! CALL INITIALIZE
    ! ! Write the initial conditions
    ! CALL TECPLOT(Timestep, Initial_Time)
    IF (Myrank == 0) THEN
       CALL OUTPUT_TO_SCREEN
    END IF
! !-----------------------------------------------------------------------
! !*******************MAIN LOOP******************************
! !-----------------------------------------------------------------------
!     DO
!        ! Increment the total time
!        tot_time = tot_time + delta_time
!        timestep = timestep + 1
!        ! Solve/Form Matrices Using PETSc
!        CALL PETSC
! !-----------------------------------------------------------------------
! !*******************EXPLICIT EVALUATION********************
! !----------------------------------------------------------------------- 
!        CALL EXPLICIT_SOLVE
!        n_e_minus = n_e_plus
!        n_e_minus_orig = n_e_plus
!        n_i_minus = n_i_plus
!        n_i_minus_orig = n_i_plus
! !-----------------------------------------------------------------------
! !***********************  TECPLOT  ******************************
! !-----------------------------------------------------------------------
!        if (float(timestep/print) .eq. float(timestep)/float(print))then
!           if(myrank.eq.0)then
!              WRITE(*,*) 'Step #', Timestep, 'DT: ', delta_time, &
!                   'Present Time: ', tot_time, 'Final Time: ', final_time
!           end if
!           call TECPLOT (timestep,tot_time)
!        end if
!        IF(tot_time >= final_time) EXIT
!     END DO
!-----------------------------------------------------------------------
!*******************END MPI*********************************
!-----------------------------------------------------------------------
    WRITE(*,*) '------ done ------'

    STOP
END PROGRAM MAIN


SUBROUTINE OUTPUT_TO_SCREEN
    USE CONTAIN
    IMPLICIT NONE
    ! This Suboutine Prints Out Relevant Domain Info
    WRITE(*,*)
    WRITE(*,*) ' ********   ANALYSIS INFORMATION ********'
    WRITE(*,*) ' processor rank =', myrank
    WRITE(*,*) ' number of nodes in x =', nnodes_x
    WRITE(*,*) ' number of nodes in y=', nnodes_y
    WRITE(*,*) ' total number of nodes =', nnodes
    WRITE(*,*) ' number of bcs   =', nboundary
    WRITE(*,*) ' number of ics   =', ninitial
    WRITE(*,*) ' number of flux bcs   =', nflux_part+nflux_phi
    WRITE(*,*) ' number of equations =', neqn
    WRITE(*,*) ' number of dofs =', ndof
    WRITE(*,*) 'Initial_time = ',Initial_time
    WRITE(*,*) 'Final_time = ',Final_time
    WRITE(*,*)
END SUBROUTINE OUTPUT_TO_SCREEN







SUBROUTINE INPUT_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE
    CHARACTER(80):: Header_name(19), line
!-----------------------------------------------------------------------
!*******************INITIALIZE CHARACTER*******************
!-----------------------------------------------------------------------
    ! Header_name = (/ '*STENCIL_WIDTH', '*RF', '*NPARTS  ', &
    !      '*PRINT', '*TIME_ORDER', '*EPSILON_ABS', &
    !      '*INIT_TIME   ', '*FINALTIME   ', '*DELTA_TIME '/)
    ! Header_name = (/    '*STENCIL_WIDTH_', &
    !                     '*RF____________', &
    !                     '*NPARTS________', &
    !                     '*PRINT_________', &
    !                     '*TIME_ORDER____', &
    !                     '*EPSILON_ABS___', &
    !                     '*INIT_TIME_____', &
    !                     '*FINALTIME_____', &
    !                     '*DELTA_TIME____'/)

    Header_name = (/    '*DELTA_TIME____', &
                        '*N_INITIAL_____', &
                        '*FN____________', &
                        '*N_CELLS_X_____', &
                        '*N_CELLLS_Y____', &
                        '*T_FINAL_______', &
                        '*DT____________', &
                        '*DT_TO_SAVE____', &
                        '*INCLUDE_SOURCE', &
                        '*CLOSE_INLET___', &
                        '*CLOSE_OUTLET__', &
                        '*INCLUDE_GUN_BC', &
                        '*USE_HOMOG_GRID', &
                        '*DX_INIT_______', &
                        '*DX_FACTOR_____', &
                        '*DY_FACTOR_____', &
                        '*RESTART_SIM___', &
                        '*DIR_RESTART___', &
                        '*RESTART_INDEX_'/)

    ! Header_name = (/    '*STENCIL_WIDTH', &
    !                     '*RF', &
    !                     '*NPARTS', &
    !                     '*PRINT', &
    !                     '*TIME_ORDER', &
    !                     '*EPSILON_ABS', &
    !                     '*INIT_TIME', &
    !                     '*FINALTIME', &
    !                     '*DELTA_TIME'/)

!-----------------------------------------------------------------------
!*******************OPEN FILE******************************
!-----------------------------------------------------------------------
    OPEN(UNIT=100,FILE='Input/Input_Parameters')
!-----------------------------------------------------------------------
!*******************READ CONTENTS**************************
!-----------------------------------------------------------------------
    WRITE(*,*) 'got here'
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
          ! CALL MPI_FINALIZE(ierr)
          STOP
       END IF
    END DO
    CLOSE(100)
END SUBROUTINE INPUT_PARAMETERS_READIN
