MODULE CONTAIN
    IMPLICIT NONE
    SAVE
!-----------------------------------------------------------------------
!*******************PARALLEL VARIBLES**********************
!-----------------------------------------------------------------------
    ! Parallel Assembly/Solver Variables
    ! INTEGER:: Nparts, Neqn_Global
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




!-----------------------------------------------------------------------
!*******************INITIAL PARAMETERS*************************
!-----------------------------------------------------------------------
    REAL(8)::n,ns,Fn
    INTEGER::n_cells_x,n_cells_y
    INTEGER::dt_to_save
    REAL(8)::t_final,dt
    REAL(8)::dx_0,dx_factor,dy_factor
    LOGICAL::include_source,close_inlet,close_outlet,include_gun_boundaries,use_homogenous_grid,restart_simulation
    CHARACTER(80)::dir_cur
    INTEGER::it_restart

!-----------------------------------------------------------------------
!*******************DYNAMIC ARRAYS*************************
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!*******************MISC. VARIBLES*************************
!-----------------------------------------------------------------------

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


! positive ion properties
    REAL(8), PARAMETER:: m_p= 1.67262178d-27
    REAL(8), PARAMETER:: mu_i = 1.45d3/p*1.d-4*phi_0*tau_0/(L_0**2)
    REAL(8), PARAMETER:: T_i = 350.d0
    REAL(8), PARAMETER:: v_i = Sqrt( (8.d0*k_b*T_i)/(Pi*m_i) )*Tau_0/L_0
    REAL(8), PARAMETER:: D_i = mu_i*(k_b*T_e/E_charge)/Phi_0


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
    ! WRITE(*,*) ' processor rank =', myrank
    ! WRITE(*,*) ' number of nodes in x =', nnodes_x
    ! WRITE(*,*) ' number of nodes in y=', nnodes_y
    ! WRITE(*,*) ' total number of nodes =', nnodes
    ! WRITE(*,*) ' number of bcs   =', nboundary
    ! WRITE(*,*) ' number of ics   =', ninitial
    ! WRITE(*,*) ' number of flux bcs   =', nflux_part+nflux_phi
    ! WRITE(*,*) ' number of equations =', neqn
    ! WRITE(*,*) ' number of dofs =', ndof


    WRITE(*,*) 'n=',n
    WRITE(*,*) 't_final=',t_final
    WRITE(*,*) 'include_source=',include_source
    WRITE(*,*) 'use_homogenous_grid=',use_homogenous_grid
    WRITE(*,*) 'dir_cur = ',dir_cur

END SUBROUTINE OUTPUT_TO_SCREEN







SUBROUTINE INPUT_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE
    CHARACTER(80):: Header_name(18), line
!-----------------------------------------------------------------------
!*******************INITIALIZE CHARACTER*******************
!-----------------------------------------------------------------------

    ! this isn't strictly necessary, I'm listing the headers directly
    Header_name = (/    '*N_INITIAL_____', &
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
        IF      (line == '*N_INITIAL_____') THEN
          ! initial density
          READ(100,*) n
        ELSE IF (line == '*NS_INITIAL____') THEN
          ! density of the source
          READ(100,*) Fn
        ELSE IF (line == '*FN____________') THEN
          ! number of particles represented by each superparticle
          READ(100,*) Fn
        ELSE IF (line == '*N_CELLS_X_____') THEN
          ! number of cells in the x direction (if using homogeneous grid)
          READ(100,*) n_cells_x
        ELSE IF (line == '*N_CELLLS_Y____') THEN
          ! number of cells in the y direction (if using homogeneous grid)
          READ(100,*) n_cells_y
        ELSE IF (line == '*T_FINAL_______') THEN
          ! end time of the simulation
          READ(100,*) t_final
        ELSE IF (line == '*DT____________') THEN
          ! timestep
          READ(100,*) dt
        ELSE IF (line == '*DT_TO_SAVE____') THEN
          ! how often to save simulation data (every 'dt-to-save'th timestep)
          READ(100,*) dt_to_save
        ELSE IF (line == '*INCLUDE_SOURCE') THEN
          ! whether to include source at inlet
          READ(100,*) include_source
        ELSE IF (line == '*CLOSE_INLET___') THEN
          ! whether to put boundary at inlet (after any source is turned off)
          READ(100,*) close_inlet
        ELSE IF (line == '*CLOSE_OUTLET__') THEN
          ! whether to put boundary at the outlet
          READ(100,*) close_outlet
        ELSE IF (line == '*INCLUDE_GUN_BC') THEN
          ! include the boundaries of the gun geometry
          READ(100,*) include_gun_boundaries
        ELSE IF (line == '*USE_HOMOG_GRID') THEN
          ! whether to use a homogeneous grid or the geometric one
          READ(100,*) use_homogenous_grid
        ELSE IF (line == '*DX_INIT_______') THEN
          ! initial cell-size (at center of inlet)
          READ(100,*) dx_0
        ELSE IF (line == '*DX_FACTOR_____') THEN
          ! relative cell-size at the outlet in the x-direction
          READ(100,*) dx_factor
        ELSE IF (line == '*DY_FACTOR_____') THEN
          ! relative cell-size at edges in the y-direction
          READ(100,*) dy_factor
        ELSE IF (line == '*RESTART_SIM___') THEN
          ! whether or not to restart a partially-completed simulation
          READ(100,*) restart_simulation
        ELSE IF (line == '*DIR_RESTART___') THEN
          ! directory of current data
          READ(100,*) dir_cur
        ELSE IF (line == '*RESTART_INDEX_') THEN
          ! index to restart from
          READ(100,*) it_restart
          EXIT
       ELSE
          ! Error Analyzing Input
          WRITE(*,*) 'Error Analyzing Input_Parameters by Processor' &
               , Myrank
          WRITE(*,*) line
          STOP
       END IF
    END DO
    CLOSE(100)
END SUBROUTINE INPUT_PARAMETERS_READIN



SUBROUTINE INITIALIZE
    USE CONTAIN
    IMPLICIT NONE
    ! This suboutine initializes various matrices/parameters in the simulation
    

END SUBROUTINE INITIALIZE




