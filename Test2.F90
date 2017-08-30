MODULE CONTAIN
    IMPLICIT NONE
    SAVE
!-----------------------------------
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

  PROGRAM Test2
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
  END PROGRAM Test2


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


SUBROUTINE MESH_READIN
    USE CONTAIN
    IMPLICIT NONE
    CHARACTER(80):: Header_name(7), line
    INTEGER:: I, J, Status
    INTEGER:: index_x, index_y, type_x, type_y
!-----------------------------------------------------------------------
!*******************INITIALIZE CHARACTER*******************
!-----------------------------------------------------------------------
    Header_name = (/ '*DIMENSIONS ', '*NODE ', '*BOUNDARY', &
         '*INITIAL ', '*FLUX_PART', '*FLUX_PHI', '*STOP '/)
!-----------------------------------------------------------------------
!*******************OPEN FILE******************************
!-----------------------------------------------------------------------
    OPEN(UNIT=200,FILE='Input/Domain_Input')
!-----------------------------------------------------------------------
!*******************READ CONTENTS**************************
!-----------------------------------------------------------------------
    DO
       READ(200,*) line
       IF (line == Header_name(1)) THEN
          ! Number of Dimensions
          READ(200,*) Ndimensions
       ELSE IF (line == Header_name(2)) THEN
          ! Readin Global Node Data
          READ(200, *) Nnodes_x
          READ(200, *) Nnodes_y
          Nnodes = Nnodes_x * Nnodes_y
          ALLOCATE(x(Nnodes_x), y(Nnodes_y), STAT=STATUS)
          IF (STATUS == 0) THEN
             ! Allocation was Successful
             DO 	I=1, Nnodes_x
                READ(200,*) index_x, x(i)
             END DO
             DO  I=1, Nnodes_y
                READ(200,*) index_y, y(i)
             END DO
          ELSE
             ! Allocation Failed
             WRITE(*,*) 'Allocation Error With: X and Y'
          END IF
       ELSE IF (line == Header_name(3)) THEN
          ! Readin Boundary Condition Data
          READ(200, *) Nboundary
          ALLOCATE( Nodeboundary_x(nboundary), Nodeboundary_y(nboundary), &
               Dofboundary(nboundary),valboundary(nboundary), STAT=STATUS)
          IF (STATUS == 0) THEN
             ! Allocation was Successful
             DO I=1, Nboundary
                READ(200,*) nodeboundary_x(I), nodeboundary_y(I), &
                     dofboundary(I), valboundary(I)
             END DO
          ELSE
             ! Allocation Failed
             WRITE(*,*) 'Allocation Error With: Node_Boundary'
          END IF
       ELSE IF (line == Header_name(4)) THEN
          ! Readin Initial Condition Data
          READ(200, *) Ninitial
          ALLOCATE( nodeinitial(ninitial), dofinitial(ninitial), &
               valinitial(ninitial), STAT=STATUS)
          IF (STATUS == 0) THEN
             ! Allocation was Successful
             DO I=1, Ninitial
                READ(200,*) nodeinitial(I), dofinitial(I), &
                     valinitial(I)
             END DO
          ELSE
             ! Allocation Failed
             WRITE(*,*) 'Allocation Error With: NodeInitial'
          END IF
       ELSE IF (line == Header_name(5)) THEN
          ! Readin Flux Boundary Conditions for Particle Species
          READ(200, *) NFlux_Part
          ALLOCATE(Node_Type_Part_X(Nnodes_x,Nnodes_y), &
               Node_Type_Part_Y(Nnodes_x,Nnodes_y), STAT=STATUS)
          Node_Type_Part_X = 0
          Node_Type_Part_Y = 0
          IF (STATUS==0) THEN
             DO I=1, NFlux_Part
                READ(200,*) index_x, index_y, Type_X, Type_Y
                Node_Type_Part_X(index_x, index_y) = Type_X
                Node_Type_Part_Y(index_x, index_y) = Type_Y
             END DO
          ELSE
             ! Allocation Failed
             WRITE(*,*) 'Allocation Error With: Particle Flux'
          END IF
       ELSE IF (line == Header_name(6)) THEN
          ! Readin Flux Boundary Conditions for Electric Potential
          READ(200, *) NFlux_Phi
          ALLOCATE(Node_Type_Phi_X(Nnodes_x,Nnodes_y), &
               Node_Type_Phi_Y(Nnodes_x,Nnodes_y), STAT=STATUS)
          Node_Type_Phi_X = 0
          Node_Type_Phi_Y = 0
          IF (STATUS==0) THEN
             DO I=1, NFlux_Phi
                READ(200,*) index_x, index_y, Type_X, Type_Y
                Node_Type_Phi_X(index_x, index_y) = Type_X
                Node_Type_Phi_Y(index_x, index_y) = Type_Y
             END DO
          ELSE
             ! Allocation Failed
             WRITE(*,*) 'Allocation Error With: Phi Flux'
          END IF
       ELSE IF (line == Header_name(7)) THEN
          ! Exit Input File Flag
          EXIT
       ELSE
          ! Error Analyzing Input
          WRITE(*,*) 'Error Analyzing Domain_Input by Processor'
          CALL MPI_FINALIZE(ierr)
          STOP
       END IF
    END DO
    CLOSE(200)
    IF (SIZE>1) THEN
        ! Send Input Data to Other Processors
        CALL INPUT_BROADCAST
    END IF
END SUBROUTINE MESH_READIN


SUBROUTINE INPUT_BROADCAST
    USE CONTAIN
    IMPLICIT NONE
!-----------------------------------------------------------------------
!*******************MPI BROADCAST*******************
!-----------------------------------------------------------------------
    ! Broadcast Dimensions/ Input Parameters
    CALL MPI_BCAST(Ndimensions,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Stencil_Width,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(RF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(NParts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Print,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Runge_Kutta_Order,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Time_Error_Threshold,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Initial_Time,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Final_Time,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Delta_Time,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    ! Broadcast Node Information
    CALL MPI_BCAST(Nnodes_x,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Nnodes_y,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Nnodes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(x,Nnodes_x,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(y,Nnodes_y,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    ! Broadcast Boundary Conditions Information
    CALL MPI_BCAST(Nboundary,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nodeboundary_x,Nboundary,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nodeboundary_y,Nboundary,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(dofboundary,Nboundary,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(valboundary,Nboundary,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    ! Broadcast Initial Conditions
    CALL MPI_BCAST(Ninitial,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Nodeinitial,Ninitial,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(dofinitial,Ninitial,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(valinitial,Ninitial,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    ! Broadcast Flux Boundary Conditions (Particles & Phi)
    CALL MPI_BCAST(NFlux_Phi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(NFlux_Part,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Node_Type_Part_X,NNodes,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Node_Type_Part_Y,NNodes,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Node_Type_Phi_X,NNodes,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(Node_Type_Phi_Y,NNodes,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
END SUBROUTINE INPUT_BROADCAST
