SUBROUTINE INITIALIZE
    ! This Subroutine sets up the degree of freedom array, finds
    ! the number of dofs, and the number of equations.
    USE CONTAIN
    IMPLICIT NONE
    INTEGER:: I, J, iter_temp
!-----------------------------------------------------------------------
!*******************SETUP ARRAYS***************************
!-----------------------------------------------------------------------
    ! Setup 2-D N_i
    Allocate(n_i_plus(Nnodes_x, Nnodes_y), n_i_minus(Nnodes_x, Nnodes_y), &
         n_i_minus_orig(Nnodes_x, Nnodes_y))
    ! Setup 2-D N_e
    Allocate(n_e_plus(Nnodes_x, Nnodes_y), n_e_minus(Nnodes_x, Nnodes_y), &
         n_e_minus_orig(Nnodes_x, Nnodes_y))
    ! Setup 2-D Phi
    Allocate(phi_plus(Nnodes_x, Nnodes_y))
    ! Setup Surface_Charge
    Allocate(Sigma_Charge(Nnodes_x, Nnodes_y))
    ! Temporal Integration Terms
    Allocate(K_i(Runge_Kutta_Order,Nnodes_x, Nnodes_y), &
         K_e(Runge_Kutta_Order,Nnodes_x, Nnodes_y))
    Allocate(Error_Offset_n_i(Nnodes_x, Nnodes_y), &
         Error_Offset_n_e(Nnodes_x, Nnodes_y))
!-----------------------------------------------------------------------
!*******************SETUP DOF MATRICES*********************
!-----------------------------------------------------------------------
    ! Setup Matrices for First Timestep
    Ndof = Nnodes
    ALLOCATE(localdof(Nnodes_x, Nnodes_y), globaldof(Nnodes_x, Nnodes_y))
    localdof = 0
    globaldof = 0
    iter_temp = 0
    Sigma_Charge = 0
    K_i = 0
    K_e = 0
    Error_Offset_n_i = 0
    Error_Offset_n_e = 0
!-----------------------------------------------------------------------
!*******************BOUNDARY NODES*************************
!-----------------------------------------------------------------------		
    DO J=1, Nnodes_y
       DO I=1, Nnodes_x
          ! Specify Initial Conditions Via Input (ValInitial)
          iter_temp = iter_temp + 1
          n_i_minus(I,J) = Valinitial(iter_temp)
          n_i_minus_orig(I,J) = Valinitial(iter_temp)                
          n_i_plus(I,J) = Valinitial(iter_temp)
       END DO
    END DO
    DO J=1, Nnodes_y
       DO I=1, Nnodes_x
          ! Specify Initial Conditions Via Input (ValInitial)
          iter_temp = iter_temp + 1
          n_e_minus(I,J) = Valinitial(iter_temp)
          n_e_minus_orig(I,J) = Valinitial(iter_temp)
          n_e_plus(I,J) = Valinitial(iter_temp)
       END DO
    END DO
    DO J=1, Nnodes_y
       DO I=1, Nnodes_x
          ! Specify Initial Conditions Via Input (ValInitial)
          iter_temp = iter_temp + 1
          phi_plus(I,J) = Valinitial(iter_temp)
       END DO
    END DO
    DO J=1, Nnodes_y
       DO I=1, Nnodes_x
          ! Specify Initial Conditions Via Input (ValInitial)
          iter_temp = iter_temp + 1
          Sigma_Charge(I,J) = Valinitial(iter_temp)
       END DO
    END DO
    DO I = 1, Nboundary 
       ! Setup Globaldof array (Boundary Conditions Set As Neg.)
       Localdof(Nodeboundary_x(I), Nodeboundary_y(I)) = -I
       phi_plus(Nodeboundary_x(I), Nodeboundary_y(I)) = Valboundary(I)
    END DO
    ! Calculate Number of Eqns. (No Boundary Conditions)
    Neqn = 0
    DO J = 1, Nnodes_y
       DO I = 1, Nnodes_x
          IF (Localdof(I,J) == 0) THEN
             ! Ignores Negatives (BCs) In Calculating # of Eqns.
             Neqn = Neqn + 1
             Localdof(I,J) = Neqn
          END IF
       END DO
    END DO
!-----------------------------------------------------------------------
!*******************NONZEROS PER ROW***********************
!-----------------------------------------------------------------------
    ! Assuming Maximum DOF is on Node 1
    Max_Dof_of_Node = 1
    Nonzeros_per_row = Stencil_Width
    Jacobian_Size = Stencil_Width
!-----------------------------------------------------------------------
!*******************GLOBAL DOF CALCULATION*****************
!-----------------------------------------------------------------------
    Globaldof = Localdof
    Neqn_Global = Neqn
END SUBROUTINE INITIALIZE

SUBROUTINE PETSC
    ! This Subroutine Assembles Matrices in PETSc and Solves Them
    USE CONTAIN
    IMPLICIT NONE
!-----------------------------------------------------------------------
!*******************RELEVANT VARIABLES*********************
!-----------------------------------------------------------------------
    INTEGER:: Low_Value, High_Value
    INTEGER:: I, J, K, Inl, jnl, side_y
    REAL(8):: X_s(neqn), Val, Nodex(Ndimensions)
    LOGICAL:: Steady_State
    REAL(8):: f_temp
    REAL(8):: A_temp(1, Jacobian_Size)
    INTEGER:: Nodedof, Jacobian_Location(Jacobian_Size)
    LOGICAL:: Flux_Check

!-----------------------------------------------------------------------
!*******************A MATRIX SETUP*************************
!-----------------------------------------------------------------------
    IF (assembly_check) THEN
       CALL MatCreate(comm, A_pet, ierr)
       CALL MatSetSizes(A_pet, neqn, neqn, neqn_global, neqn_global, ierr)
       CALL MatSetFromOptions(A_pet, ierr)
       !IF (SIZE>1) THEN
       ! Parallel Solver
       !ELSE
       ! Single Processor Solver
       CALL MatSeqAIJSetPreallocation(A_pet, nonzeros_per_row, &
            PETSC_NULL_INTEGER, ierr)
       I_start = 0
       I_end = neqn
       !END IF
       CALL MatSetOption(A_pet, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)
    END IF
!-----------------------------------------------------------------------
!*******************GLOBALF VECTOR SETUP*******************
!-----------------------------------------------------------------------	
    ! Note VecSetValues() Uses 0-Base Row/Column Numbers
    CALL VecCreate(comm, globalf_pet, ierr)
    CALL VecSetSizes(globalf_pet, neqn, neqn_global, ierr)
    CALL VecSetFromOptions(globalf_pet, ierr)
    !IF (SIZE>1) THEN
    ! Parallel Solver
    !ELSE
    ! Single Processor Solver
    low_value = 0
    high_value = neqn
    !END IF 
    CALL VecSetOption(globalf_pet, VEC_IGNORE_NEGATIVE_INDICES, &
         PETSC_TRUE, ierr)
!-----------------------------------------------------------------------
!*******************MATRIX ASSEMBLY************************
!-----------------------------------------------------------------------
    DO jnl = 1, Nnodes_y
       DO inl =1, Nnodes_x
          Nodedof = Globaldof(Inl,Jnl)-1
          IF (Nodedof >= 0) THEN
             ! Initialize Arrays  
             A_temp = 0
             f_temp = 0
             Jacobian_Location = -5

             side_y = Node_type_phi_Y(inl,jnl)
             IF (side_y == 2 .or. side_y == -2) THEN
                CALL SURFACE_CHARGE_EVALUATION(inl,jnl)
             END IF
             
             CALL EQUATION_LIBRARY(Nodedof, A_temp, f_temp, &
                  inl, jnl, Jacobian_Location) 
             IF (assembly_check) THEN
                CALL MatSetValues(A_pet, 1, Nodedof, Jacobian_Size, &
                     Jacobian_Location, A_temp(1,:), INSERT_VALUES, ierr)
             END IF
             CALL VecSetValues(globalf_pet, 1, Nodedof, &
                  f_temp, INSERT_VALUES, ierr)
          END IF
       END DO
    END DO
    IF (assembly_check) THEN
       ! Assemble A Matrix
       CALL MatAssemblyBegin(A_pet, MAT_FINAL_ASSEMBLY, ierr)
       CALL MatAssemblyEnd(A_pet, MAT_FINAL_ASSEMBLY, ierr)
    END IF
    ! Assemble Globalf Matrix
    CALL VecAssemblyBegin(globalf_pet, ierr)
    CALL VecAssemblyEnd(globalf_pet, ierr)
    CALL MatView(A_pet,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !CALL VecView(globalf_pet,PETSC_VIEWER_STDOUT_WORLD,ierr)
    read(*,*) inl
!-----------------------------------------------------------------------
!*******************SOLVER*********************************
!-----------------------------------------------------------------------
    IF (assembly_check) THEN
       CALL KSPCreate(comm, ksp, ierr)
       ! Set Options From Runtime (e.g. -ksp_monitor, -ksp_rtol <rtol>, ...
       CALL KSPSetFromOptions(ksp, ierr)
       ! Solve Newton System
       CALL KSPSetOperators(ksp, A_pet, A_pet, SAME_NONZERO_PATTERN, ierr)
    END IF
    CALL KSPSolve(ksp, globalf_pet, globalf_pet, ierr)
!-----------------------------------------------------------------------
!*******************READ INTO X_S**************************
!-----------------------------------------------------------------------
    x_s = 0.d0
    J = 0
    DO I = low_value, high_value-1
       J = J + 1
       CALL VecGetValues(globalf_pet, 1, I, x_s(J), ierr)
    END DO
!-----------------------------------------------------------------------
!*******************BOUNDARY CONDITIONS********************
!-----------------------------------------------------------------------
    ! RF Boundary Conditions (Insert If Needed)
    ! Prescribed Boundary Conditions
    DO I = 1, Nboundary
       phi_plus(nodeboundary_x(i),nodeboundary_y(i)) =  ValBoundary(I)
    END DO
!-----------------------------------------------------------------------
!*******************ROOT UPDATE****************************
!-----------------------------------------------------------------------
    DO J = 1, Nnodes_y
       DO I = 1, Nnodes_x
          IF (Localdof(I,J) > 0) THEN
             phi_plus(I,J) = Phi_plus(I,J) + x_s(localdof(I,J))
          END IF
       END DO
    END DO
    !CALL SendSolution
!-----------------------------------------------------------------------
!*******************DESTROY OBJECTS************************
!-----------------------------------------------------------------------
    !CALL KSPDestroy(ksp, ierr)
    !CALL MatDestroy(A_pet, ierr)
    CALL VecDestroy(globalf_pet, ierr)
    assembly_check = .false.
END SUBROUTINE PETSC

SUBROUTINE EXPLICIT_SOLVE
    ! This Subroutine Establishes Explicit Variables to be Solved
    USE CONTAIN
    IMPLICIT NONE
    INTEGER:: inl, jnl, stage
!-----------------------------------------------------------------------
!*******************SOLVE FLUX*****************************
!----------------------------------------------------------------------- 
    DO stage = 1, Runge_Kutta_Order
!-----------------------------------------------------------------------
!*******************SOLVE DOMAIN*****************************
!-----------------------------------------------------------------------
       DO jnl = 26, Nnodes_y
          DO inl = 1, Nnodes_x
             CALL EQUATION_LIBRARY_EXPLICIT(inl, jnl,stage)
          END DO
       END DO
       CALL RUNGE_KUTTA_EVALUATION(stage)
    END DO
END SUBROUTINE EXPLICIT_SOLVE

