SUBROUTINE EQUATION_LIBRARY(Nodedof, A_temp, &
     f_temp, inl, jnl, Jacobian_Location)
    USE CONTAIN
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: inl, jnl, Nodedof
    INTEGER, INTENT(INOUT) :: Jacobian_Location(Jacobian_Size)
    REAL(8), INTENT(INOUT) :: A_temp(1, Jacobian_Size)
    REAL(8), INTENT(INOUT) :: f_temp
    
!-----------------------------------------------------------------------
!*******************CALCULATE F VECTOR*********************
!-----------------------------------------------------------------------
    CALL FUNCTION_EVALUATION(inl, jnl, f_temp)
!-----------------------------------------------------------------------
!*******************CALCULATE JACOBIAN*********************
!-----------------------------------------------------------------------
    IF (assembly_check) THEN
       ! Calculate Jacobian First Time
       CALL NUMERICAL_JACOBIAN(inl, jnl, A_temp, &
            f_temp, Jacobian_Location)
    END IF
!-----------------------------------------------------------------------
!*******************UPDATE F VECTOR************************
!-----------------------------------------------------------------------
    ! F Vector Negative Sign for Newton Iteration
    f_temp = -f_temp
END SUBROUTINE EQUATION_LIBRARY


SUBROUTINE FUNCTION_EVALUATION(inl, jnl, f_temp)
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER, INTENT(IN):: inl, jnl
    REAL(8), INTENT(INOUT):: f_temp
    REAL(8):: delta_x_iplushalf, delta_x_iminushalf, delta_x
    REAL(8):: delta_y_jplushalf, delta_y_jminushalf, delta_y
    REAL(8):: Term_1, Term_2, Term_3
    REAL(8):: Ex_minus, Ey_minus, E_minus
    REAL(8):: del_phi_plus, del_Flux_minus, alpha_minus, Term_2_minus_S
    REAL(8):: Flux_i_iplushalf, Flux_i_iminushalf, Flux_i_jplushalf, Flux_i_jminushalf
    REAL(8):: Flux_e_iplushalf, Flux_e_iminushalf, Flux_e_jplushalf, Flux_e_jminushalf
    REAL(8):: Flux_e, Flux_i, Flux_coeff
    INTEGER:: side_x, side_y
!-----------------------------------------------------------------------
!*******************CASE NOTES*****************************
!-----------------------------------------------------------------------
    ! Assuming STRUCTURED MESH
!-----------------------------------------------------------------------
!*******************POISSON EQUATION***********************
!-----------------------------------------------------------------------
    side_x = Node_Type_Phi_X(inl,jnl)
    side_y = Node_Type_Phi_Y(inl,jnl)

    ! x-direction
    IF (side_x == 0) THEN
       ! Grid spacing
       delta_x_iplushalf = x(inl+1)-x(inl)
       delta_x_iminushalf = x(inl)-x(inl-1)
       delta_x = delta_x_iplushalf

       ! Term: d^2(phi)/dx^2 (x-term)
       Term_1 = 2.d0*phi_plus(inl+1,jnl)/(delta_x_iplushalf*(delta_x_iplushalf+delta_x_iminushalf)) &
            - 2.d0*phi_plus(inl,jnl)/(delta_x_iplushalf*delta_x_iminushalf) &
            + 2.d0*phi_plus(inl-1,jnl)/(delta_x_iminushalf*(delta_x_iplushalf+delta_x_iminushalf))

       IF (jnl > 25) THEN
          ! Charge Projection (x-direction)
          del_phi_plus = phi_plus(inl+1,jnl) - phi_plus(inl,jnl)
          CALL ION_FLUX(Flux_i_iplushalf,del_phi_plus,delta_x_iplushalf,mu_i,D_i,&
               n_i_minus(inl,jnl),n_i_minus(inl+1,jnl))
          CALL ELECTRON_FLUX(Flux_e_iplushalf,del_phi_plus,delta_x_iplushalf,mu_e,D_e,&
               n_e_minus(inl,jnl),n_e_minus(inl+1,jnl))
          del_phi_plus = phi_plus(inl,jnl) - phi_plus(inl-1,jnl)
          CALL ION_FLUX(Flux_i_iminushalf,del_phi_plus,delta_x_iminushalf,mu_i,D_i,&
               n_i_minus(inl-1,jnl),n_i_minus(inl,jnl))
          CALL ELECTRON_FLUX(Flux_e_iminushalf,del_phi_plus,delta_x_iminushalf,mu_e,D_e,&
               n_e_minus(inl-1,jnl),n_e_minus(inl,jnl))
       END IF

    ELSE IF (side_x == 1) THEN
       ! Grid spacing
       delta_x_iminushalf = x(inl)-x(inl-1)
       delta_x = delta_x_iminushalf

       ! Term: d^2(phi)/dx^2 (x-term)
       Term_1 = 2.0*phi_plus(inl-1,jnl)/(delta_x_iminushalf*delta_x_iminushalf) &
            - 2.0*phi_plus(inl,jnl)/(delta_x_iminushalf*delta_x_iminushalf)

       IF (jnl > 25) THEN
          ! Charge Projection (x-direction)
          Flux_i_iminushalf = D_i*(n_i_minus(inl-1,jnl)-n_i_minus(inl,jnl))/&
               delta_x_iminushalf
          Flux_e_iminushalf = D_e*(n_e_minus(inl-1,jnl)-n_e_minus(inl,jnl))/&
               delta_x_iminushalf
          Flux_i_iplushalf = -Flux_i_iminushalf
          Flux_e_iplushalf = -Flux_e_iminushalf
       END IF

    ELSE IF (side_x == -1) THEN
       ! Grid spacing
       delta_x_iplushalf = x(inl+1)-x(inl)
       delta_x = delta_x_iplushalf

       ! Term: d^2(phi)/dx^2 (x-term)
       Term_1 = 2.0*phi_plus(inl+1,jnl)/(delta_x_iplushalf*delta_x_iplushalf) &
            - 2.0*phi_plus(inl,jnl)/(delta_x_iplushalf*delta_x_iplushalf)

       IF (jnl > 25) THEN
       ! Charge Projection (x-direction)
          Flux_i_iplushalf = D_i*(n_i_minus(inl,jnl)-n_i_minus(inl+1,jnl))/&
               delta_x_iplushalf
          Flux_e_iplushalf = D_e*(n_e_minus(inl,jnl)-n_e_minus(inl+1,jnl))/&
               delta_x_iplushalf
          Flux_i_iminushalf = -Flux_i_iplushalf
          Flux_e_iminushalf = -Flux_e_iplushalf
       END IF

    ELSE IF (side_x == -2) THEN
       delta_x = 1.d0
    END IF

    ! y-direction
    IF (side_y == 0) THEN
       ! Grid spacing
       delta_y_jplushalf = y(jnl+1)-y(jnl)
       delta_y_jminushalf = y(jnl)-y(jnl-1)
       delta_y = delta_y_jplushalf

       ! Term: d^2(phi)/dy^2 (y-term)
       Term_2 = 2.d0*phi_plus(inl,jnl+1)/(delta_y_jplushalf*(delta_y_jplushalf+delta_y_jminushalf)) &
            - 2.d0*phi_plus(inl,jnl)/(delta_y_jplushalf*delta_y_jminushalf) &
            + 2.d0*phi_plus(inl,jnl-1)/(delta_y_jminushalf*(delta_y_jplushalf+delta_y_jminushalf))
       
       IF (jnl > 25) THEN
          ! Charge Projection (y-direction)
          del_phi_plus = phi_plus(inl,jnl+1) - phi_plus(inl,jnl)
          CALL ION_FLUX(Flux_i_jplushalf,del_phi_plus,delta_y_jplushalf,mu_i,D_i,&
               n_i_minus(inl,jnl),n_i_minus(inl,jnl+1))
          CALL ELECTRON_FLUX(Flux_e_jplushalf,del_phi_plus,delta_y_jplushalf,mu_e,D_e,&
               n_e_minus(inl,jnl),n_e_minus(inl,jnl+1))
          del_phi_plus = phi_plus(inl,jnl) - phi_plus(inl,jnl-1)
          CALL ION_FLUX(Flux_i_jminushalf,del_phi_plus,delta_y_jminushalf,mu_i,D_i,&
               n_i_minus(inl,jnl-1),n_i_minus(inl,jnl))
          CALL ELECTRON_FLUX(Flux_e_jminushalf,del_phi_plus,delta_y_jminushalf,mu_e,D_e,&
               n_e_minus(inl,jnl-1),n_e_minus(inl,jnl))    
       END IF

    ELSE IF (side_y == 1) THEN
       ! Grid spacing
       delta_y_jminushalf = y(jnl)-y(jnl-1)
       delta_y = delta_y_jminushalf

       ! Term: d^2(phi)/dy^2 (y-term)
       Term_2 = 2.0*phi_plus(inl,jnl-1)/(delta_y_jminushalf*delta_y_jminushalf) &
            - 2.0*phi_plus(inl,jnl)/(delta_y_jminushalf*delta_y_jminushalf) 

       IF (jnl > 25) THEN
          ! Charge Projection (y-direction)
          Flux_i_jminushalf = D_i*(n_i_minus(inl,jnl-1)-n_i_minus(inl,jnl))/&
               delta_y_jminushalf
          Flux_e_jminushalf = D_e*(n_e_minus(inl,jnl-1)-n_e_minus(inl,jnl))/&
               delta_y_jminushalf
          Flux_i_jplushalf = -Flux_i_jminushalf
          Flux_e_jplushalf = -Flux_e_jminushalf
       END IF

    ELSE IF (side_y == 2) THEN
       ! Grid spacing
       delta_y_jplushalf = y(jnl+1)-y(jnl)
       delta_y_jminushalf = y(jnl)-y(jnl-1)
       delta_y = delta_y_jplushalf

       Term_2 = -(eps_d/delta_y_jminushalf)*phi_plus(inl,jnl-1) &
            +(eps_p/delta_y_jplushalf+eps_d/delta_y_jminushalf)*phi_plus(inl,jnl) &
            -(eps_p/delta_y_jplushalf)*phi_plus(inl,jnl+1)
    END IF

    ! Source Term: e/epsilon*(n_i-n_e)
    IF (jnl == 26) THEN
       Term_3 = -(L_0**2.)*n_0*E_charge/(Epsilon*Phi_0)*sigma_charge(inl,jnl)       
    ELSE IF (jnl > 25) THEN
       Term_3 = (L_0**2.)*n_0*E_charge/(Epsilon*Phi_0)*(n_i_plus(inl,jnl)-n_e_plus(inl,jnl) + &
            delta_time*(2.0*(Flux_e_iplushalf - Flux_e_iminushalf)/(delta_x_iplushalf+delta_x_iminushalf)+&
            2.0*(Flux_e_jplushalf - Flux_e_jminushalf)/(delta_y_jplushalf+delta_y_jminushalf)-&
            2.0*(Flux_i_iplushalf - Flux_i_iminushalf)/(delta_x_iplushalf+delta_x_iminushalf)-&
            2.0*(Flux_i_jplushalf - Flux_i_jminushalf)/(delta_y_jplushalf+delta_y_jminushalf)))
    ELSE
       Term_3 = 0.d0
    END IF

    ! Calculate F
    f_temp = (Term_1+Term_2+Term_3)*delta_x*delta_y

END SUBROUTINE FUNCTION_EVALUATION


SUBROUTINE NUMERICAL_JACOBIAN(inl, jnl, A_temp, f_temp_initial, &
     Jacobian_Location)
    USE CONTAIN
    IMPLICIT NONE
    INTEGER, INTENT(IN):: inl, jnl
    INTEGER, INTENT(INOUT) :: Jacobian_Location(Jacobian_Size)
    REAL(8), INTENT(IN):: f_temp_initial
    REAL(8), INTENT(INOUT):: A_temp(1, Jacobian_Size)
    LOGICAL:: Zero_Perturb
    REAL(8):: Perturb, f_temp_perturb
    REAL(8):: Temporary_Storage
    INTEGER:: I, J, K

    ! Initialize
    Jacobian_index = 0
    Temporary_Storage = 0
    Perturb = 0.0001
    DO K=1, stencil_width
       IF (K==1) THEN
          I = inl-1
          J = jnl
       ELSE IF (K==2) THEN
          I = inl+1
          J = jnl           
       ELSE IF (K==3) THEN
          I = inl
          J = jnl
       ELSE IF (K==4) THEN
          I = inl
          J = jnl+1
       ELSE IF (K==5) THEN
          I = inl
          J = jnl-1
       END IF
       Zero_Perturb = .False.
!-----------------------------------------------------------------------
!*******************JACOBIAN SETUP*************************
!-----------------------------------------------------------------------
       Jacobian_Index = Jacobian_Index + 1
       Temporary_Storage = phi_plus(I,J)
       IF (Abs(phi_plus(I,J)) > 0) THEN
          phi_plus(I,J) = phi_plus(I,J) + &
               phi_plus(I,J)*Perturb
       ELSE
          Zero_Perturb = .True.
          phi_plus(I,J) = Perturb
       END IF
       CALL FUNCTION_EVALUATION(inl, jnl, f_temp_perturb)
       IF (.NOT. Zero_Perturb) THEN
          phi_plus(I,J) = Temporary_Storage
       ELSE
          phi_plus(I,J) = 1
       END IF
       CALL JACOBIAN_EVALUATION(I, J, f_temp_perturb, &
            f_temp_initial, Jacobian_Location, A_temp, &
            Perturb)
       IF (Zero_Perturb) THEN
          phi_plus(I,J) = Temporary_Storage
       END IF
    END DO
END SUBROUTINE NUMERICAL_JACOBIAN


SUBROUTINE JACOBIAN_EVALUATION(I, J, f_temp_perturb, f_temp_initial, &
     Jacobian_Location, A_temp, Perturb)
    USE CONTAIN
    IMPLICIT NONE
    INTEGER, INTENT(IN):: I, J
    REAL(8), INTENT(IN):: f_temp_perturb, f_temp_initial, Perturb
    REAL(8), INTENT(INOUT):: A_temp(1, Jacobian_Size)
    INTEGER, INTENT(INOUT):: Jacobian_Location(Jacobian_Size)
    
    A_temp(1, Jacobian_Index) = (f_temp_perturb - &
         f_temp_initial)/(phi_plus(I,J)*Perturb)
    Jacobian_Location(Jacobian_Index) = &
         Globaldof(I,J)-1
END SUBROUTINE JACOBIAN_EVALUATION
