SUBROUTINE EQUATION_LIBRARY_EXPLICIT(inl, jnl, stage)
    USE CONTAIN
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: inl, jnl, stage

!-----------------------------------------------------------------------
!*******************CALL EVALUATION************************
!-----------------------------------------------------------------------
    CALL FUNCTION_EVALUATION_EXPLICIT(inl, jnl, stage)
    
END SUBROUTINE EQUATION_LIBRARY_EXPLICIT


SUBROUTINE RUNGE_KUTTA_EVALUATION(stage)
    USE CONTAIN
    IMPLICIT NONE
    INTEGER, INTENT(IN):: stage

    IF (Runge_Kutta_Order == 1) THEN
       ! Euler Scheme
       n_i_plus = n_i_minus_orig + K_i(1,:,:)*delta_time
       n_e_plus = n_e_minus_orig + K_e(1,:,:)*delta_time
    ELSE IF (Runge_Kutta_Order == 2) THEN
       ! 2nd Order
       IF (stage == 1) THEN
          n_i_plus = n_i_minus_orig + K_i(1,:,:)*delta_time/2.d0
          n_e_plus = n_e_minus_orig + K_e(1,:,:)*delta_time/2.d0
          n_i_minus = n_i_minus_orig + K_i(1,:,:)*delta_time
          n_e_minus = n_e_minus_orig + K_e(1,:,:)*delta_time
       ELSE
          n_i_plus = n_i_plus + K_i(2,:,:)*delta_time/2.d0
          n_e_plus = n_e_plus + K_e(2,:,:)*delta_time/2.d0
       END IF
    ELSE IF (Runge_Kutta_Order == 3) THEN
       ! 3rd Order
       IF (stage == 1) THEN
       ELSE IF (stage == 2) THEN
       ELSE 
       END IF
    ELSE IF (Runge_Kutta_Order == 4) THEN
       ! 4th Order
       IF (stage == 1) THEN
          n_i_plus = n_i_minus_orig + K_i(1,:,:)*delta_time/6.d0
          n_e_plus = n_e_minus_orig + K_e(1,:,:)*delta_time/6.d0
          n_i_minus = n_i_minus_orig + K_i(1,:,:)*delta_time/2.d0
          n_e_minus = n_e_minus_orig + K_e(1,:,:)*delta_time/2.d0               
       ELSE IF (stage == 2) THEN
          n_i_plus = n_i_plus + K_i(2,:,:)*delta_time/3.d0
          n_e_plus = n_e_plus + K_e(2,:,:)*delta_time/3.d0
          n_i_minus = n_i_minus_orig + K_i(2,:,:)*delta_time/2.d0
          n_e_minus = n_e_minus_orig + K_e(2,:,:)*delta_time/2.d0        
       ELSE IF (stage == 3) THEN
          n_i_plus = n_i_plus + K_i(3,:,:)*delta_time/3.d0
          n_e_plus = n_e_plus + K_e(3,:,:)*delta_time/3.d0
          n_i_minus = n_i_minus_orig + K_i(3,:,:)*delta_time
          n_e_minus = n_e_minus_orig + K_e(3,:,:)*delta_time           
       ELSE
          n_i_plus = n_i_plus + K_i(4,:,:)*delta_time/6.d0
          n_e_plus = n_e_plus + K_e(4,:,:)*delta_time/6.d0       
       END IF
    ELSE IF (Runge_Kutta_Order == 7) THEN
       ! Dormand-Prince Method (Adaptive Algorithm)
       IF (stage == 1) THEN
          n_i_minus = n_i_minus_orig + K_i(1,:,:)*delta_time/5.d0
          n_e_minus = n_e_minus_orig + K_e(1,:,:)*delta_time/5.d0          
       ELSE IF (stage == 2) THEN
          n_i_minus = n_i_minus_orig + delta_time*(3.d0*K_i(1,:,:)/40.d0&
               +9.d0/40.d0*K_i(2,:,:))
          n_e_minus = n_e_minus_orig + delta_time*(3.d0*K_e(1,:,:)/40.d0&
               +9.d0/40.d0*K_e(2,:,:))        
       ELSE IF (stage == 3) THEN
          n_i_minus = n_i_minus_orig + delta_time*(44.d0*K_i(1,:,:)/45.d0&
               -56.d0/15.d0*K_i(2,:,:)+32.d0/9.d0*K_i(3,:,:))
          n_e_minus = n_e_minus_orig + delta_time*(44.d0*K_e(1,:,:)/45.d0&
               -56.d0/15.d0*K_e(2,:,:)+32.d0/9.d0*K_e(3,:,:))          
       ELSE IF (stage == 4) THEN
          n_i_minus = n_i_minus_orig + delta_time*(19372.d0*K_i(1,:,:)/6561.d0&
               -25360.d0/2187.d0*K_i(2,:,:) + 64448.d0/6561.d0*K_i(3,:,:)&
               -212.d0/729.d0*K_i(4,:,:))
          n_e_minus = n_e_minus_orig + delta_time*(19372.d0*K_e(1,:,:)/6561.d0&
               -25360.d0/2187.d0*K_e(2,:,:) + 64448.d0/6561.d0*K_e(3,:,:)&
               -212.d0/729.d0*K_e(4,:,:))           
       ELSE IF (stage == 5) THEN
          n_i_minus = n_i_minus_orig + delta_time*(9017.d0*K_i(1,:,:)/3168.d0&
               -355.d0/33.d0*K_i(2,:,:) + 46732.d0/5247.d0*K_i(3,:,:)&
               +49.d0/176.d0*K_i(4,:,:) - 5103.d0/18656.d0*K_i(5,:,:))
          n_e_minus = n_e_minus_orig + delta_time*(9017.d0*K_e(1,:,:)/3168.d0&
               -355.d0/33.d0*K_e(2,:,:) + 46732.d0/5247.d0*K_e(3,:,:)&
               +49.d0/176.d0*K_e(4,:,:) - 5103.d0/18656.d0*K_e(5,:,:))        
       ELSE IF (stage == 6) THEN
          n_i_minus = n_i_minus_orig + delta_time*(35.d0*K_i(1,:,:)/384.d0&
               +500.d0/1113.d0*K_i(3,:,:) + 125.d0/192.d0*K_i(4,:,:)&
               -2187.d0/6784.d0*K_i(5,:,:) + 11.d0/84.d0*K_i(6,:,:))
          n_e_minus = n_e_minus_orig + delta_time*(35.d0*K_e(1,:,:)/384.d0&
               +500.d0/1113.d0*K_e(3,:,:) + 125.d0/192.d0*K_e(4,:,:)&
               -2187.d0/6784.d0*K_e(5,:,:) + 11.d0/84.d0*K_e(6,:,:))           
       ELSE
          ! 4th Order
          n_i_plus = n_i_minus_orig + delta_time*(35.d0*K_i(1,:,:)/384.d0&
               +500.d0/1113.d0*K_i(3,:,:) + 125.d0/192.d0*K_i(4,:,:)&
               -2187.d0/6784.d0*K_i(5,:,:) + 11.d0/84.d0*K_i(6,:,:))
          n_e_plus = n_e_minus_orig + delta_time*(35.d0*K_e(1,:,:)/384.d0&
               +500.d0/1113.d0*K_e(3,:,:) + 125.d0/192.d0*K_e(4,:,:)&
               -2187.d0/6784.d0*K_e(5,:,:) + 11.d0/84.d0*K_e(6,:,:))
          ! 5th Order
          Error_Offset_n_i = n_i_minus_orig + delta_time*(5179.d0*K_i(1,:,:)/57600.d0&
               +7571.d0/16695.d0*K_i(3,:,:) + 393.d0/640.d0*K_i(4,:,:)&
               -92097.d0/339200.d0*K_i(5,:,:) +  187.d0/2100.d0*K_i(6,:,:)&
               +1.d0/40.d0*K_i(7,:,:))
          Error_Offset_n_e = n_e_minus_orig + delta_time*(5179.d0*K_e(1,:,:)/57600.d0&
               +7571.d0/16695.d0*K_e(3,:,:) + 393.d0/640.d0*K_e(4,:,:)&
               -92097.d0/339200.d0*K_e(5,:,:) + 187.d0/2100.d0*K_e(6,:,:)&
               +1.d0/40.d0*K_e(7,:,:))
          ! Find Maximum Error Offset
          Error_Offset_n_i = abs(Error_Offset_n_i - n_i_plus)
          Error_Offset_n_e = abs(Error_Offset_n_e - n_e_plus)
          Error_Offset = MAXVAL(Error_Offset_n_e)
          IF (Error_Offset < MAXVAL(Error_Offset_n_i)) THEN
             Error_Offset = MAXVAL(Error_Offset_n_i)
          END IF
          ! Scale Adaptive Time Step
          Scaling_Factor = (time_error_threshold*delta_time/(2*&
               Error_Offset))**(0.2d0)
          !WRITE(*,*) 'Initial delta_time:', delta_time
          delta_time = Scaling_Factor*delta_time
          !WRITE(*,*) 'Second delta_time:', delta_time            
          !READ(*,*) iter
       END IF
    END IF
END SUBROUTINE RUNGE_KUTTA_EVALUATION


SUBROUTINE FUNCTION_EVALUATION_EXPLICIT(inl, jnl, stage)
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER, INTENT(IN):: inl, jnl, stage
    REAL(8):: delta_x_iplushalf, delta_x_iminushalf
    REAL(8):: delta_y_jplushalf, delta_y_jminushalf
    REAL(8):: Ex_plus, Ey_plus
    REAL(8):: Flux_i_iplushalf, Flux_i_iminushalf, Flux_i_jplushalf, Flux_i_jminushalf
    REAL(8):: Flux_e_iplushalf, Flux_e_iminushalf, Flux_e_jplushalf, Flux_e_jminushalf
    REAL(8):: Flux_e, Flux_i, Flux_coeff
    REAL(8):: del_phi_plus, E_plus, alpha_minus, Term_2_minus_S
    REAL(8):: del_Flux_minus_i, del_Flux_minus_e
    INTEGER:: side_x, side_y

!-----------------------------------------------------------------------
!*******************CASE NOTES*****************************
!-----------------------------------------------------------------------
    ! Assuming STRUCTURED MESH
!-----------------------------------------------------------------------
!*******************PARTICLE CONTINUITY********************
!-----------------------------------------------------------------------
    side_x = Node_Type_Part_X(inl,jnl)
    side_y = Node_Type_Part_Y(inl,jnl)


    ! x-direction
    IF (side_x == 0) THEN
       ! Grid Spacing
       delta_x_iplushalf = x(inl+1)-x(inl)
       delta_x_iminushalf = x(inl)-x(inl-1)

       ! Evaluate Fluxes (in x-direction)
       del_phi_plus = phi_plus(inl+1,jnl) - phi_plus(inl,jnl)
       Ex_plus = -del_phi_plus/(delta_x_iplushalf)
       CALL ION_FLUX(Flux_i_iplushalf,del_phi_plus,delta_x_iplushalf,mu_i,D_i,&
            n_i_minus(inl,jnl),n_i_minus(inl+1,jnl))
       CALL ELECTRON_FLUX(Flux_e_iplushalf,del_phi_plus,delta_x_iplushalf,mu_e,D_e,&
            n_e_minus(inl,jnl),n_e_minus(inl+1,jnl))
    
       del_phi_plus = phi_plus(inl,jnl) - phi_plus(inl-1,jnl)
       Ex_plus = (Ex_plus + -del_phi_plus/delta_x_iminushalf)/2
       CALL ION_FLUX(Flux_i_iminushalf,del_phi_plus,delta_x_iminushalf,mu_i,D_i,&
            n_i_minus(inl-1,jnl),n_i_minus(inl,jnl))
       CALL ELECTRON_FLUX(Flux_e_iminushalf,del_phi_plus,delta_x_iminushalf,mu_e,D_e,&
            n_e_minus(inl-1,jnl),n_e_minus(inl,jnl))
    
       del_Flux_minus_i = 2.0*(Flux_i_iplushalf - Flux_i_iminushalf)/&
            (delta_x_iplushalf+delta_x_iminushalf)
       del_Flux_minus_e = 2.0*(Flux_e_iplushalf - Flux_e_iminushalf)/&
            (delta_x_iplushalf+delta_x_iminushalf)

    ELSE IF (side_x == 1) THEN
       ! Grid Spacing
       ! Right
       delta_x_iminushalf = x(inl)-x(inl-1)

       ! Evaluate Fluxes (in x-direction)
       Ex_plus = 0.d0
       Flux_i_iminushalf = D_i*(n_i_minus(inl-1,jnl)-n_i_minus(inl,jnl))/&
            delta_x_iminushalf
       Flux_e_iminushalf = D_e*(n_e_minus(inl-1,jnl)-n_e_minus(inl,jnl))/&
            delta_x_iminushalf

       Flux_i_iplushalf = -Flux_i_iminushalf
       Flux_e_iplushalf = -Flux_e_iminushalf

       del_Flux_minus_i = (Flux_i_iplushalf - Flux_i_iminushalf)/delta_x_iminushalf
       del_Flux_minus_e = (Flux_e_iplushalf - Flux_e_iminushalf)/delta_x_iminushalf
    ELSE IF (side_x == -1) THEN
       ! Grid Spacing
       ! Left
       delta_x_iplushalf = x(inl+1)-x(inl)

       ! Evaluate Fluxes (in x-direction)
       Ex_plus = 0.d0
       Flux_i_iplushalf = D_i*(n_i_minus(inl,jnl)-n_i_minus(inl+1,jnl))/&
            delta_x_iplushalf
       Flux_e_iplushalf = D_e*(n_e_minus(inl,jnl)-n_e_minus(inl+1,jnl))/&
            delta_x_iplushalf

       Flux_i_iminushalf = -Flux_i_iplushalf
       Flux_e_iminushalf = -Flux_e_iplushalf
       
       del_Flux_minus_i = (Flux_i_iplushalf - Flux_i_iminushalf)/delta_x_iplushalf
       del_Flux_minus_e = (Flux_e_iplushalf - Flux_e_iminushalf)/delta_x_iplushalf       
    END IF


    ! y-direction
    IF (side_y == 0) THEN
       ! Grid Spacing
       delta_y_jplushalf = y(jnl+1)-y(jnl)
       delta_y_jminushalf = y(jnl)-y(jnl-1)

       ! Evaluate Fluxes (in y-direction)
       del_phi_plus = phi_plus(inl,jnl+1) - phi_plus(inl,jnl)
       Ey_plus = -del_phi_plus/(delta_y_jplushalf)
       CALL ION_FLUX(Flux_i_jplushalf,del_phi_plus,delta_y_jplushalf,mu_i,D_i,&
            n_i_minus(inl,jnl),n_i_minus(inl,jnl+1))
       CALL ELECTRON_FLUX(Flux_e_jplushalf,del_phi_plus,delta_y_jplushalf,mu_e,D_e,&
            n_e_minus(inl,jnl),n_e_minus(inl,jnl+1))    
       
       del_phi_plus = phi_plus(inl,jnl) - phi_plus(inl,jnl-1)
       Ey_plus = (Ey_plus + -del_phi_plus/delta_y_jminushalf)/2
       CALL ION_FLUX(Flux_i_jminushalf,del_phi_plus,delta_y_jminushalf,mu_i,D_i,&
            n_i_minus(inl,jnl-1),n_i_minus(inl,jnl))
       CALL ELECTRON_FLUX(Flux_e_jminushalf,del_phi_plus,delta_y_jminushalf,mu_e,D_e,&
            n_e_minus(inl,jnl-1),n_e_minus(inl,jnl))
    
       del_Flux_minus_i = del_Flux_minus_i + &
            2.0*(Flux_i_jplushalf - Flux_i_jminushalf)/&
            (delta_y_jplushalf+delta_y_jminushalf)
       del_Flux_minus_e = del_Flux_minus_e + &
            2.0*(Flux_e_jplushalf - Flux_e_jminushalf)/&
            (delta_y_jplushalf+delta_y_jminushalf)
    ELSE IF (side_y == 1) THEN
       ! Grid Spacing
       delta_y_jminushalf = y(jnl)-y(jnl-1)

       ! Evaluate Fluxes (in y-direction)
       Ey_plus = 0.d0
       Flux_i_jminushalf = D_i*(n_i_minus(inl,jnl-1)-n_i_minus(inl,jnl))/&
            delta_y_jminushalf
       Flux_e_jminushalf = D_e*(n_e_minus(inl,jnl-1)-n_e_minus(inl,jnl))/&
            delta_y_jminushalf

       Flux_i_jplushalf = -Flux_i_jminushalf
       Flux_e_jplushalf = -Flux_e_jminushalf

       del_Flux_minus_i = del_Flux_minus_i + &
            (Flux_i_jplushalf - Flux_i_jminushalf)/delta_y_jminushalf
       del_Flux_minus_e = del_Flux_minus_e + &
            (Flux_e_jplushalf - Flux_e_jminushalf)/delta_y_jminushalf
    ELSE IF (side_y == -2) THEN
       ! Grid Spacing
       delta_y_jplushalf = y(jnl+1)-y(jnl)

       ! Evaluate Fluxes (in y-direction)
       del_phi_plus = phi_plus(inl,jnl+1) - phi_plus(inl,jnl)
       Ey_plus = -del_phi_plus/delta_y_jplushalf
       CALL ION_FLUX(Flux_i_jplushalf,del_phi_plus,delta_y_jplushalf,mu_i,D_i,&
            n_i_minus(inl,jnl),n_i_minus(inl,jnl+1))
       CALL ELECTRON_FLUX(Flux_e_jplushalf,del_phi_plus,delta_y_jplushalf,mu_e,D_e,&
            n_e_minus(inl,jnl),n_e_minus(inl,jnl+1)) 

       ! Dot product: n.Gamma
       ! Assuming Anode at bottom side
       IF (del_phi_plus/delta_y_jplushalf > 0) THEN
          ! Defines alpha in flux BC
          Flux_Coeff = 1.d0
       ELSE
          Flux_Coeff = 0.d0
       END IF
       Flux_i = -1/4.*v_i*n_i_minus(inl,jnl)-Flux_Coeff*mu_i*&
            n_i_minus(inl,jnl)*del_phi_plus/delta_y_jplushalf
       Flux_e = -1/4.*v_e*n_e_minus(inl,jnl)+(1.d0-Flux_Coeff)*mu_e*&
            n_e_minus(inl,jnl)*del_phi_plus/delta_y_jplushalf
       Flux_e_jminushalf = 2.d0*Flux_e - Flux_e_jplushalf

       del_Flux_minus_i = del_Flux_minus_i + &
            (Flux_i_jplushalf-Flux_i)/(1/2.*delta_y_jplushalf)
       del_Flux_minus_e = del_Flux_minus_e + &
            (Flux_e_jplushalf-Flux_e)/(1/2.*delta_y_jplushalf)
    ELSE IF (side_y == -3) THEN
       ! Grid Spacing
       delta_y_jplushalf = y(jnl+1)-y(jnl)

       ! Evaluate Fluxes (in y-direction)
       del_phi_plus = phi_plus(inl,jnl+1) - phi_plus(inl,jnl)
       Ey_plus = -del_phi_plus/delta_y_jplushalf
       CALL ION_FLUX(Flux_i_jplushalf,del_phi_plus,delta_y_jplushalf,mu_i,D_i,&
            n_i_minus(inl,jnl),n_i_minus(inl,jnl+1))
       CALL ELECTRON_FLUX(Flux_e_jplushalf,del_phi_plus,delta_y_jplushalf,mu_e,D_e,&
            n_e_minus(inl,jnl),n_e_minus(inl,jnl+1)) 
       
       ! Dot product: n.Gamma
       ! Assuming Anode at bottom side
       IF (del_phi_plus/delta_y_jplushalf > 0) THEN
          ! Defines alpha in flux BC
          Flux_Coeff = 1.d0
       ELSE
          Flux_Coeff = 0.d0
       END IF
       Flux_i = -1/4.*v_i*n_i_minus(inl,jnl)-Flux_Coeff*mu_i*&
            n_i_minus(inl,jnl)*del_phi_plus/delta_y_jplushalf
       Flux_e = -1/4.*v_e*n_e_minus(inl,jnl)- gamma*Flux_i + &
            (1.d0-Flux_Coeff)*mu_e*n_e_minus(inl,jnl)*&
            del_phi_plus/delta_y_jplushalf
       Flux_e_jminushalf = 2.d0*Flux_e - Flux_e_jplushalf

       del_Flux_minus_i = del_Flux_minus_i + &
            (Flux_i_jplushalf-Flux_i)/(1/2.*delta_y_jplushalf)
       del_Flux_minus_e = del_Flux_minus_e + &
            (Flux_e_jplushalf-Flux_e)/(1/2.*delta_y_jplushalf)
    END IF
    
    ! Evaluate Source Term
    E_plus = Sqrt(Ex_plus**2+Ey_plus**2)
    Alpha_minus = A*p*EXP(-B*p/abs(E_plus))
    IF (Alpha_minus > 3.0d5 * L_0) Alpha_minus = 3.0d5 * L_0
    Flux_e = Sqrt(((Flux_e_iplushalf+Flux_e_iminushalf)/2)**2 + &
         ((Flux_e_jplushalf+Flux_e_jminushalf)/2)**2)
    Term_2_minus_S = Alpha_minus*Flux_e- Beta*n_e_minus(inl,jnl)*n_i_minus(inl,jnl)
    
    ! Evaluate Expression
    K_i(stage,inl,jnl) = (-del_Flux_minus_i+Term_2_minus_S)
    K_e(stage,inl,jnl) = (-del_Flux_minus_e+Term_2_minus_S)
END SUBROUTINE FUNCTION_EVALUATION_EXPLICIT


SUBROUTINE SURFACE_CHARGE_EVALUATION(inl,jnl)
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER, INTENT(IN):: inl, jnl
    REAL(8):: Flux_e, Flux_i, Flux_coeff, del_phi_plus
    REAL(8):: delta_y_jplushalf, delta_y_jminushalf
    
    ! Term: -n.Gamma_i
    delta_y_jplushalf = y(jnl+1)-y(jnl)
    delta_y_jminushalf = y(jnl)-y(jnl-1)

    del_phi_plus = phi_plus(inl,jnl+1) - phi_plus(inl,jnl)
    IF (del_phi_plus/delta_y_jplushalf > 0) THEN
       ! Defines alpha in flux BC
       Flux_Coeff = 1.d0
    ELSE
       Flux_Coeff = 0.d0
    END IF
    Flux_i = -1/4.*v_i*n_i_minus(inl,jnl)-Flux_Coeff*mu_i*&
         n_i_minus(inl,jnl)*del_phi_plus/delta_y_jplushalf
    ! Term: +n.Gamma_e
    Flux_e = -1/4.*v_e*n_e_minus(inl,jnl) - gamma*Flux_i + &
         (1.d0-Flux_Coeff)*mu_e*n_e_minus(inl,jnl)*del_phi_plus/delta_y_jplushalf

    sigma_charge(inl,jnl) = sigma_charge(inl,jnl)+delta_time*(-Flux_i+Flux_e)
END SUBROUTINE SURFACE_CHARGE_EVALUATION


SUBROUTINE ION_FLUX(Flux, del_phi, del_h, mu_i, D_i, ni_left, ni_right)
    USE CONTAIN
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: Flux
    REAL(8), INTENT(IN) :: del_phi, del_h, mu_i, D_i, ni_left, ni_right
    REAL(8) :: tolerance
    
    tolerance = 1.d-12
    
    IF (abs(del_phi) < tolerance) THEN
       Flux = D_i*(ni_left-ni_right)/del_h
    ELSE IF (del_phi < 0) THEN
       Flux = -mu_i*(del_phi/del_h)*(ni_left-ni_right*EXP(mu_i*del_phi/D_i))/&
            (1.d0-EXP(mu_i*del_phi/D_i))
    ELSE
       Flux = -mu_i*(del_phi/del_h)*(ni_left*EXP(-mu_i*del_phi/D_i)-ni_right)/&
            (EXP(-mu_i*del_phi/D_i)-1.d0)
    END IF
END SUBROUTINE ION_FLUX


SUBROUTINE ELECTRON_FLUX(Flux, del_phi, del_h, mu_e, D_e, ne_left, ne_right)
    USE CONTAIN
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: Flux
    REAL(8), INTENT(IN) :: del_phi, del_h, mu_e, D_e, ne_left, ne_right
    REAL(8) :: tolerance

    tolerance = 1.d-12

    IF (abs(del_phi) < tolerance) THEN
       Flux = D_e*(ne_left-ne_right)/del_h
    ELSE IF (del_phi > 0) THEN
       Flux = mu_e*(del_phi/del_h)*(ne_left-ne_right*EXP(-mu_e*del_phi/D_e))/&
            (1.d0-EXP(-mu_e*del_phi/D_e))
    ELSE
       Flux = mu_e*(del_phi/del_h)*(ne_left*EXP(mu_e*del_phi/D_e)-ne_right)/&
            (EXP(mu_e*del_phi/D_e)-1.d0)
    END IF
END SUBROUTINE ELECTRON_FLUX
