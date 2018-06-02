SUBROUTINE ROTATE_TO_AXIAL_PLANE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::i,j

    ! procedure from GA Bird chapter 15
    ! calculate new radius, then rotate back to axial plane


    x_vec(1:N_simulated,2) =  SQRT( (x_vec_prev(1:N_simulated,2)+dt*v_vec(1:N_simulated,2))**2   &
                                 + (dt*v_vec(1:N_simulated,3))**2 )

    ! DO i = 1,N_simulated
    !     IF (x_vec(i,2) > ymax) THEN
    !         WRITE(*,*) "y_prev,y,v =",x_vec_prev(i,2),x_vec(i,2),v_vec(i,:)
    !     END IF
    ! END DO

    ! v = (v1(y1+v1*dt) + w1^2*dt)/y
    v_vec(1:N_simulated,2) =  (v_vec_prev(1:N_simulated,2)*(x_vec_prev(1:N_simulated,2)+dt*v_vec(1:N_simulated,2))   &
                                 + (dt*v_vec(1:N_simulated,3)**2) ) / x_vec(1:N_simulated,2)
    ! w = (w1(y1+v1*dt) - v1*w1*dt)/y
    v_vec(1:N_simulated,3) =  (v_vec_prev(1:N_simulated,3)*(x_vec_prev(1:N_simulated,2)+dt*v_vec(1:N_simulated,2))   &
                                 - (dt*v_vec(1:N_simulated,2)*v_vec(1:N_simulated,3)) ) / x_vec(1:N_simulated,2)


END SUBROUTINE ROTATE_TO_AXIAL_PLANE


SUBROUTINE ROTATE_SOURCE_TO_AXIAL_PLANE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::i,j

    ! procedure from GA Bird chapter 15
    ! calculate new radius, then rotate back to axial plane


    xs_vec(1:Num_s,2) =  SQRT( (xs_vec_prev(1:Num_s,2)+dt*vs_vec(1:Num_s,2))**2   &
                                 + (dt*vs_vec(1:Num_s,3))**2 )

    ! v = (v1(y1+v1*dt) + w1^2*dt)/y
    vs_vec(1:Num_s,2) =  (vs_vec_prev(1:Num_s,2)*(xs_vec_prev(1:Num_s,2)+dt*vs_vec(1:Num_s,2))   &
                                 + (dt*vs_vec(1:Num_s,3)**2) ) / xs_vec(1:Num_s,2)
    ! w = (w1(y1+v1*dt) - v1*w1*dt)/y
    vs_vec(1:Num_s,3) =  (vs_vec_prev(1:Num_s,3)*(xs_vec_prev(1:Num_s,2)+dt*vs_vec(1:Num_s,2))   &
                                 - (dt*vs_vec(1:Num_s,2)*vs_vec(1:Num_s,3)) ) / xs_vec(1:Num_s,2)


END SUBROUTINE ROTATE_SOURCE_TO_AXIAL_PLANE






SUBROUTINE UPDATE_WEIGHTS
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::i,j
    REAL(8)::rn,w1,w2



    ! weight_factor_vec_old(1:N_simulated) = weight_factor_vec(1:N_simulated)
    ! weight_factor_vec(1:N_simulated) = 1 + RWF * x_vec(1:N_simulated,2) / ymax

    i_duplicated = N_simulated + 1
    N_duplicated = 0
    N_deleted = 0
    DO i = 1,N_simulated
        ! w1 = weight_factor_vec_old(i)
        ! w2 = weight_factor_vec(i)
        w1 = weight_factor_vec(i)
        w2 = 1 + RWF*x_vec(i,2)/ymax

        weight_factor_vec(i) = w2

        
        IF (w2 > w1) THEN      
            ! particle moving out / increasing in weight factor
            ! chance of destruction
            CALL RANDOM_NUMBER(rn)
            ! IF (rn < w1/w2) THEN
            IF (rn < (w2-w1)/w2) THEN
                ! Destroy particle
                
                removed_from_sim(i) = .true.

                ! WRITE(*,*)"Removing: i,w1,w2, p=",i,w1,w2,(w2-w1)/w2
                ! WRITE(*,*)"y_prev,y_cur=",x_vec_prev(i,2),x_vec(i,2)
                N_deleted = N_deleted + 1
            END IF

        ELSE
            ! particle moving in / decreasing in weight factor
            ! chance of duplication
            CALL RANDOM_NUMBER(rn)
            IF (rn < (w1-w2)/w2) THEN
                ! Duplicate particle


                x_vec(i_duplicated,:) = x_vec(i,:)
                v_vec(i_duplicated,1:2) = v_vec(i,1:2)
                v_vec(i_duplicated,3) = -v_vec(i,3)        ! flip z-component of velocity to add diversity while maintaining symmetry
                ! i_cell_vec(i_duplicated,:) = i_cell_vec(i,:)


                ! weight_factor_vec(i) = w2
                weight_factor_vec(i_duplicated) = w2
                ! weight_factor_vec_prev(i) = w1
                ! weight_factor_vec_prev(i_duplicated) = w1
                

                ! WRITE(*,*) "w1,w2,i_dup,i,x_vec,=",w1,w2,i_duplicated,i,x_vec(i,:)
                

                i_duplicated = i_duplicated + 1
                N_duplicated = N_duplicated + 1
            END IF

        END IF



        ! DO j = 1,N_simulated
        !     IF (x_vec(j,2) > ymax) THEN
        !         WRITE(*,*) "After loop: i,j,i_duplicated=,y = ",i,j,i_duplicated,x_vec(j,2)
        !     END IF
        ! END DO

    END DO
    N_simulated = N_simulated + N_duplicated
    ! N_simulated = N_simulated + N_duplicated - N_deleted


    ! IF (N_duplicated > 0) THEN
    !     WRITE(*,*) "N_duplicated = ",N_duplicated
    ! END IF

    ! WRITE(*,*) "N_deleted = ",N_deleted

    ! WRITE(*,*) "N_deleted,N_duplicated,N_simulated,N_array = ",-N_deleted,N_duplicated,N_simulated,N_array
    


END SUBROUTINE UPDATE_WEIGHTS

