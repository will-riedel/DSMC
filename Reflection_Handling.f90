SUBROUTINE SPECULAR_REFLECTION
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    REAL(8):: temp
    INTEGER:: i,j


    CALL CPU_TIME(t0_BC)


    xr_vec(1:Num_r,:) = x_vec(1:Num_r,:)
    xr_vec_prev(1:Num_r,:) = x_vec_prev(1:Num_r,:)
    vr_vec(1:Num_r,:) = v_vec(1:Num_r,:)
    vr_vec_prev(1:Num_r,:) = v_vec_prev(1:Num_r,:)
    xr_walls = x_walls

    xr_vec_new(1:Num_r,:) = xr_vec(1:Num_r,:)
    vr_vec_new(1:Num_r,:) = vr_vec(1:Num_r,:)
    collision_occured(1:Num_r,:) = .false.
    reflected_in(1:Num_r) = .false.
    reflected_out(1:Num_r) = .false.
    removed_from_sim(1:Num_r) = .false.
    collision_dt(1:Num_r,:) = 1.d10

    DO i = 1,num_walls
        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)

        x0(1:Num_r) = xr_vec_prev(1:Num_r,1)
        y0(1:Num_r) = xr_vec_prev(1:Num_r,2)
        xt(1:Num_r) = xr_vec(1:Num_r,1)
        yt(1:Num_r) = xr_vec(1:Num_r,2)
        m(1:Num_r) = vr_vec_prev(1:Num_r,2)/vr_vec_prev(1:Num_r,1)
        b(1:Num_r) = y0(1:Num_r)-m(1:Num_r)*x0(1:Num_r)

        ! if (xw1=xw2)
        ! vertical boundary
        IF (yw2 > yw1) THEN
            temp = yw1
            yw1 = yw2
            yw2 = temp
        END IF

        yc(1:Num_r) = m(1:Num_r)*xw1+b(1:Num_r)
        xc(1:Num_r) = xw1

        ! crossed(1:Num_r) = ((x0(1:Num_r)<xw1).and.(xw1<xt(1:Num_r))) .or. ((x0(1:Num_r)>xw1).and.(xw1>xt(1:Num_r)))
        ! collision_occured(1:Num_r,i) = ((x0(1:Num_r)<xw1).and.(xw1<xt(1:Num_r))) .or. ((x0(1:Num_r)>xw1).and.(xw1>xt(1:Num_r)))
        collision_occured(1:Num_r,i) = ( (x0(1:Num_r)<xw1) .neqv. (xt(1:Num_r)<xw1) )
        dt_cross(1:Num_r) = (xc(1:Num_r)-x0(1:Num_r)) / vr_vec_prev(1:Num_r,1)
        
        ! ########################################################################################
        ! ########################################################################################
        ! ########################################################################################
        ! is this valid? Double-check, might save time
        ! collision_occured(1:Num_r,i) = ( (y0(1:Num_r)<yw1) .neqv. (yt(1:Num_r)<yw1) )
        ! collision_occured(1:Num_r,i) = ( (x0(1:Num_r)<xw1) .neqv. (xt(1:Num_r)<xw1) )
        ! ########################################################################################
        ! ########################################################################################
        ! ########################################################################################



        ! collision_occured(1:Num_r,i) = crossed(1:Num_r)
        ! WHERE (crossed(1:Num_r))
        WHERE (collision_occured(1:Num_r,i))
            collision_dt(1:Num_r,i) = dt_cross(1:Num_r)
        ELSEWHERE
        END WHERE

    END DO

    min_collision_dt(1:Num_r) = MINVAL(collision_dt(1:Num_r,:),2)

    DO i = 1,num_walls
        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)

        x0(1:Num_r) = xr_vec_prev(1:Num_r,1)
        y0(1:Num_r) = xr_vec_prev(1:Num_r,2)
        xt(1:Num_r) = xr_vec(1:Num_r,1)
        yt(1:Num_r) = xr_vec(1:Num_r,2)
        m(1:Num_r) = vr_vec_prev(1:Num_r,2)/vr_vec_prev(1:Num_r,1)
        b(1:Num_r) = y0(1:Num_r)-m(1:Num_r)*x0(1:Num_r)

        ! if (xw1=xw2)
        ! vertical boundary
        IF (yw2 > yw1) THEN
            temp = yw1
            yw1 = yw2
            yw2 = temp
        END IF

        yc(1:Num_r) = m(1:Num_r)*xw1+b(1:Num_r)
        xc(1:Num_r) = xw1

        first_collision(1:Num_r) =  (collision_occured(1:Num_r,i) .eqv. .true.) .and. &
                                    (collision_dt(1:Num_r,i) == min_collision_dt(1:Num_r))


        ! WHERE (first_collision) 
        WHERE (first_collision(1:Num_r) .eqv. .true.) 
            xr_vec_new(1:Num_r,1) = 2*xw1 - xr_vec(1:Num_r,1)
            vr_vec_new(1:Num_r,1) = -vr_vec(1:Num_r,1)
            x_coll(1:Num_r,1) = xc(1:Num_r)
            x_coll(1:Num_r,2) = yc(1:Num_r)
        ELSEWHERE
        END WHERE
        

        IF ((xw1==xmin).and.(xw2==xmin)) THEN
            WHERE(first_collision(1:Num_r))
                reflected_in(1:Num_r) = .true.
            ELSEWHERE
            END WHERE
        ELSE IF ((xw1==xmax).and.(xw2==xmax)) THEN
            WHERE(first_collision(1:Num_r))
                reflected_out(1:Num_r) = .true.
            ELSEWHERE
            END WHERE
        END IF

        ! other angles/wall-types go here



    END DO
    


    x_vec(1:Num_r,:) = xr_vec_new(1:Num_r,:)
    v_vec(1:Num_r,:) = vr_vec_new(1:Num_r,:)
    


    N_removed = 0
    ! remove(ignore) exiting particles from simulation
    IF (close_outlet .EQV. .false.) THEN
        ! WHERE( reflected_out .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim(1:Num_r) = (removed_from_sim(1:Num_r) .or. reflected_out(1:Num_r))
        N_removed = N_removed + COUNT(reflected_out(1:Num_r))
    END IF
    IF (close_inlet .EQV. .false.) THEN
        ! WHERE( reflected_in .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim(1:Num_r) = (removed_from_sim(1:Num_r) .or. reflected_in(1:Num_r))
        N_removed = N_removed + COUNT(reflected_in(1:Num_r))
    END IF
    N_simulated = N_simulated - N_removed

    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    

END SUBROUTINE SPECULAR_REFLECTION




! ########################################################################################
! ########################################################################################
! ########################################################################################
! ########################################################################################
! ########################################################################################
! ########################################################################################
! ########################################################################################
! ########################################################################################





SUBROUTINE SPECULAR_REFLECTION_SOURCE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    REAL(8):: temp
    INTEGER:: i,j

    ! This just reflects any source particles that pass above or below the inlet edges before entering the sim


    CALL CPU_TIME(t0_BC)

    Num_r = Num_s

    collision_occured(1:Num_r,1) = (xs_vec(1:Num_r,2)<y_inlet(1))
    collision_occured(1:Num_r,2) = (xs_vec(1:Num_r,2)>y_inlet(2))


    DO i = 1,2
        yw1 = y_inlet(i)
        WHERE (collision_occured(1:Num_r,i) .eqv. .true.) 
            xs_vec(1:Num_r,2) = 2*yw1 - xs_vec(1:Num_r,2)
            vs_vec(1:Num_r,2) = -vs_vec(1:Num_r,2)
        ELSEWHERE
        END WHERE
    END DO

    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    

END SUBROUTINE SPECULAR_REFLECTION_SOURCE




