SUBROUTINE SPECULAR_REFLECTION
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    REAL(8):: temp
    INTEGER:: i,j


    CALL CPU_TIME(t0_BC)


    IF (string_in == 'SOURCE_PARTICLES') THEN
        xr_vec(1:Num_r,:) = xs_vec(1:Num_r,:)
        xr_vec_prev(1:Num_r,:) = xs_vec_prev(1:Num_r,:)
        vr_vec(1:Num_r,:) = vs_vec(1:Num_r,:)
        vr_vec_prev(1:Num_r,:) = vs_vec_prev(1:Num_r,:)
        ! xr_walls = x_walls(5:)
    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        xr_vec(1:Num_r,:) = x_vec(1:Num_r,:)
        xr_vec_prev(1:Num_r,:) = x_vec_prev(1:Num_r,:)
        vr_vec(1:Num_r,:) = v_vec(1:Num_r,:)
        vr_vec_prev(1:Num_r,:) = v_vec_prev(1:Num_r,:)
        xr_walls = x_walls
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF

    xr_vec_new(1:Num_r,:) = xr_vec(1:Num_r,:)
    vr_vec_new(1:Num_r,:) = vr_vec(1:Num_r,:)
    ! i_refl_out/in?
    ! collision_occured(:,:) = .false.
    ! reflected_in(:) = .false.
    ! reflected_out(:) = .false.
    ! collision_dt(:,:) = 1.d10
    collision_occured(1:Num_r,:) = .false.
    reflected_in(1:Num_r) = .false.
    reflected_out(1:Num_r) = .false.
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

        !i_cross =  np.logical_and( (yw1<yc)&(yc<yw2) , np.logical_or( (x0<xw1)&(xw1<xt) , (x0>xw1)&(xw1>xt) ) )
        ! crossed(1:Num_r) = ((x0(1:Num_r)<xw1).and.(xw1<xt(1:Num_r))) .or. ((x0(1:Num_r)>xw1).and.(xw1>xt(1:Num_r)))
        collision_occured(1:Num_r,i) = ((x0(1:Num_r)<xw1).and.(xw1<xt(1:Num_r))) .or. ((x0(1:Num_r)>xw1).and.(xw1>xt(1:Num_r)))
        dt_cross(1:Num_r) = (xc(1:Num_r)-x0(1:Num_r)) / vr_vec_prev(1:Num_r,1)
        

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
                                    (collision_dt(1:Num_r,i) == min_collision_dt(1:Num_r,:))


        ! WHERE (first_collision) 
        WHERE (first_collision(1:Num_r) .eqv. .true.) 
            xr_vec_new(1:Num_r,1) = 2*xw1 - xr_vec(1:Num_r,1)
            vr_vec_new(1:Num_r,1) = -vr_vec(1:Num_r,1)
            x_coll(1:Num_r,1) = xc(1:Num_r)
            x_coll(1:Num_r,2) = yc(1:Num_r)
        ELSEWHERE
        END WHERE
        

        IF (string_in == 'SIM_PARTICLES___') THEN
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
        END IF

        ! other angles/wall-types go here



    END DO
    


    IF (string_in == 'SOURCE_PARTICLES') THEN
        xs_vec(1:Num_r,:) = xr_vec_new(1:Num_r,:)
        vs_vec(1:Num_r,:) = vr_vec_new(1:Num_r,:)
    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        x_vec(1:Num_r,:) = xr_vec_new(1:Num_r,:)
        v_vec(1:Num_r,:) = vr_vec_new(1:Num_r,:)
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF



    ! SHOULD ONLY REMOVE IF SIM PARTICLES
    N_removed = 0
    ! remove(ignore) exiting particles from simulation
    IF (close_outlet .EQV. .false.) THEN
        ! WHERE( reflected_out .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim = (removed_from_sim .or. reflected_out)
        N_removed = N_removed + COUNT(reflected_out)
    END IF
    IF (close_inlet .EQV. .false.) THEN
        ! WHERE( reflected_in .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim = (removed_from_sim .or. reflected_in)
        N_removed = N_removed + COUNT(reflected_in)
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


    CALL CPU_TIME(t0_BC)


    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        xs_vec(1:Num_s,:) = x_vec(1:Num_s,:)
        xs_vec_prev(1:Num_s,:) = x_vec_prev(1:Num_s,:)
        vs_vec(1:Num_s,:) = v_vec(1:Num_s,:)
        vs_vec_prev(1:Num_s,:) = v_vec_prev(1:Num_s,:)
        xs_walls = x_walls
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF

    xs_vec_new(1:Num_s,:) = xs_vec(1:Num_s,:)
    vs_vec_new(1:Num_s,:) = vs_vec(1:Num_s,:)
    ! i_refl_out/in?
    ! collision_occured(:,:) = .false.
    ! reflected_in(:) = .false.
    ! reflected_out(:) = .false.
    ! collision_dt(:,:) = 1.d10
    collision_occured(1:Num_s,:) = .false.
    reflected_in(1:Num_s) = .false.
    reflected_out(1:Num_s) = .false.
    collision_dt(1:Num_s,:) = 1.d10

    DO i = 1,num_walls
        xw1 = xs_walls(1,i)
        yw1 = xs_walls(2,i)
        xw2 = xs_walls(3,i)
        yw2 = xs_walls(4,i)

        x0(1:Num_s) = xs_vec_prev(1:Num_s,1)
        y0(1:Num_s) = xs_vec_prev(1:Num_s,2)
        xt(1:Num_s) = xs_vec(1:Num_s,1)
        yt(1:Num_s) = xs_vec(1:Num_s,2)
        m(1:Num_s) = vs_vec_prev(1:Num_s,2)/vs_vec_prev(1:Num_s,1)
        b(1:Num_s) = y0(1:Num_s)-m(1:Num_s)*x0(1:Num_s)

        ! if (xw1=xw2)
        ! vertical boundary
        IF (yw2 > yw1) THEN
            temp = yw1
            yw1 = yw2
            yw2 = temp
        END IF

        yc(1:Num_s) = m(1:Num_s)*xw1+b(1:Num_s)
        xc(1:Num_s) = xw1

        !i_cross =  np.logical_and( (yw1<yc)&(yc<yw2) , np.logical_or( (x0<xw1)&(xw1<xt) , (x0>xw1)&(xw1>xt) ) )
        ! crossed(1:Num_s) = ((x0(1:Num_s)<xw1).and.(xw1<xt(1:Num_s))) .or. ((x0(1:Num_s)>xw1).and.(xw1>xt(1:Num_s)))
        collision_occured(1:Num_s,i) = ((x0(1:Num_s)<xw1).and.(xw1<xt(1:Num_s))) .or. ((x0(1:Num_s)>xw1).and.(xw1>xt(1:Num_s)))
        dt_cross(1:Num_s) = (xc(1:Num_s)-x0(1:Num_s)) / vs_vec_prev(1:Num_s,1)
        

        ! collision_occured(1:Num_s,i) = crossed(1:Num_s)
        ! WHERE (crossed(1:Num_s))
        WHERE (collision_occured(1:Num_s,i))
            collision_dt(1:Num_s,i) = dt_cross(1:Num_s)
        ELSEWHERE
        END WHERE

    END DO

    min_collision_dt(1:Num_s) = MINVAL(collision_dt(1:Num_s,:),2)

    DO i = 1,num_walls
        xw1 = xs_walls(1,i)
        yw1 = xs_walls(2,i)
        xw2 = xs_walls(3,i)
        yw2 = xs_walls(4,i)

        x0(1:Num_s) = xs_vec_prev(1:Num_s,1)
        y0(1:Num_s) = xs_vec_prev(1:Num_s,2)
        xt(1:Num_s) = xs_vec(1:Num_s,1)
        yt(1:Num_s) = xs_vec(1:Num_s,2)
        m(1:Num_s) = vs_vec_prev(1:Num_s,2)/vs_vec_prev(1:Num_s,1)
        b(1:Num_s) = y0(1:Num_s)-m(1:Num_s)*x0(1:Num_s)

        ! if (xw1=xw2)
        ! vertical boundary
        IF (yw2 > yw1) THEN
            temp = yw1
            yw1 = yw2
            yw2 = temp
        END IF

        yc(1:Num_s) = m(1:Num_s)*xw1+b(1:Num_s)
        xc(1:Num_s) = xw1

        first_collision(1:Num_s) =  (collision_occured(1:Num_s,i) .eqv. .true.) .and. &
                                    (collision_dt(1:Num_s,i) == min_collision_dt(1:Num_s,:))


        ! WHERE (first_collision) 
        WHERE (first_collision(1:Num_s) .eqv. .true.) 
            xs_vec_new(1:Num_s,1) = 2*xw1 - xs_vec(1:Num_s,1)
            vs_vec_new(1:Num_s,1) = -vs_vec(1:Num_s,1)
            x_coll(1:Num_s,1) = xc(1:Num_s)
            x_coll(1:Num_s,2) = yc(1:Num_s)
        ELSEWHERE
        END WHERE
        

        IF (string_in == 'SIM_PARTICLES___') THEN
            IF ((xw1==xmin).and.(xw2==xmin)) THEN
                WHERE(first_collision(1:Num_s))
                    reflected_in(1:Num_s) = .true.
                ELSEWHERE
                END WHERE
            ELSE IF ((xw1==xmax).and.(xw2==xmax)) THEN
                WHERE(first_collision(1:Num_s))
                    reflected_out(1:Num_s) = .true.
                ELSEWHERE
                END WHERE
            END IF
        END IF

        ! other angles/wall-types go here



    END DO
    


    IF (string_in == 'SOURCE_PARTICLES') THEN
        xs_vec(1:Num_s,:) = xs_vec_new(1:Num_s,:)
        vs_vec(1:Num_s,:) = vs_vec_new(1:Num_s,:)
    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        x_vec(1:Num_s,:) = xs_vec_new(1:Num_s,:)
        v_vec(1:Num_s,:) = vs_vec_new(1:Num_s,:)
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF



    ! SHOULD ONLY REMOVE IF SIM PARTICLES
    N_removed = 0
    ! remove(ignore) exiting particles from simulation
    IF (close_outlet .EQV. .false.) THEN
        ! WHERE( reflected_out .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim = (removed_from_sim .or. reflected_out)
        N_removed = N_removed + COUNT(reflected_out)
    END IF
    IF (close_inlet .EQV. .false.) THEN
        ! WHERE( reflected_in .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim = (removed_from_sim .or. reflected_in)
        N_removed = N_removed + COUNT(reflected_in)
    END IF
    N_simulated = N_simulated - N_removed


    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    

END SUBROUTINE SPECULAR_REFLECTION_SOURCE




