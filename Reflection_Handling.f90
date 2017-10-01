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

    ! Find where collisions occur, and intersection between each point and each wall
    DO i = 1,num_walls
        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)
        m_w = (yw2-yw1)/(xw2-xw1)
        b_w = yw1 - m_w*xw1

        xy0(1:Num_r,:) = xr_vec_prev(1:Num_r,:)
        xyt(1:Num_r,:) = xr_vec(1:Num_r,:)
        m(1:Num_r) = vr_vec_prev(1:Num_r,2)/vr_vec_prev(1:Num_r,1)
        b(1:Num_r) = xy0(1:Num_r,2)-m(1:Num_r)*xy0(1:Num_r,1)
        
        IF (xw1 == xw2) THEN    ! (infinite slopes, slightly different calculation of crossing point)
            xc(1:Num_r) = xw1
            yc(1:Num_r) = m(1:Num_r)*xw1 + b(1:Num_r)
            collision_occured(1:Num_r,i) = ( (xy0(1:Num_r,1)<xc(1:Num_r))   .neqv. (xyt(1:Num_r,1)<xc(1:Num_r)) ) &
                                 .and. ( (yc(1:Num_r)<yw1)              .neqv. (yc(1:Num_r)<yw2) )
        ELSE
            xc(1:Num_r) = (b_w-b(1:Num_r)) / (m(1:Num_r) - m_w)
            yc(1:Num_r) = m_w*xc(1:Num_r) + b_w
            collision_occured(1:Num_r,i) = ( (xy0(1:Num_r,1)<xc(1:Num_r))   .neqv. (xyt(1:Num_r,1)<xc(1:Num_r)) ) &
                                 .and. ( (xc(1:Num_r)<xw1)              .neqv. (xc(1:Num_r)<xw2) )
        END IF

        DO j = 1,Num_r
            IF (collision_occured(j,i) .eqv. .true.) THEN
                collision_dt(j,i) = (xc(j)-xy0(j,1)) / vr_vec_prev(j,1)
            END IF
        END DO

    END DO


    min_collision_dt(1:Num_r) = MINVAL(collision_dt(1:Num_r,:),2)
    DO j = 1,Num_r
        IF (min_collision_dt(j) == 1.d10) THEN
            min_collision_dt(j) = 2.d10 !  collision_dt = min_collision_dt only where a collision actually occurred
        END IF
    END DO



    
    ! Process reflection/scattering
    CALL RANDOM_NUMBER(rn_vec)
    DO i = 1,num_walls
        ! first_collision(1:Num_r) =  (collision_occured(1:Num_r,i) .eqv. .true.) .and. &
        !                             (collision_dt(1:Num_r,i) == min_collision_dt(1:Num_r))
        first_collision(1:Num_r) =  (collision_dt(1:Num_r,i) == min_collision_dt(1:Num_r))

        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)
        m_w = (yw2-yw1)/(xw2-xw1)
        b_w = yw1 - m_w*xw1

        xy0(1:Num_r,:) = xr_vec_prev(1:Num_r,:)
        xyt(1:Num_r,:) = xr_vec(1:Num_r,:)
        m(1:Num_r) = vr_vec_prev(1:Num_r,2)/vr_vec_prev(1:Num_r,1)
        b(1:Num_r) = xy0(1:Num_r,2)-m(1:Num_r)*xy0(1:Num_r,1)
        

        IF ((xw1 /= xw2) .and. (yw1 /= yw2) ) THEN
            ! set up extra constants for the angled wall
            ! find angle to rotate
            Theta = wall_angle_vec(i)            
            Rotation_mat_neg(:,1) = (/  COS(-Theta), -SIN(-Theta) /)
            Rotation_mat_neg(:,2) = (/ SIN(-Theta),  COS(-Theta) /)
            Rotation_mat_pos(:,1) = (/  COS(Theta), -SIN(Theta) /)
            Rotation_mat_pos(:,2) = (/ SIN(Theta),  COS(Theta) /)




            ! xy_w(1) = xw2
            ! xy_w(2) = yw2
            ! WRITE(*,*) "Theta=",Theta
            ! WRITE(*,*) "xyw orig =",xy_w
            ! xy_w(1) = xy_w(1) - xw1_0
            ! xy_w(2) = xy_w(2) - yw1_0
            ! WRITE(*,*) "xyw translated =",xy_w
            ! xy_w = MATMUL(xy_w,Rotation_mat_neg)
            ! WRITE(*,*) "xyw_rotated away =",xy_w
            ! xy_w = MATMUL(xy_w,Rotation_mat_pos)
            ! WRITE(*,*) "xyw_rotated back =",xy_w
            ! xy_w(1) = xy_w(1) + xw1_0
            ! xy_w(2) = xy_w(2) + yw1_0
            ! WRITE(*,*) "xyw translated back =",xy_w



            ! translate origin to xw1,yw1
            xw1_0 = xw1
            yw1_0 = yw1
            xw2_0 = xw2
            yw2_0 = yw2
            xw2 = xw2-xw1_0
            yw2 = yw2-yw1_0
            xw1 = 0
            yw1 = 0


        END IF

        ! WRITE(*,*) "walls=",xr_walls(:,1:num_walls)
        ! WRITE(*,*) "xw1,xw2=",xw1,xw2
        ! WRITE(*,*) "yw1,yw2=",yw1,yw2


        DO j = 1,Num_r
            IF (first_collision(j) .eqv. .true.) THEN
                IF (xw1 == xw2) THEN
                    ! vertical boundary
                    ! WRITE(*,*) "vertical boundary: ",xw1,xw2,yw1,yw2
                    IF (rn_vec(j) > accommodation) THEN
                        ! specular reflection
                        N_specular = N_specular+1

                        xr_vec_new(j,1) = 2*xw1 - xyt(j,1)
                        vr_vec_new(j,1) = -vr_vec(j,1)
                    ELSE
                        ! diffuse reflection
                        N_diffuse = N_diffuse+1

                        xc(j) = xw1
                        yc(j) = m(j)*xc(j) + b(j)

                        CALL RANDOM_NUMBER(alpha_vec)
                        IF (vr_vec_new(j,1) > 0) THEN
                            temp = -1
                        ELSE
                            temp = 1
                        END IF
                        vr_vec_new(j,1) = temp*SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(3)))
                        vr_vec_new(j,2) = SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(1)))*COS(2*Pi*alpha_vec(2))
                        vr_vec_new(j,3) = SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(1)))*SIN(2*Pi*alpha_vec(2))
                        xr_vec_new(j,1) = xc(j) + vr_vec_new(j,1)*(dt - collision_dt(j,i))
                        xr_vec_new(j,2) = yc(j) + vr_vec_new(j,2)*(dt - collision_dt(j,i))

                    END IF

                ELSE IF (yw1 == yw2) THEN
                    ! horizontal boundary
                    ! WRITE(*,*) "horizontal boundary: ",xw1,xw2,yw1,yw2
                    IF (rn_vec(j) > accommodation) THEN
                        ! specular reflection
                        N_specular = N_specular+1

                        xr_vec_new(j,2) = 2*yw1 - xyt(j,2)
                        vr_vec_new(j,2) = -vr_vec(j,2)
                    ELSE
                        ! diffuse reflection
                        N_diffuse = N_diffuse+1

                        ! xyc(j,1) = (yw1-b(j))/m(j)
                        ! xyc(j,2) = yw1
                        xc(j) = (yw1-b(j))/m(j)
                        yc(j) = yw1

                        CALL RANDOM_NUMBER(alpha_vec)
                        IF(vr_vec_new(j,2) > 0) THEN
                            temp = -1
                        ELSE
                            temp = 1
                        END IF
                        vr_vec_new(j,1) = SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(1)))*COS(2*Pi*alpha_vec(2))
                        vr_vec_new(j,2) = temp*SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(3)))
                        vr_vec_new(j,3) = SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(1)))*SIN(2*Pi*alpha_vec(2))
                        xr_vec_new(j,1) = xc(j) + vr_vec_new(j,1)*(dt - collision_dt(j,i))
                        xr_vec_new(j,2) = yc(j) + vr_vec_new(j,2)*(dt - collision_dt(j,i))

                    END IF

                ELSE 
                    ! angled boundary
                    ! WRITE(*,*) "angled boundary: ",xw1,xw2,yw1,yw2

                    xr_vec_new(j,1) = xr_vec_new(j,1)-xw1_0
                    xr_vec_new(j,2) = xr_vec_new(j,2)-yw1_0
                    xyt(j,1) = xyt(j,1)-xw1_0
                    xyt(j,2) = xyt(j,2)-yw1_0
                    xy0(j,1) = xy0(j,1)-xw1_0
                    xy0(j,2) = xy0(j,2)-yw1_0
                    ! xc(1:Num_r) = xw1
                    ! yc(1:Num_r) = m(1:Num_r)*xw1 + b(1:Num_r)

                    ! rotate by -T (to bring to horizontal line)
                    xr_vec_new(j,1:2) = MATMUL(xr_vec_new(j,1:2),Rotation_mat_neg)
                    xyt(j,1:2) = MATMUL(xyt(j,1:2),Rotation_mat_neg)
                    xy0(j,1:2) = MATMUL(xy0(j,1:2),Rotation_mat_neg)
                    vr_vec_new(j,1:2) = MATMUL(vr_vec_new(j,1:2),Rotation_mat_neg)

                    IF (rn_vec(j) > accommodation) THEN
                        ! specular reflection (for horizontal boundary)
                        N_specular = N_specular+1

                        xr_vec_new(j,2) = 2*yw1 - xyt(j,2)
                        vr_vec_new(j,2) = -vr_vec_new(j,2)
                    ELSE
                        ! diffuse reflection
                        N_diffuse = N_diffuse+1

                        m(j) = vr_vec_new(j,2)/vr_vec_new(j,1)
                        b(j) = xy0(j,2)-m(j)*xy0(j,1)
                        xc(j) = (yw1-b(j))/m(j)
                        yc(j) = yw1

                        CALL RANDOM_NUMBER(alpha_vec)
                        IF(vr_vec_new(j,2) > 0) THEN
                            temp = -1
                        ELSE
                            temp = 1
                        END IF
                        vr_vec_new(j,1) = SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(1)))*COS(2*Pi*alpha_vec(2))
                        vr_vec_new(j,2) = temp*SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(3)))
                        vr_vec_new(j,3) = SQRT(-2*k_b*T_g/m_g*LOG(alpha_vec(1)))*SIN(2*Pi*alpha_vec(2))
                        xr_vec_new(j,1) = xc(j) + vr_vec_new(j,1)*(dt - collision_dt(j,i))
                        xr_vec_new(j,2) = yc(j) + vr_vec_new(j,2)*(dt - collision_dt(j,i))
                    END IF

                    ! rotate back by +T (to bring to horizontal line)
                    xr_vec_new(j,1:2) = MATMUL(xr_vec_new(j,1:2),Rotation_mat_pos)
                    vr_vec_new(j,1:2) = MATMUL(vr_vec_new(j,1:2),Rotation_mat_pos)

                    ! translate origin back to 0,0
                    xr_vec_new(j,1) = xr_vec_new(j,1)+xw1_0
                    xr_vec_new(j,2) = xr_vec_new(j,2)+yw1_0
                END IF

            END IF
        END DO


        IF ((xw1==xmin).and.(xw2==xmin)) THEN
            DO j = 1,Num_r
                IF (first_collision(j) .eqv. .true.) THEN
                    reflected_in(j) = .true.
                END IF
            END DO
        ELSE IF ((xw1==xmax).and.(xw2==xmax)) THEN
            DO j = 1,Num_r
                IF (first_collision(j) .eqv. .true.) THEN
                    reflected_out(j) = .true.
                END IF
            END DO
        END IF



    END DO
    


    x_vec(1:Num_r,:) = xr_vec_new(1:Num_r,:)
    v_vec(1:Num_r,:) = vr_vec_new(1:Num_r,:)
    


    ! remove(ignore) exiting particles from simulation

    IF (close_outlet .EQV. .false.) THEN
        DO j = 1,Num_r
                IF (reflected_out(j) .eqv. .true.) THEN
                    removed_from_sim(j) = .true.
                END IF
            END DO
    END IF
    IF (close_inlet .EQV. .false.) THEN
        DO j = 1,Num_r
                IF (reflected_in(j) .eqv. .true.) THEN
                    removed_from_sim(j) = .true.
                END IF
            END DO
    END IF



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

SUBROUTINE SPECULAR_REFLECTION_1D
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

        ! IF (xw1=xw2) THEN
        ! vertical boundary
        IF (yw2 > yw1) THEN
            ! temp = yw1
            ! yw1 = yw2
            ! yw2 = temp
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
    


    ! N_removed = 0
    ! remove(ignore) exiting particles from simulation
    IF (close_outlet .EQV. .false.) THEN
        ! WHERE( reflected_out .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim(1:Num_r) = (removed_from_sim(1:Num_r) .or. reflected_out(1:Num_r))
        ! N_removed = N_removed + COUNT(reflected_out(1:Num_r))
    END IF
    IF (close_inlet .EQV. .false.) THEN
        ! WHERE( reflected_in .EQV. .true.)
        !     removed_from_sim = .true.
        ! ELSEWHERE
        ! END WHERE
        removed_from_sim(1:Num_r) = (removed_from_sim(1:Num_r) .or. reflected_in(1:Num_r))
        ! N_removed = N_removed + COUNT(reflected_in(1:Num_r))
    END IF
    ! N_simulated = N_simulated - N_removed

    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    

END SUBROUTINE SPECULAR_REFLECTION_1D




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
        ! WHERE (collision_occured(1:Num_r,i) .eqv. .true.) 
        !     xs_vec(1:Num_r,2) = 2*yw1 - xs_vec(1:Num_r,2)
        !     vs_vec(1:Num_r,2) = -vs_vec(1:Num_r,2)
        ! ELSEWHERE
        ! END WHERE

        DO j = 1,Num_r
            IF (collision_occured(j,i) .eqv. .true.) THEN
                xs_vec(j,2) = 2*yw1 - xs_vec(j,2)
                vs_vec(j,2) = -vs_vec(j,2)
            END IF
        END DO
    END DO

    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    

END SUBROUTINE SPECULAR_REFLECTION_SOURCE




