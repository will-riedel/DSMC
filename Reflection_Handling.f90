SUBROUTINE COMPUTE_REFLECTION
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    REAL(8):: temp
    INTEGER:: i,j,i_min,i_max


    CALL CPU_TIME(t0_BC)

    xr_walls = x_walls

    xr_vec(1:Num_r,:) = x_vec(1:Num_r,:)
    xr_vec_prev(1:Num_r,:) = x_vec_prev(1:Num_r,:)
    vr_vec(1:Num_r,:) = v_vec(1:Num_r,:)
    vr_vec_prev(1:Num_r,:) = v_vec_prev(1:Num_r,:)

    xr_vec_new(1:Num_r,:) = xr_vec(1:Num_r,:)
    vr_vec_new(1:Num_r,:) = vr_vec(1:Num_r,:)
    collision_occured(1:Num_r,:) = .false.
    ! collision_occured_any(1:Num_r) = .false.
    ! reflected_in(1:Num_r) = .false.
    ! reflected_out(1:Num_r) = .false.
    ! removed_from_sim(1:Num_r) = .false.
    collision_dt(1:Num_r,:) = 1.d10

    ! Find where collisions occur, and intersection between each point and each wall
    ! CALL CPU_TIME(t0_BC1)
    DO i = 1,num_walls
        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)

        m_w = (yw2-yw1)/(xw2-xw1)
        b_w = yw1 - m_w*xw1

        ! i_min = starting_index(i_cell_lim_x(i,1),1)
        ! i_max = starting_index(i_cell_lim_x(i,2),ny) + Npc_slice(i_cell_lim_x(i,2),ny)
        DO cx = i_cell_lim_x(i,1) , i_cell_lim_x(i,2)
            DO cy = i_cell_lim_y(i,1) , i_cell_lim_y(i,2)
                i_min = starting_index(cx,cy)
                ! i_max = starting_index(cx,cy) + Npc_slice(cx,cy)
                i_max = final_index(cx,cy)



                ! CALL CPU_TIME(t0_BC2)

                xy0(i_min:i_max,:) = xr_vec_prev(i_min:i_max,:)
                xyt(i_min:i_max,:) = xr_vec(i_min:i_max,:)
                m(i_min:i_max) = vr_vec_prev(i_min:i_max,2)/vr_vec_prev(i_min:i_max,1)
                b(i_min:i_max) = xy0(i_min:i_max,2)-m(i_min:i_max)*xy0(i_min:i_max,1)

                IF (xw1 == xw2) THEN    ! (vertical walls = infinite slopes, slightly different calculation of crossing point)
                    xc(i_min:i_max) = xw1
                    yc(i_min:i_max) = m(i_min:i_max)*xw1 + b(i_min:i_max)
                    collision_occured(i_min:i_max,i) = &
                                    ( ((xy0(i_min:i_max,1)-xc(i_min:i_max))*(xyt(i_min:i_max,1)-xc(i_min:i_max))) < 0 ) &
                             .and.  ( ((yc(i_min:i_max)-yw1)*(yc(i_min:i_max)-yw2)) < 0 )
                ELSE
                    xc(i_min:i_max) = (b_w-b(i_min:i_max)) / (m(i_min:i_max) - m_w)
                    yc(i_min:i_max) = m_w*xc(i_min:i_max) + b_w
                    collision_occured(i_min:i_max,i) = &
                                    ! ( ((xy0(i_min:i_max,1)-xc(i_min:i_max))*(xyt(i_min:i_max,1)-xc(i_min:i_max))) < 0 ) &
                                    ( ((xy0(i_min:i_max,2)-yc(i_min:i_max))*(xyt(i_min:i_max,2)-yc(i_min:i_max))) < 0 ) &
                             .and.  ( ((xc(i_min:i_max)-xw1)*(xc(i_min:i_max)-xw2)) < 0 )
                END IF

                ! IF (i == 8) THEN
                !     ! WRITE(*,*) "this horizontal boundary!"
                !     WRITE(*,*) "i=",i
                !     ! WRITE(*,*) "xy1,xy2=",xw1,yw1,xw2,yw2
                !     ! WRITE(*,*) (yw2-yw1),yw1==yw2
                !     ! WRITE(*,*) "i_cell_lim_x=",i_cell_lim_x(i,:)
                !     ! WRITE(*,*) "i_cell_lim_y=",i_cell_lim_y(i,:)
                !     DO j = i_min,i_max
                !         WRITE(*,*)"---"
                !         ! WRITE(*,*) "xy0,xyt=",xy0(j,:),xyt(j,:)
                !         WRITE(*,*) "y0,yt,yc=",xy0(j,2),xyt(j,2),yc(j)
                !         WRITE(*,*) "y_particle, iy =",xr_vec(j,2),i_cell_vec(j,2)
                !         WRITE(*,*) "yw1, i_cell_lim_y=",yw1,i_cell_lim_y(i,:)
                !         WRITE(*,*) "y_cells_vec=",y_cells_vec

                !         WRITE(*,*) "crossed-y ,on_target-x,collision_occured=", &
                !                                                                 (xy0(j,2)-yc(j)), &
                !                                                                 (xyt(j,2)-yc(j)), &
                !                                                                 ( ( (xy0(j,2)-yc(j))*(xyt(j,2)-yc(j)) ) < 0 ) ,&
                !                                                                 ( ( (xc(j)-xw1)*(xc(j)-xw2) ) < 0 ),&
                !                                                                 collision_occured(j,i)
                !         ! WRITE(*,*) "xyc,c-O=",xc(j),yc(j),collision_occured(j,i)
                !     END DO
                ! END IF

                DO j = i_min,i_max
                    IF (collision_occured(j,i) .eqv. .true.) THEN
                        collision_dt(j,i) = (xc(j)-xy0(j,1)) / vr_vec_prev(j,1)
                        ! collision_occured_any(j) = .true.
                    END IF
                END DO

            END DO
        END DO

    END DO
    ! CALL CPU_TIME(t_temp)
    ! t_BC1 = t_BC1 + (t_temp-t0_BC1)

    
    ! Process reflection/scattering
    ! CALL CPU_TIME(t0_BC4)
    CALL RANDOM_NUMBER(rn_vec)
    DO i = 1,num_walls

        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)

        m_w = (yw2-yw1)/(xw2-xw1)
        b_w = yw1 - m_w*xw1

        IF ((xw1 /= xw2) .and. (yw1 /= yw2) ) THEN
            ! set up extra constants for the angled wall
            ! find angle to rotate
            Theta = wall_angle_vec(i)            
            Rotation_mat_neg(:,1) = (/  COS(-Theta), -SIN(-Theta) /)
            Rotation_mat_neg(:,2) = (/ SIN(-Theta),  COS(-Theta) /)
            Rotation_mat_pos(:,1) = (/  COS(Theta), -SIN(Theta) /)
            Rotation_mat_pos(:,2) = (/ SIN(Theta),  COS(Theta) /)

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



        ! i_min = starting_index(i_cell_lim_x(i,1),1)
        ! i_max = starting_index(i_cell_lim_x(i,2),ny) + Npc_slice(i_cell_lim_x(i,2),ny)
        DO cx = i_cell_lim_x(i,1) , i_cell_lim_x(i,2)
            DO cy = i_cell_lim_y(i,1) , i_cell_lim_y(i,2)
                i_min = starting_index(cx,cy)
                i_max = starting_index(cx,cy) + Npc_slice(cx,cy) - 1    

                DO j = i_min,i_max
                    IF (collision_occured(j,i) .eqv. .true.) THEN
                        IF (collision_dt(j,i) == MINVAL(collision_dt(j,:))) THEN

                            
                            ! WRITE(*,*) "got here"
                            ! WRITE(*,*) "CI=",close_inlet
                            ! IF  (close_inlet .EQV. .false.) THEN
                            !     WRITE(*,*) "open inlet"
                            ! END IF

                            ! remove(ignore) exiting particles from simulation
                            IF ( (close_outlet .EQV. .false.) .and. ((xw1==xmax).and.(xw2==xmax)) ) THEN
                                    removed_from_sim(j) = .true.
                                ! END IF
                            ! ! ELSE IF ( (close_inlet .EQV. .false.) .and. (x_vec(j,2) < ((-1)*x_vec(j,1) + b_source_barrier)) ) THEN
                            !         removed_from_sim(j) = .true.
                            !     ! END IF
                            ! ELSE IF ( (close_inlet .EQV. .false.) .and. (x_vec(j,2) < ((-1)*x_vec(j,1) + b_source_barrier)) ) THEN
                            
                            ELSE IF ( (close_inlet .EQV. .false.) .and. ((xw1==xmin).and.(xw2==xmin)) ) THEN
                            ! ELSE IF ( (close_inlet .EQV. .false.) .and. (xr_vec(j,2) < ((-1)*xr_vec(j,1) + b_source_barrier)) ) THEN
                                    removed_from_sim(j) = .true.
                                ! END IF
                            ELSE


                            ! Otherwise, process collision
                                xy0(j,:) = xr_vec_prev(j,:)
                                xyt(j,:) = xr_vec(j,:)
                                m(j) = vr_vec_prev(j,2)/vr_vec_prev(j,1)
                                b(j) = xy0(j,2)-m(j)*xy0(j,1)


                                IF (xw1 == xw2) THEN
                                    ! vertical boundary
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

                                    ! WRITE(*,*) "horizontal boundary!"
                                    


                                    IF (rn_vec(j) > accommodation) THEN
                                        ! specular reflection
                                        N_specular = N_specular+1

                                        xr_vec_new(j,2) = 2*yw1 - xyt(j,2)
                                        vr_vec_new(j,2) = -vr_vec(j,2)
                                    ELSE
                                        ! diffuse reflection
                                        N_diffuse = N_diffuse+1

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
                                    xr_vec_new(j,1) = xr_vec_new(j,1)-xw1_0
                                    xr_vec_new(j,2) = xr_vec_new(j,2)-yw1_0
                                    xyt(j,1) = xyt(j,1)-xw1_0
                                    xyt(j,2) = xyt(j,2)-yw1_0
                                    xy0(j,1) = xy0(j,1)-xw1_0
                                    xy0(j,2) = xy0(j,2)-yw1_0
                                    ! xc,yc only calculated if diffuse reflection

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

                        END IF
                    END IF
                END DO
 

            END DO
        END DO


    END DO
    ! CALL CPU_TIME(t_temp)
    ! t_BC4 = t_BC4 + (t_temp-t0_BC4)


    IF (initial_distribution(1:5) == "SPLIT") THEN
        ! x_vec is wrong here, it's actually after movement but before reflection
        flux_downstream_total(ii) = COUNT( (xr_vec_prev(1:Num_r,1) < x_split) .and. (xr_vec_new(1:Num_r,1) > x_split) )
        flux_upstream_total(ii) = COUNT( (xr_vec_prev(1:Num_r,1) > x_split) .and. (xr_vec_new(1:Num_r,1) < x_split) )
    END IF

    x_vec(1:Num_r,:) = xr_vec_new(1:Num_r,:)
    v_vec(1:Num_r,:) = vr_vec_new(1:Num_r,:)
    


    ! ! remove(ignore) exiting particles from simulation

    ! IF (close_outlet .EQV. .false.) THEN
    !     DO j = 1,Num_r
    !             IF (reflected_out(j) .eqv. .true.) THEN
    !                 removed_from_sim(j) = .true.
    !             END IF
    !         END DO
    ! END IF
    ! IF (close_inlet .EQV. .false.) THEN
    !     DO j = 1,Num_r
    !             IF (reflected_in(j) .eqv. .true.) THEN
    !                 removed_from_sim(j) = .true.
    !             END IF
    !         END DO
    ! END IF


    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    

END SUBROUTINE COMPUTE_REFLECTION




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

    ! removed_from_sim(1:Num_r) = .false.

    DO j = 1,Num_r

        ! IF (close_outlet .EQV. .false.) THEN
        !     IF (x_vec(j,1)>xmax) THEN
        !         removed_from_sim(j) = .true.
        !     END IF
        ! END IF
        ! IF (close_inlet .EQV. .false.) THEN
        !     IF (x_vec(j,1)<xmin) THEN
        !         removed_from_sim(j) = .true.
        !     END IF
        ! END IF


        IF (x_vec(j,1) > xmax) THEN
            IF (close_outlet .EQV. .true.) THEN
                x_vec(j,1) = 2*xmax - x_vec(j,1)
                v_vec(j,1) = -v_vec(j,1)
            ELSE
                removed_from_sim(j) = .true.
            END IF

        ELSE IF (x_vec(j,1) < xmin) THEN
            IF (close_inlet .EQV. .true.) THEN
                x_vec(j,1) = 2*xmin - x_vec(j,1)
                v_vec(j,1) = -v_vec(j,1)
            ELSE
                removed_from_sim(j) = .true.
            END IF

        END IF


    END DO


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



! SUBROUTINE SPECULAR_REFLECTION_SOURCE_TEMP
!     USE CONTAIN
!     USE PROPERTIES
!     IMPLICIT NONE

!     REAL(8):: temp,xA0,yA0
!     INTEGER:: i,j

!     ! This just reflects any source particles that pass above or below the inlet edges before entering the sim

!     CALL CPU_TIME(t0_BC)

!     Num_r = Num_s

!     ! collision_occured(1:Num_r,1) = (xs_vec(1:Num_r,2)<y_inlet(1))
!     ! collision_occured(1:Num_r,2) = (xs_vec(1:Num_r,2)>y_inlet(2))
!     collision_occured(1:Num_r,1) = ( xs_vec(1:Num_r,2) < ( (1)*xs_vec(1:Num_r,1) + b_source_A) )   ! bottom
!     collision_occured(1:Num_r,2) = ( xs_vec(1:Num_r,2) > ( (1)*xs_vec(1:Num_r,1) + b_source_B) )   ! top


!     Theta = Theta_source
            
!     Rotation_mat_neg(:,1) = (/  COS(-Theta), -SIN(-Theta) /)
!     Rotation_mat_neg(:,2) = (/ SIN(-Theta),  COS(-Theta) /)
!     Rotation_mat_pos(:,1) = (/  COS(Theta), -SIN(Theta) /)
!     Rotation_mat_pos(:,2) = (/ SIN(Theta),  COS(Theta) /)

!     xA0 = x_source_corners(1,1)
!     yA0 = x_source_corners(1,2)

!     ! translate to xy_A
!     xs_vec(1:Num_r,1) = xs_vec(1:Num_r,1) - xA0
!     xs_vec(1:Num_r,2) = xs_vec(1:Num_r,1) - yA0

!     ! rotate by negative theta


!     DO i = 1,2
!         yw1 = y_inlet(i)

!         DO j = 1,Num_r
!             IF (collision_occured(j,i) .eqv. .true.) THEN
!                 xs_vec(j,2) = 2*yw1 - xs_vec(j,2)
!                 vs_vec(j,2) = -vs_vec(j,2)
!             END IF
!         END DO
!     END DO



!     CALL CPU_TIME(t_temp)
!     t_BC = t_BC + (t_temp-t0_BC)
    

! END SUBROUTINE SPECULAR_REFLECTION_SOURCE_TEMP



