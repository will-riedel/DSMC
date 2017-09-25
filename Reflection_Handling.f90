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
        m_w = (yw2-yw1)/(xw2-xw1)
        b_w = yw1 - m_w*xw1

        xy0(1:Num_r,:) = xr_vec_prev(1:Num_r,:)
        xyt(1:Num_r,:) = xr_vec(1:Num_r,:)
        m(1:Num_r) = vr_vec_prev(1:Num_r,2)/vr_vec_prev(1:Num_r,1)
        b(1:Num_r) = xy0(1:Num_r,2)-m(1:Num_r)*xy0(1:Num_r,1)
        
        IF (xw1 == xw2) THEN
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

        ! where they are parallel, dt_cross *should* go to infinity (because xc goes to infinity)
        ! collision_occured(1:Num_r,i) = ( (x0(1:Num_r)<xc(1:Num_r))  .neqv. (xt(1:Num_r)<xc(1:Num_r)) ) &
        !                          .and. ( (xc(1:Num_r)<xw1)          .neqv. (xc(1:Num_r)<xw2) )
        ! collision_occured(1:Num_r,i) = ( (xy0(1:Num_r,1)<xc(1:Num_r))   .neqv. (xyt(1:Num_r,1)<xc(1:Num_r)) ) &
        !                          .and. ( (xc(1:Num_r)<xw1)              .neqv. (xc(1:Num_r)<xw2) )
        ! collision_dt(1:Num_r,i) = (xc(1:Num_r)-x0(1:Num_r)) / vr_vec_prev(1:Num_r,1)
        WHERE (collision_occured(1:Num_r,i))
            ! collision_dt(1:Num_r,i) = (xc(1:Num_r)-x0(1:Num_r)) / vr_vec_prev(1:Num_r,1)
            collision_dt(1:Num_r,i) = (xc(1:Num_r)-xy0(1:Num_r,1)) / vr_vec_prev(1:Num_r,1)
        END WHERE

        ! IF (i == 5) THEN
        !     WRITE(*,*) "---"
        !     WRITE(*,*) i,COUNT(collision_occured(1:Num_r,i)), & 
        !                 COUNT( xy0(1:Num_r,1)<xc(1:Num_r) ), COUNT( xyt(1:Num_r,1)<xc(1:Num_r) ) 

        !     ! WRITE(*,*) i,COUNT( xy0(1:Num_r,1)<xc(1:Num_r) ), COUNT( xyt(1:Num_r,1)<xc(1:Num_r) ) 
        !     ! WRITE(*,*) i,xc(1),xy0(1,:)
        !     ! WRITE(*,*) i,m(1),b(1),m_w,b_w
        ! ENDIF

    END DO



    min_collision_dt(1:Num_r) = MINVAL(collision_dt(1:Num_r,:),2)
    WHERE(min_collision_dt(1:Num_r) == 1.d10) 
        min_collision_dt(1:Num_r) = 2.d10 !  collision_dt = min_collision_dt only where a collision actually occurred
    END WHERE



    

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


        ! x0(1:Num_r) = xr_vec_prev(1:Num_r,1)
        ! y0(1:Num_r) = xr_vec_prev(1:Num_r,2)
        ! xt(1:Num_r) = xr_vec(1:Num_r,1)
        ! yt(1:Num_r) = xr_vec(1:Num_r,2)
        ! m(1:Num_r) = vr_vec_prev(1:Num_r,2)/vr_vec_prev(1:Num_r,1)
        ! b(1:Num_r) = y0(1:Num_r)-m(1:Num_r)*x0(1:Num_r)
        xy0(1:Num_r,:) = xr_vec_prev(1:Num_r,:)
        xyt(1:Num_r,:) = xr_vec(1:Num_r,:)
        m(1:Num_r) = vr_vec_prev(1:Num_r,2)/vr_vec_prev(1:Num_r,1)
        b(1:Num_r) = xy0(1:Num_r,2)-m(1:Num_r)*xy0(1:Num_r,1)

        
        ! xc(1:Num_r) = (b_w-b(1:Num_r)) / (m(1:Num_r) - m_w)
        ! yc(1:Num_r) = m_w*xc(1:Num_r) + b_w

        
        IF (xw1==xw2) THEN
            ! vertical boundary
            WHERE (first_collision(1:Num_r) .eqv. .true.) 
                xr_vec_new(1:Num_r,1) = 2*xw1 - xr_vec(1:Num_r,1)
                vr_vec_new(1:Num_r,1) = -vr_vec(1:Num_r,1)
            ELSEWHERE
            END WHERE

        ELSE IF (yw1==yw2) THEN
            ! horizontal boundary
            WHERE (first_collision(1:Num_r) .eqv. .true.) 
                xr_vec_new(1:Num_r,2) = 2*yw1 - xr_vec(1:Num_r,2)
                vr_vec_new(1:Num_r,2) = -vr_vec(1:Num_r,2)
            ELSEWHERE
            END WHERE
        ELSE 
            ! angled boundary

            ! find angle to rotate
            Theta = wall_angle_vec(i)
            ! Theta = 0.
            ! Rotation_mat_neg(:,1) = (/  COS(-Theta), SIN(-Theta) /)
            ! Rotation_mat_neg(:,2) = (/ -SIN(-Theta),  COS(-Theta) /)
            ! Rotation_mat_pos(:,1) = (/  COS(Theta), SIN(Theta) /)
            ! Rotation_mat_pos(:,2) = (/ -SIN(Theta),  COS(Theta) /)
            
            Rotation_mat_neg(:,1) = (/  COS(-Theta), -SIN(-Theta) /)
            Rotation_mat_neg(:,2) = (/ SIN(-Theta),  COS(-Theta) /)
            Rotation_mat_pos(:,1) = (/  COS(Theta), -SIN(Theta) /)
            Rotation_mat_pos(:,2) = (/ SIN(Theta),  COS(Theta) /)

            ! xy0(1:Num_r) = MATMUL(xy0(1:Num_r),Rotation_mat)
            ! xyt(1:Num_r) = MATMUL(xyt(1:Num_r),Rotation_mat)
            ! xy_w = (/ xw2,yw2 /)
            ! xy_w = MATMUL(xy_w,Rotation_mat)
            ! xw2 = xy_w(1)
            ! yw2 = xy_w(2)
            ! vn(1:Num_r,:) = vr_vec_prev(1:Num_r,:)
            ! vn(1:Num_r,1:2) = MATMUL(vn(1:Num_r,1:2),Rotation_mat)

            ! ! horizontal boundary
            ! m(1:Num_r) = vn(1:Num_r,2)/vn(1:Num_r,1)
            ! b(1:Num_r) = xy0(1:Num_r,2) - m(1:Num_r)*xy0(1:Num_r,1)

            ! xr_vec(1:Num_r) = MATMUL(xr_vec(1:Num_r),Rotation_mat)
            ! xr_vec_new(1:Num_r) = MATMUL(xr_vec_new(1:Num_r),Rotation_mat)
            ! vr_vec(1:Num_r,1:2) = MATMUL(vr_vec(1:Num_r,1:2),Rotation_mat)
            ! vr_vec_new(1:Num_r,1:2) = MATMUL(vr_vec_new(1:Num_r,1:2),Rotation_mat)


            ! translate origin to xw1,yw1
            xw1_0 = xw1
            yw1_0 = yw1
            xw2_0 = xw2
            yw2_0 = yw2
            ! x0 = x0-xw1_0
            ! y0 = y0-yw1_0
            ! xt = xt-xw1_0
            ! yt = yt-yw1_0
            ! xy0(1:Num_r,1) = xy0(1:Num_r,1)-xw1_0
            ! xy0(1:Num_r,2) = xy0(1:Num_r,2)-yw1_0
            xyt(1:Num_r,1) = xyt(1:Num_r,1)-xw1_0
            xyt(1:Num_r,2) = xyt(1:Num_r,2)-yw1_0
            ! xr_vec(1:Num_r,1) = xr_vec(1:Num_r,1)-xw1_0
            ! xr_vec(1:Num_r,2) = xr_vec(1:Num_r,2)-yw1_0
            ! xr_vec_new(1:Num_r,1) = xr_vec_new(1:Num_r,1)-xw1_0
            ! xr_vec_new(1:Num_r,2) = xr_vec_new(1:Num_r,2)-yw1_0
            xw2 = xw2-xw1_0
            yw2 = yw2-yw1_0
            xw1 = 0
            yw1 = 0
            xr_vec_new(1:Num_r,1) = xr_vec_new(1:Num_r,1)-xw1_0
            xr_vec_new(1:Num_r,2) = xr_vec_new(1:Num_r,2)-yw1_0

            ! rotate by -T (to bring to horizontal line)
            xr_vec_new(1:Num_r,1:2) = MATMUL(xr_vec_new(1:Num_r,1:2),Rotation_mat_neg)
            xyt(1:Num_r,1:2) = MATMUL(xyt(1:Num_r,1:2),Rotation_mat_neg)
            vr_vec_new(1:Num_r,1:2) = MATMUL(vr_vec_new(1:Num_r,1:2),Rotation_mat_neg)
            
            ! xy_w = (/ xw2,yw2 /)
            ! xy_w = MATMUL(xy_w,Rotation_mat_neg)
            ! xw2 = xy_w(1)
            ! yw2 = xy_w(2)


            ! DO j=1,Num_r
            !     IF (first_collision(j) .eqv. .true.) THEN
            !         WRITE(*,*) "-------------"
            !         WRITE(*,*) xw1,xw2,xw1_0,xw2_0
            !         WRITE(*,*) yw1,yw2,yw1_0,yw2_0
            !         WRITE(*,*) "theta     =",Theta*180./Pi
            !         WRITE(*,*) "xr_vec    =",xr_vec(j,:)
            !         WRITE(*,*) "xyt       =",xyt(j,:)
            !         WRITE(*,*) "xr_vec_new=",xr_vec_new(j,:)
            !     END IF
            ! END DO


            ! specular condition for horizontal boundary
            WHERE (first_collision(1:Num_r) .eqv. .true.) 
                ! xr_vec_new(1:Num_r,2) = 2*yw1 - xr_vec_new(1:Num_r,2) ! needs to be xt
                xr_vec_new(1:Num_r,2) = 2*yw1 - xyt(1:Num_r,2)
                vr_vec_new(1:Num_r,2) = -vr_vec_new(1:Num_r,2)
            ELSEWHERE
            END WHERE


            ! DO j=1,Num_r
            !     IF (first_collision(j) .eqv. .true.) THEN
            !         WRITE(*,*) "---"
            !         WRITE(*,*) "xr_vec    =",xr_vec(j,:)
            !         WRITE(*,*) "xr_vec_new=",xr_vec_new(j,:)
            !     END IF
            ! END DO

            ! rotate back by +T (to bring to horizontal line)
            xr_vec_new(1:Num_r,1:2) = MATMUL(xr_vec_new(1:Num_r,1:2),Rotation_mat_pos)
            vr_vec_new(1:Num_r,1:2) = MATMUL(vr_vec_new(1:Num_r,1:2),Rotation_mat_pos)
            ! xy_w = MATMUL(xy_w,Rotation_mat_pos)
            ! xw2 = xy_w(1)
            ! yw2 = xy_w(2)

            ! translate origin back to 0,0
            ! xw1 = xw1_0
            ! yw1 = yw1_0
            ! xw2 = xw2_0
            ! yw2 = yw2_0
            xr_vec_new(1:Num_r,1) = xr_vec_new(1:Num_r,1)+xw1_0
            xr_vec_new(1:Num_r,2) = xr_vec_new(1:Num_r,2)+yw1_0

            
        END IF
        




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

    END DO
    


    x_vec(1:Num_r,:) = xr_vec_new(1:Num_r,:)
    v_vec(1:Num_r,:) = vr_vec_new(1:Num_r,:)
    


    ! N_removed = 0
    ! remove(ignore) exiting particles from simulation
    IF (close_outlet .EQV. .false.) THEN
        ! removed_from_sim(1:Num_r) = (removed_from_sim(1:Num_r) .or. reflected_out(1:Num_r))
        WHERE (reflected_out(1:Num_r))
            removed_from_sim(1:Num_r) = .true.
        ENDWHERE
        ! N_removed = N_removed + COUNT(reflected_out(1:Num_r))
    END IF
    IF (close_inlet .EQV. .false.) THEN
        ! removed_from_sim(1:Num_r) = (removed_from_sim(1:Num_r) .or. reflected_in(1:Num_r))
        WHERE (reflected_in(1:Num_r))
            removed_from_sim(1:Num_r) = .true.
        ENDWHERE
        ! N_removed = N_removed + COUNT(reflected_in(1:Num_r))
    END IF
    ! N_simulated = N_simulated - N_removed

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
        WHERE (collision_occured(1:Num_r,i) .eqv. .true.) 
            xs_vec(1:Num_r,2) = 2*yw1 - xs_vec(1:Num_r,2)
            vs_vec(1:Num_r,2) = -vs_vec(1:Num_r,2)
        ELSEWHERE
        END WHERE
    END DO

    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    

END SUBROUTINE SPECULAR_REFLECTION_SOURCE




