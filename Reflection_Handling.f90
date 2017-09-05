SUBROUTINE SPECULAR_REFLECTION(string_in,Num_r,counter)
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    ! REAL(8),DIMENSION(Num_r,ndim):: xr_vec,xr_vec_prev, x_coll, xr_vec_new
    ! REAL(8),DIMENSION(Num_r,3):: vr_vec,vr_vec_prev, vr_vec_new
    ! REAL(8),DIMENSION(num_walls,4)::xr_walls
    ! LOGICAL,DIMENSION(Num_r,num_walls):: collision_occured
    ! REAL(8),DIMENSION(Num_r,num_walls):: collision_dt
    ! REAL(8),DIMENSION(Num_r):: x0,y0,xt,yt,m,b,xc,yc, dt_cross
    ! LOGICAL,DIMENSION(Num_r):: reflected_in, reflected_out, first_collision
    ! INTEGER,DIMENSION(Num_r):: i_cross, i_first
    ! REAL(8):: xw1,xw2,yw1,yw2,xmin,xmax, temp
    ! INTEGER:: Num_r,num_walls,counter, i

    CHARACTER(16):: string_in
    REAL(8),DIMENSION(Num_r,ndim):: xr_vec,xr_vec_prev, x_coll, xr_vec_new
    REAL(8),DIMENSION(Num_r,3):: vr_vec,vr_vec_prev, vr_vec_new
    REAL(8),DIMENSION(4,num_walls)::xr_walls
    LOGICAL,DIMENSION(Num_r,num_walls):: collision_occured
    REAL(8),DIMENSION(Num_r,num_walls):: collision_dt
    REAL(8),DIMENSION(Num_r):: x0,y0,xt,yt,m,b,xc,yc, dt_cross
    LOGICAL,DIMENSION(Num_r):: first_collision,crossed
    INTEGER,DIMENSION(Num_r):: i_cross, i_first
    REAL(8):: xw1,xw2,yw1,yw2, temp
    INTEGER:: Num_r,counter, i,j

    CALL CPU_TIME(t0_BC)

    IF (string_in == 'SOURCE_PARTICLES') THEN
        xr_vec = xs_vec
        xr_vec_prev = xs_vec_prev
        vr_vec = vs_vec
        vr_vec_prev = vs_vec_prev
        ! xr_walls = x_walls(5:)
    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        xr_vec = x_vec
        xr_vec_prev = x_vec_prev
        vr_vec = v_vec
        vr_vec_prev = v_vec_prev
        xr_walls = x_walls
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF

    xr_vec_new = xr_vec
    vr_vec_new = vr_vec
    ! i_refl_out/in?
    collision_occured(:,:) = .false.
    reflected_in(:) = .false.
    reflected_out(:) = .false.

    collision_dt = collision_dt + 1e10

    DO i = 1,num_walls
        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)

        x0 = xr_vec_prev(:,1)
        y0 = xr_vec_prev(:,2)
        xt = xr_vec(:,1)
        yt = xr_vec(:,2)
        m = vr_vec_prev(:,2)/vr_vec_prev(:,1)
        b = y0-m*x0

        ! if (xw1=xw2)
        ! vertical boundary
        IF (yw2 > yw1) THEN
            temp = yw1
            yw1 = yw2
            yw2 = temp
        END IF

        yc = m*xw1+b
        xc(:) = xw1

        ! i_cross = INT( ((x0<xw1) .and. (xw1<xt)) .or. ((x0>xw1).and.(xw1>xt)) )
        crossed = ((x0<xw1).and.(xw1<xt)) .or. ((x0>xw1).and.(xw1>xt))
        dt_cross = (xc-x0)/vr_vec_prev(:,1)
        !crossed = (ABS(dt_cross) < dt) .and. (dt_cross > 0)
        
        ! DO j = 1,N_simulated
        !     IF (xr_vec(j,1) > .61) THEN
        !         WRITE(*,*) " --- (within wall-loop) --"
        !         WRITE(*,*) "x_wall=",xr_walls(i,:)
        !         WRITE(*,*) "x0,xw,xt=",x0(j),xw1,xt(j)
        !         WRITE(*,*) "x = ",xr_vec(j,1)
        !         WRITE(*,*) "x_prev = ",xr_vec_prev(j,1)
        !         WRITE(*,*) "collision_occured=",collision_occured(j,:)
        !         WRITE(*,*) "collision_dt=",collision_dt(j,:)
        !     ENDIF
        ! ENDDO

        ! collision_occured(i_cross,i) = .true.
        ! collision_dt(i_cross,i) = dt_cross(i_cross)
        ! WRITE(*,*) "shape(crossed) = ", SHAPE(crossed)
        ! WRITE(*,*) "shape(coll_occured) = ", SHAPE(collision_occured)
        ! WRITE(*,*) "shape(colision_dt) = ", SHAPE(collision_dt)
        ! WRITE(*,*) "shape(dt_cross) = ", SHAPE(dt_cross)
        collision_occured(:,i) = crossed
        WHERE (crossed)
            ! collision_occured(:,i) = .true.
            collision_dt(:,i) = dt_cross
        ELSEWHERE
        END WHERE

    END DO

    DO i = 1,num_walls
        xw1 = xr_walls(1,i)
        yw1 = xr_walls(2,i)
        xw2 = xr_walls(3,i)
        yw2 = xr_walls(4,i)

        x0 = xr_vec_prev(:,1)
        y0 = xr_vec_prev(:,2)
        xt = xr_vec(:,1)
        yt = xr_vec(:,2)
        m = vr_vec_prev(:,2)/vr_vec_prev(:,1)
        b = y0-m*x0

        ! if (xw1=xw2)
        ! vertical boundary
        IF (yw2 > yw1) THEN
            temp = yw1
            yw1 = yw2
            yw2 = temp
        END IF

        yc = m*xw1+b
        xc(:) = xw1

        ! CALL LOG2INT( i_first, (collision_occured(:,i) .eqv. .true.) .and. (collision_dt(:,i) == MINVAL(collision_dt,2)) , Num_r)
        ! xr_vec_new(i_first,1) = 2*xw1 - xr_vec(i_first,1)
        ! vr_vec_new(i_first,1) = -vr_vec(i_first,1)
        ! x_coll(i_first,1) = xc(i_first)
        ! x_coll(i_first,2) = yc(i_first)

        ! IF ((xw1==xmin).and.(xw2==xmin)) THEN
        !     reflected_in(i_first) = .true.
        ! ELSE IF (xw1==xmax).and.(xw2==xmax) THEN
        !     reflected_out(i_first) = .true.
        ! END IF


        first_collision = (collision_occured(:,i) .eqv. .true.) .and. (collision_dt(:,i) == MINVAL(collision_dt,2))
        WHERE (first_collision) 
            xr_vec_new(:,1) = 2*xw1 - xr_vec(:,1)
            vr_vec_new(:,1) = -vr_vec(:,1)
            x_coll(:,1) = xc
            x_coll(:,2) = yc
        ELSEWHERE
        END WHERE

        IF (string_in == 'SIM_PARTICLES___') THEN
            IF ((xw1==xmin).and.(xw2==xmin)) THEN
                WHERE(first_collision)
                    reflected_in = .true.
                ELSEWHERE
                END WHERE
            ELSE IF ((xw1==xmax).and.(xw2==xmax)) THEN
                WHERE(first_collision)
                    reflected_out = .true.
                ELSEWHERE
                END WHERE
            END IF
        END IF

        ! other angles/wall-types go here



    END DO

    ! xr_vec = xr_vec_new
    ! vr_vec = vr_vec_new

    ! DO i = 1,N_simulated
    !     IF (xr_vec(i,1) > .61) THEN
    !         WRITE(*,*) " --- (after wall-loop) --"
    !         WRITE(*,*) "x = ",xr_vec(i,1)
    !         WRITE(*,*) "x_prev = ",xr_vec_prev(i,1)
    !         WRITE(*,*) "collision_occured=",collision_occured(i,:)
    !         WRITE(*,*) "collision_dt=",collision_dt(i,:)
    !     ENDIF
    ! ENDDO


    IF (string_in == 'SOURCE_PARTICLES') THEN
        xs_vec = xr_vec_new
        vs_vec = vr_vec_new
    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        x_vec = xr_vec_new
        v_vec = vr_vec_new
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF



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





