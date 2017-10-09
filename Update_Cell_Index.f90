SUBROUTINE UPDATE_CELL_INDEX
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::i,j


    ! note: if (x,y) == (0,0), then just set index to -1000 or something (or 0, since indices start with 1 here)
    IF (N_simulated > 0) THEN
        i_cell_vec_prev = i_cell_vec


        IF (use_homogenous_grid .EQV. .true.) THEN

            i_cell_vec(1:N_simulated,1) = FLOOR( (x_vec(1:N_simulated,1)-xmin)/(xmax-xmin)*nx ) + 1

            

            IF (ny>1) THEN
                i_cell_vec(1:N_simulated,2) = FLOOR( (x_vec(1:N_simulated,2)-ymin)/(ymax-ymin)*ny ) + 1
            ELSE
                i_cell_vec(1:N_simulated,2) = 1
            END IF

            ! DO i=1,N_simulated
            !     IF (i_cell_vec(i,2) > ny) THEN
            !         WRITE(*,*) "x_vec= ",x_vec(i,:)
            !         ! WRITE(*,*) "val1 = ",(x_vec(i,1)-xmin)/(xmax-xmin)
            !         WRITE(*,*) "val1 = ",(x_vec(i,2)-ymin)/(ymax-ymin)
            !         WRITE(*,*) "val2 = ",FLOOR( (x_vec(i,2)-ymin)/(ymax-ymin)*ny )
            !         WRITE(*,*) "N_simulated=", N_simulated

            !         WRITE(*,*) "shape(Npc_slice)=", SHAPE(Npc_slice)
            !         WRITE(*,*) "shape(Npc_added)=", SHAPE(Npc_added)
            !         ! WRITE(*,*) "max(i_cell_vec)=", MAXVAL(i_cell_vec(1:N_simulated,1))
            !     ENDIF
            ! END DO

        ELSE
            alpha_x  = -LOG(1/dx_factor)/xmax
            n_inf = 1/(dx_0*alpha_x)
            i_cell_vec(1:N_simulated,1) = FLOOR(n_inf*(1-EXP( -alpha_x*x_vec(1:N_simulated,1) )))+1

            IF (ny > 1) THEN
                alpha_y = -LOG(1/dy_factor) / (ymax-ymid)
                n_inf = 1/(dy_0*alpha_y)
                IF (MOD(ny,2)==0) THEN
                    pos_offset = ny/2. - 0
                    neg_offset = ny/2. - 1
                ELSE
                    pos_offset = ny/2.-.5
                    neg_offset = ny/2.-.5
                ENDIF

                DO i = 1,N_simulated
                    IF (x_vec(i,2) < ymid) THEN
                        i_cell_vec(i,2) = &
                        FLOOR( neg_offset - FLOOR( n_inf*(1-EXP(-alpha_y*(ymid - x_vec(i,2))))) ) 
                    ELSE 
                        i_cell_vec(i,2) = &
                        CEILING( pos_offset + FLOOR( n_inf*(1-EXP(-alpha_y*(x_vec(i,2) - ymid)))) ) 
                    END IF
                END DO

                i_cell_vec(1:N_simulated,2) = i_cell_vec(1:N_simulated,2) + 1
            ELSE
                i_cell_vec(1:N_simulated,2) = 1
            ENDIF

        END IF

        DO i = 1,N_simulated
            IF (i_cell_vec(i,1) > nx) THEN
                i_cell_vec(i,1) = nx
            END IF
            IF (i_cell_vec(i,2) > ny) THEN
                i_cell_vec(i,2) = ny
            END IF
        END DO


        IF (finding_wall_cells .EQV. .false.) THEN

            DO i = 1,N_simulated
                IF (removed_from_sim(i) .eqv. .true.) THEN
                    i_cell_vec(i,1) = 0
                    i_cell_vec(i,2) = 0
                END IF 
            END DO

            CALL SORT_ARRAYS


        END IF

    ENDIF


END SUBROUTINE UPDATE_CELL_INDEX



SUBROUTINE SORT_ARRAYS
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::i,j,i_sorted,current_sum

    ! sort array that stores which particles are in each cell
    CALL CPU_TIME(t0_test)

    ! bucket method based on NASA paper (https://www.nas.nasa.gov/assets/pdf/techreports/1990/rnr-90-017.pdf)
    ! 1. one  pass - count Npc_slice, index spots in sorted cell array  (starting indices)
    ! 2. one pass through cells - find starting index for each cell
    ! 3. one pass - assign each particle to it's sorted place
    
    Npc_slice(:,:) = 0
    Npc_added(:,:) = 0
    
    ! count how many particles in each cell
    DO i = 1,N_simulated
        cx = i_cell_vec(i,1)
        cy = i_cell_vec(i,2)
        IF (cx > 0) THEN
            Npc_slice(cx,cy) = Npc_slice(cx,cy) + 1
        END IF
    END DO

    ! find starting index for each cell in sorted arrays
    current_sum = 0
    DO cx = 1,nx
        DO cy = 1,ny
            starting_index(cx,cy) = current_sum+1
            current_sum = current_sum + Npc_slice(cx,cy)
        END DO
    END DO

    ! generate sorted arrays
    x_vec_unsorted(1:N_simulated,:) = x_vec(1:N_simulated,:)
    v_vec_unsorted(1:N_simulated,:) = v_vec(1:N_simulated,:)
    i_cell_vec_unsorted(1:N_simulated,:) = i_cell_vec(1:N_simulated,:)
    ! removed_from_sim_unsorted = removed_from_sim
    DO i = 1,N_simulated
        cx = i_cell_vec_unsorted(i,1)
        cy = i_cell_vec_unsorted(i,2)
        IF ( (cx > 0) .and. (cy > 0) ) THEN
            i_sorted = starting_index(cx,cy) + Npc_added(cx,cy)
            Npc_added(cx,cy) = Npc_added(cx,cy) + 1

            x_vec(i_sorted,:) = x_vec_unsorted(i,:)
            v_vec(i_sorted,:) = v_vec_unsorted(i,:)
            i_cell_vec(i_sorted,:) = i_cell_vec_unsorted(i,:)
        END IF
    END DO

    removed_from_sim(1:N_simulated) = .false.   
    N_simulated = current_sum
    N_total(ii-1) = N_simulated

    CALL CPU_TIME(t_temp)
    t_test = t_test + (t_temp-t0_test)



END SUBROUTINE SORT_ARRAYS




SUBROUTINE FIND_WALL_CELLS
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::i,j


    ! ! find cell range for each wall
    N_simulated = num_walls*2
    ALLOCATE(x_vec(N_simulated,ndim))
    ALLOCATE(i_cell_vec(N_simulated,2))
    ALLOCATE(i_cell_vec_prev(N_simulated,2))

    DO i = 1,num_walls
        xw1 = x_walls(1,i)
        yw1 = x_walls(2,i)
        xw2 = x_walls(3,i)
        yw2 = x_walls(4,i)
        x_vec(2*i-1,1) = xw1
        x_vec(2*i-1,2) = yw1
        x_vec(2*i,1) = xw2
        x_vec(2*i,2) = yw2
    END DO


    finding_wall_cells = .true.
    CALL UPDATE_CELL_INDEX
    finding_wall_cells = .false.

    DO i = 1,num_walls
        i_cell_lim(i,1) = MINVAL( i_cell_vec( (2*i-1):(2*i) , 1 ) ,1 ) - 2
        i_cell_lim(i,2) = MAXVAL( i_cell_vec( (2*i-1):(2*i) , 1 ) ,1 ) + 2
    END DO
    DO i = 1,num_walls
        IF (i_cell_lim(i,1) < 1) THEN
            i_cell_lim(i,1) = 1
        END IF
        IF (i_cell_lim(i,2) > nx) THEN
            i_cell_lim(i,2) = nx
        END IF
    END DO

    ! WRITE(*,*) "x_vec="
    ! WRITE(*,*) x_vec(1:N_simulated,:)
    ! WRITE(*,*) "i_cell_vec="
    ! WRITE(*,*) i_cell_vec(1:N_simulated,:)
    ! WRITE(*,*) "i_cell_lim(:,1)="
    ! WRITE(*,*) i_cell_lim(1:num_walls,1)
    ! WRITE(*,*) "i_cell_lim(:,2)="
    ! WRITE(*,*) i_cell_lim(1:num_walls,2)

    DEALLOCATE(x_vec)
    DEALLOCATE(i_cell_vec)
    DEALLOCATE(i_cell_vec_prev)
    N_simulated = 0

    ! i_cell_lim(:,1) = 1
    ! i_cell_lim(:,2) = nx



END SUBROUTINE FIND_WALL_CELLS


