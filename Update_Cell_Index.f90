SUBROUTINE UPDATE_CELL_INDEX
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::i,cx_test,i_test,i_sorted,cx_cur,current_sum

    ! note: if (x,y) == (0,0), then just set index to -1000 or something (or 0, since indices start with 1 here)
    IF (N_all > 0) THEN

        i_cell_vec_prev = i_cell_vec


        IF (use_homogenous_grid .EQV. .true.) THEN

            i_cell_vec(1:N_all,1) = FLOOR( (x_vec(1:N_all,1)-xmin)/(xmax-xmin)*nx ) + 1

            ! DO i=1,N_all
            !     IF (i_cell_vec(i,1) > nx) THEN
            !         WRITE(*,*) "x_vec= ",x_vec(i,1)
            !         WRITE(*,*) "val1 = ",(x_vec(i,1)-xmin)/(xmax-xmin)
            !         WRITE(*,*) "val2 = ",FLOOR( (x_vec(i,1)-xmin)/(xmax-xmin)*nx )
            !         WRITE(*,*) "N_all=", N_all

            !         WRITE(*,*) "shape(Npc_slice)=", SHAPE(Npc_slice)
            !         WRITE(*,*) "shape(Npc_added)=", SHAPE(Npc_added)
            !         WRITE(*,*) "max(i_cell_vec)=", MAXVAL(i_cell_vec(1:N_all,1))
            !     ENDIF
            ! END DO

            IF (ny>1) THEN
                i_cell_vec(1:N_all,2) = FLOOR( (x_vec(1:N_all,2)-xmin)/(xmax-xmin)*ny ) + 1
            ELSE
                i_cell_vec(1:N_all,2) = 1
            END IF

        ELSE
            alpha_x  = -LOG(1/dx_factor)/xmax
            n_inf = 1/(dx_0*alpha_x)
            i_cell_vec(1:N_all,1) = FLOOR(n_inf*(1-EXP( -alpha_x*x_vec(1:N_all,1) )))+1

            IF (ny > 1) THEN
                alpha_y = -LOG(1/dy_factor) / (ymax-ymid)
                n_inf = 1/(dx_0*alpha_y)
                IF (MOD(ny,2)==0) THEN
                    pos_offset = ny/2. - 0
                    neg_offset = ny/2. - 1
                ELSE
                    pos_offset = ny/2.-.5
                    neg_offset = ny/2.-.5
                ENDIF

                WHERE (x_vec(1:N_all,2) < ymid)
                    i_cell_vec(1:N_all,2) = FLOOR( neg_offset - FLOOR( n_inf*(1-EXP(-alpha_y*(ymid - x_vec(1:N_all,2))))) )
                ELSEWHERE
                    i_cell_vec(1:N_all,2) = CEILING( pos_offset - FLOOR( n_inf*(1-EXP(-alpha_y*(x_vec(1:N_all,2) - ymid)))) )
                END WHERE

                i_cell_vec(1:N_all,2) = i_cell_vec(1:N_all,2) + 1
            ELSE
                i_cell_vec(1:N_all,2) = 1
            ENDIF

        END IF

        WHERE (removed_from_sim .EQV. .true.)
            i_cell_vec(:,1) = 0
            i_cell_vec(:,2) = 0
        ELSEWHERE
        END WHERE



        ! sort array that stores which particles are in each cell
        CALL CPU_TIME(t0_test)

        ! bucket method based on NASA paper (https://www.nas.nasa.gov/assets/pdf/techreports/1990/rnr-90-017.pdf)
        ! 1. one  pass - count Npc_slice, index spots in sorted cell array  (starting indices)
        ! 2. one pass through cells - find starting index for each cell
        ! 3. one pass - assign each particle to it's sorted place


        
        ! WRITE(*,*) "cx_cur=", cx_cur
! problem happens here because there's a cell listed as out of bouncs

        ! find how many particles in each cell
        Npc_slice(:,:) = 0
        Npc_added(:,:) = 0
        ! DO i = 1,N_all
        DO i = 1,N_simulated
            cx_cur = i_cell_vec(i,1)
            IF (cx_cur > 0) THEN
                Npc_slice(cx_cur,cy) = Npc_slice(cx_cur,cy) + 1
            END IF
        END DO


        ! find starting index for each cell in sorted arrays
        current_sum = 0
        DO cx_cur = 1,nx
            starting_index(cx_cur,cy) = current_sum+1
            current_sum = current_sum + Npc_slice(cx_cur,cy)
        END DO

        ! N_exited = N_simulated - current_sum


        ! DO i=1,N_all
        !     IF (i_cell_vec(i,1) > nx) THEN
        !         WRITE(*,*) " --- (updating cell index) --- "
        !         WRITE(*,*) "i = ",i
        !         WRITE(*,*) "rfs(i) = ",removed_from_sim(i)
        !         WRITE(*,*) "x_vec= ",x_vec(i,1)
        !         WRITE(*,*) "x_vec_prev= ",x_vec_prev(i,1)
        !         WRITE(*,*) "val1 = ",(x_vec(i,1)-xmin)/(xmax-xmin)
        !         WRITE(*,*) "val2 = ",FLOOR( (x_vec(i,1)-xmin)/(xmax-xmin)*nx )
        !         WRITE(*,*) "N_all=", N_all

        !         WRITE(*,*) "shape(Npc_slice)=", SHAPE(Npc_slice)
        !         WRITE(*,*) "shape(Npc_added)=", SHAPE(Npc_added)
        !         WRITE(*,*) "max(i_cell_vec)=", MAXVAL(i_cell_vec(1:N_all,1))
        !     ENDIF
        ! END DO

        ! generate sorted arrays
        x_vec_unsorted = x_vec
        v_vec_unsorted = v_vec
        i_cell_vec_unsorted = i_cell_vec
        ! removed_from_sim_unsorted = removed_from_sim
        ! DO i = 1,N_all
        DO i = 1,N_simulated
            cx_cur = i_cell_vec(i,1)
            IF (cx_cur > 0) THEN
                i_sorted = starting_index(cx_cur,cy) + Npc_added(cx_cur,cy)
                Npc_added(cx_cur,cy) = Npc_added(cx_cur,cy) + 1

                x_vec(i_sorted,:) = x_vec_unsorted(i,:)
                v_vec(i_sorted,:) = v_vec_unsorted(i,:)
                i_cell_vec(i_sorted,:) = i_cell_vec_unsorted(i,:)
            END IF

        END DO



        ! N_simulated = N_simulated - N_exited
        N_simulated = current_sum
        removed_from_sim(:) = .false.   


        CALL CPU_TIME(t_temp)
        t_test = t_test + (t_temp-t0_test)



    ENDIF
END SUBROUTINE UPDATE_CELL_INDEX