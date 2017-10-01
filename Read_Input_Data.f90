SUBROUTINE INPUT_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE
    CHARACTER(80):: Header_name(18), line
!-----------------------------------------------------------------------
!*******************INITIALIZE CHARACTER*******************
!-----------------------------------------------------------------------

    ! this isn't strictly necessary, I'm listing the headers directly
    Header_name = (/    '*N_INITIAL_____', &
                        '*FN____________', &
                        '*N_CELLS_X_____', &
                        '*N_CELLLS_Y____', &
                        '*T_MAX_________', &
                        '*DT____________', &
                        '*DT_TO_SAVE____', &
                        '*INCLUDE_SOURCE', &
                        '*CLOSE_INLET___', &
                        '*CLOSE_OUTLET__', &
                        '*INCLUDE_GUN_BC', &
                        '*USE_HOMOG_GRID', &
                        '*DX_INIT_______', &
                        '*DX_FACTOR_____', &
                        '*DY_FACTOR_____', &
                        '*RESTART_SIM___', &
                        '*DIR_RESTART___', &
                        '*RESTART_INDEX_'/)         

!-----------------------------------------------------------------------
!*******************OPEN FILE******************************
!-----------------------------------------------------------------------
    OPEN(UNIT=100,FILE='Input/Input_Parameters')
!-----------------------------------------------------------------------
!*******************READ CONTENTS**************************
!-----------------------------------------------------------------------
    DO
        READ(100,*) line
        IF      (line == '*N_INITIAL_____') THEN
          ! initial density
          READ(100,*) n
        ELSE IF (line == '*NS_INITIAL____') THEN
          ! density of the source
          READ(100,*) ns
        ELSE IF (line == '*FN____________') THEN
          ! number of particles represented by each superparticle
          READ(100,*) Fn
        ELSE IF (line == '*N_CELLS_X_____') THEN
          ! number of cells in the x direction (if using homogeneous grid)
          READ(100,*) n_cells_x
          n_cells_vec(1) = n_cells_x
          nx = n_cells_x
        ELSE IF (line == '*N_CELLLS_Y____') THEN
          ! number of cells in the y direction (if using homogeneous grid)
          READ(100,*) n_cells_y
          n_cells_vec(2) = n_cells_y
          ny = n_cells_y
        ELSE IF (line == '*T_MAX_________') THEN
          ! end time of the simulation
          READ(100,*) tmax
        ELSE IF (line == '*DT____________') THEN
          ! timestep
          READ(100,*) dt
        ELSE IF (line == '*DT_TO_SAVE____') THEN
          ! how often to save simulation data (every 'dt-to-save'th timestep)
          READ(100,*) dt_to_save
        ELSE IF (line == '*INCLUDE_SOURCE') THEN
          ! whether to include source at inlet
          READ(100,*) include_source
        ELSE IF (line == '*USE_TWO_BEAMS_') THEN
          ! whether to include source at inlet
          READ(100,*) include_two_beams
        ELSE IF (line == '*BEAM_VELOCITY_') THEN
          ! whether to include source at inlet
          READ(100,*) v_beam
        ELSE IF (line == '*CLOSE_INLET___') THEN
          ! whether to put boundary at inlet (after any source is turned off)
          READ(100,*) close_inlet
        ELSE IF (line == '*CLOSE_OUTLET__') THEN
          ! whether to put boundary at the outlet
          READ(100,*) close_outlet
        ELSE IF (line == '*INCLUDE_GUN_BC') THEN
          ! include the boundaries of the gun geometry
          READ(100,*) include_gun_boundaries
        ELSE IF (line == '*ACCOMMODATION_') THEN
          ! include the boundaries of the gun geometry
          READ(100,*) accommodation
        ELSE IF (line == '*USE_HOMOG_GRID') THEN
          ! whether to use a homogeneous grid or the geometric one
          READ(100,*) use_homogenous_grid
        ELSE IF (line == '*DX_INIT_______') THEN
          ! initial cell-size (at center of inlet)
          READ(100,*) dx_0
        ELSE IF (line == '*DX_FACTOR_____') THEN
          ! relative cell-size at the outlet in the x-direction
          READ(100,*) dx_factor
        ELSE IF (line == '*DY_FACTOR_____') THEN
          ! relative cell-size at edges in the y-direction
          READ(100,*) dy_factor
        ELSE IF (line == '*RESTART_SIM___') THEN
          ! whether or not to restart a partially-completed simulation
          READ(100,*) restart_simulation
        ELSE IF (line == '*DIR_RESTART___') THEN
          ! directory of current data
          READ(100,*) dir_cur
        ELSE IF (line == '*RESTART_INDEX_') THEN
          ! index to restart from
          READ(100,*) it_restart
          EXIT
       ELSE
          ! Error Analyzing Input
          WRITE(*,*) 'Error Analyzing Input_Parameters'
          WRITE(*,*) line
          STOP
       END IF
    END DO
    CLOSE(100)
END SUBROUTINE INPUT_PARAMETERS_READIN



SUBROUTINE WALL_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE

    REAL(8):: temp

    ! top/bottom boundaries for the inlet
    OPEN(UNIT=1,FILE='Input/y_inlet.txt',STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) y_inlet
    CLOSE(1)

    IF (y_inlet(1) > y_inlet(2)) THEN
        temp = y_inlet(1)
        y_inlet(1) = y_inlet(2)
        y_inlet(2) = temp
    END IF


    ! xy coordinates of the walls (1-4 are the outer boundaries of the rectangle)
    OPEN(UNIT=1,FILE='Input/num_walls.txt',STATUS='old')
    READ(1,*) num_walls
    CLOSE(1)
    num_walls = num_walls + 4
    ! nw = num_walls
    ALLOCATE(x_walls(4,num_walls))

    OPEN(UNIT=1,FILE='Input/x_walls.txt',STATUS='old')
    READ(1,*) x_walls(:,5:num_walls)
    CLOSE(1)

    ! x_walls(:,6:num_walls) = 0

    ! WRITE(*,*) "y_inlet = ",y_inlet
    ! WRITE(*,*) "x_walls=",x_walls


END SUBROUTINE WALL_PARAMETERS_READIN






SUBROUTINE RESTART_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE
    CHARACTER(80):: line
    INTEGER i
    REAL::temp_real
!-----------------------------------------------------------------------
!*******************OPEN FILE******************************
!-----------------------------------------------------------------------
    WRITE(filename,"('Output/data/data.txt')")
    OPEN(UNIT=100,FILE=filename)

!-----------------------------------------------------------------------
!*******************READ CONTENTS**************************
!-----------------------------------------------------------------------
    DO
        READ(100,*) line
        IF      (line(1:8) == '*N_total') THEN
            i = 1
            DO
                READ(100,*) line
                IF (line(1:1) == "*") THEN
                    N_simulated = N_total(i-1)
                    EXIT
                ELSE
                    READ(line,*) N_total(i)
                    i = i+1
                ENDIF
            END DO

        ELSE IF   (line(1:24) == '*N_candidate_pairs_total') THEN
            i = 1
            DO
                READ(100,*) line
                IF (line(1:1) == "*") THEN
                    EXIT
                ELSE
                    READ(line,*) N_candidate_pairs_total(i)
                    i = i+1
                ENDIF
            END DO

        ELSE IF   (line(1:23) == '*N_accepted_pairs_total') THEN
            i = 1
            DO
                READ(100,*) line
                IF (line(1:1) == "*") THEN
                    EXIT
                ELSE
                    READ(line,*) N_accepted_pairs_total(i)
                    i = i+1
                ENDIF
            END DO

        ELSE IF   (line(1:19) == '*n_collisions_total') THEN
            i = 1
            DO
                READ(100,*) line
                IF (line(1:1) == "*") THEN
                    EXIT
                ELSE
                    READ(line,*) N_collisions_total(i)
                    i = i+1
                ENDIF
            END DO

        ELSE IF   (line(1:14) == '*N_added_total') THEN
            WRITE(*,*) "recognized N_added_total"
            i = 1
            DO
                READ(100,*) line
                IF (line(1:1) == "*") THEN
                    EXIT
                ELSE
                    READ(line,*) N_added_total(i)
                    i = i+1
                ENDIF
            END DO

        ELSE IF   (line(1:19) == '*ncp_remainder') THEN
            READ(100,"(E12.5)") ncp_remainder

        ELSE IF   (line(1:5) == "*t_BC") THEN
            WRITE(*,*) "recognized t_BC"
            READ(100,"(E12.5)") t_BC

        ELSE IF   (line(1:13) == "*t_collisions") THEN
            WRITE(*,*) "recognized t_collisions"
            READ(100,"(E12.5)") t_collisions

        ELSE IF   (line(1:7) == "*t_loop") THEN
            WRITE(*,*) "recognized t_loop"
            READ(100,"(E12.5)") t_loop

        ELSE IF   (line(1:8) == "*n_saved") THEN
            WRITE(*,*) "recognized n_saved"
            READ(100,*) t_loop

        ! add in timing variables and other things that need to be carried over


        ELSE IF (line(1:4) == '*END') THEN
            EXIT
        ELSE
            ! Error Analyzing Input
            ! WRITE(*,*) 'Error Analyzing Input_Parameters'
            ! WRITE(*,*) line
            ! STOP
       END IF
    END DO
    CLOSE(100)
END SUBROUTINE RESTART_PARAMETERS_READIN







