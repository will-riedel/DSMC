SUBROUTINE INPUT_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE
    CHARACTER(80):: Header_name(19), line
!-----------------------------------------------------------------------
!*******************INITIALIZE CHARACTER*******************
!-----------------------------------------------------------------------

    ! this isn't strictly necessary, I'm listing the headers directly
    Header_name = (/    '*N_INITIAL_____', &
                        '*FN____________', &
                        '*N_CELLS_X_____', &
                        '*N_CELLS_Y_____', &
                        '*T_MAX_________', &
                        '*DT____________', &
                        '*DT_TO_SAVE____', &
                        '*INCLUDE_SOURCE', &
                        '*CLOSE_INLET___', &
                        '*CLOSE_OUTLET__', &
                        '*INCLUDE_GUN_BC', &
                        '*USE_HOMOG_GRID', &
                        '*DX_INIT_______', &
                        '*DY_INIT_______', &
                        '*DX_FACTOR_____', &
                        '*DY_FACTOR_____', &
                        '*RESTART_SIM___', &
                        '*DIRECTORY_CUR_', &
                        '*RESTART_INDEX_'/)         

!-----------------------------------------------------------------------
!*******************OPEN FILE******************************
!-----------------------------------------------------------------------
    OPEN(UNIT=100,FILE='Input/Input_Parameters')
!-----------------------------------------------------------------------
!*******************READ CONTENTS**************************
!-----------------------------------------------------------------------
    DO
        ! WRITE(*,*) "reading..."
        READ(100,*) line
        IF      (line == '*N_INITIAL_____') THEN
          ! initial density
          READ(100,*) n
        ELSE IF (line == '*NS_INITIAL____') THEN
          ! density of the source
          READ(100,*) ns
        ELSE IF (line == '*N_INITIAL_B___') THEN
          ! initial density
          READ(100,*) n_b
        ELSE IF (line == '*NS_INITIAL_B__') THEN
          ! density of the source
          READ(100,*) ns_b
        ELSE IF (line == '*FN____________') THEN
          ! number of particles represented by each superparticle
          READ(100,*) Fn
        ELSE IF (line == '*N_CELLS_X_____') THEN
          ! number of cells in the x direction (if using homogeneous grid)
          READ(100,*) n_cells_x
          n_cells_vec(1) = n_cells_x
          nx = n_cells_x
        ELSE IF (line == '*N_CELLS_Y_____') THEN
          ! number of cells in the y direction (if using homogeneous grid)
          READ(100,*) n_cells_y
          n_cells_vec(2) = n_cells_y
          ny = n_cells_y
        ELSE IF (line == '*INITIAL_DISTR_') THEN
          ! number of cells in the y direction (if using homogeneous grid)
          READ(100,*) initial_distribution
        ELSE IF (line == '*GEOMETRY_TYPE_') THEN
          ! cartesian geometry or cylindrical/axial symmetry
          ! either CYLINDRICAL or CARTESIAN right now, although really anything but cylindrical is the same
          READ(100,*) geometry_type
        ELSE IF (line == '*RADIAL_WEIGHTF') THEN
          ! cartesian geometry or cylindrical/axial symmetry
          READ(100,*) RWF_input
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
        ! ELSE IF (line == '*USE_HOMOG_GRID') THEN
        !   ! whether to use a homogeneous grid or the geometric one
        !   READ(100,*) use_homogenous_grid
        ! use_homogenous_grid = .false.
        ELSE IF (line == '*INCLUDE_COLL__') THEN
          ! whether to include source at inlet
          READ(100,*) include_collisions
        ELSE IF (line == '*DX_INIT_______') THEN
          ! initial cell-size (at center of inlet)
          READ(100,*) dx_0
        ELSE IF (line == '*DY_INIT_______') THEN
          ! initial cell-size (at center of inlet)
          READ(100,*) dy_0
        ELSE IF (line == '*DX_FACTOR_____') THEN
          ! relative cell-size at the outlet in the x-direction
          READ(100,*) dx_factor
        ELSE IF (line == '*DY_FACTOR_____') THEN
          ! relative cell-size at edges in the y-direction
          READ(100,*) dy_factor
        ELSE IF (line == '*DX_INLET______') THEN
          ! dx in high-density inlet region
          READ(100,*) dx_inlet
        ELSE IF (line == '*RESTART_SIM___') THEN
          ! whether or not to restart a partially-completed simulation
          READ(100,*) restart_simulation
        ELSE IF (line == '*X_GRID_TYPE___') THEN
          ! type of grid to use for columns
          READ(100,*) x_grid_type
        ELSE IF (line == '*Y_GRID_TYPE___') THEN
          ! type of grid to use for columns
          READ(100,*) y_grid_type
        ELSE IF (line == '*SOURCE_TYPE___') THEN
          ! type of source input
          READ(100,*) source_type
        ELSE IF (line == '*RESTART_INDEX_') THEN
          ! index to restart from
          READ(100,*) it_restart
        ELSE IF (line == '*DIRECTORY_CUR_') THEN
          ! directory of current data
          READ(100,*) dir_cur
          ! WRITE(*,*) "found directory"
          dir_cur = "Output/" // dir_cur
          dir_cur_length = LEN_TRIM(dir_cur)
          EXIT
        ELSE
          ! Error Analyzing Input
          WRITE(*,*) 'Error Analyzing Input_Parameters'
          WRITE(*,*) line
          STOP
       END IF
    END DO
    CLOSE(100)

    CALL WALL_PARAMETERS_READIN

END SUBROUTINE INPUT_PARAMETERS_READIN



SUBROUTINE WALL_PARAMETERS_READIN
    USE CONTAIN
    IMPLICIT NONE

    REAL(8):: temp

    ! top/bottom boundaries for the inlet
    OPEN(UNIT=1,FILE='Input/y_inlet.txt',STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) y_inlet
    CLOSE(1)

    OPEN(UNIT=1,FILE='Input/x_inlet.txt',STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) x_inlet
    CLOSE(1)

    IF (y_inlet(1) > y_inlet(2)) THEN
        temp = y_inlet(1)
        y_inlet(1) = y_inlet(2)
        y_inlet(2) = temp
    END IF


    ! xy coordinates of the walls (1-4 are the outer boundaries of the rectangle)
    IF (include_gun_boundaries .EQV. .true.) THEN
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

    ELSE
        num_walls = 4
        ALLOCATE(x_walls(4,4))

      END IF


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
    ! WRITE(filename,"('Output/data/data.txt')")
    ! WRITE(filename,"('/data.txt')")
    WRITE(filename,"('/data_',I7.7,'.txt')") (it_restart)

    filename = dir_cur(1:dir_cur_length) // filename
    WRITE(*,*) "got here"
    WRITE(*,*) filename
    OPEN(UNIT=100,FILE=filename)

!-----------------------------------------------------------------------
!*******************READ CONTENTS**************************
!-----------------------------------------------------------------------
    DO
        READ(100,*) line
        ! IF      (line(1:8) == '*N_total') THEN
        !     i = 1
        !     DO
        !         READ(100,*) line
        !         IF (line(1:1) == "*") THEN
        !             N_simulated = N_total(i-1)
        !             EXIT
        !         ELSE
        !             READ(line,*) N_total(i)
        !             i = i+1
        !         ENDIF
        !     END DO

        ! ELSE IF   (line(1:24) == '*N_candidate_pairs_total') THEN
        !     i = 1
        !     DO
        !         READ(100,*) line
        !         IF (line(1:1) == "*") THEN
        !             EXIT
        !         ELSE
        !             READ(line,*) N_candidate_pairs_total(i)
        !             i = i+1
        !         ENDIF
        !     END DO

        ! ELSE IF   (line(1:23) == '*N_accepted_pairs_total') THEN
        !     i = 1
        !     DO
        !         READ(100,*) line
        !         IF (line(1:1) == "*") THEN
        !             EXIT
        !         ELSE
        !             READ(line,*) N_accepted_pairs_total(i)
        !             i = i+1
        !         ENDIF
        !     END DO

        ! ELSE IF   (line(1:19) == '*n_collisions_total') THEN
        !     i = 1
        !     DO
        !         READ(100,*) line
        !         IF (line(1:1) == "*") THEN
        !             EXIT
        !         ELSE
        !             READ(line,*) N_collisions_total(i)
        !             i = i+1
        !         ENDIF
        !     END DO

        ! ELSE IF   (line(1:14) == '*N_added_total') THEN
        !     WRITE(*,*) "recognized N_added_total"
        !     i = 1
        !     DO
        !         READ(100,*) line
        !         IF (line(1:1) == "*") THEN
        !             EXIT
        !         ELSE
        !             READ(line,*) N_added_total(i)
        !             i = i+1
        !         ENDIF
        !     END DO

        ! IF   (line(1:19) == '*ncp_remainder') THEN
        !     READ(100,"(E12.5)") ncp_remainder
        IF   (line(1:8) == "*t_total") THEN
            WRITE(*,*) "recognized t_total"
            READ(100,"(E12.5)") t_total

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
            READ(100,*) n_saved

        ELSE IF   (line(1:12) == "*N_simulated") THEN
            WRITE(*,*) "recognized N_simulated"
            READ(100,*) N_simulated

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


    WRITE(filename,"('/N_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) N_total(1:(it_restart+1))
    CLOSE(1)


    ! WRITE(*,*) N_total(1:(n_saved*dt_to_save+5))
    ! WRITE(*,*) "size=",SHAPE(N_total)
    ! WRITE(*,*) "it_restart+1=",it_restart+1
    ! WRITE(*,*) "n_saved*dt_to_save+1=",n_saved*dt_to_save+1
    ! WRITE(*,*) "N_simulated=",N_simulated


    WRITE(filename,"('/n_collisions_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) n_collisions_total(1:(it_restart+1))
    CLOSE(1)

    WRITE(filename,"('/N_candidate_pairs_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) N_candidate_pairs_total(1:(it_restart+1))
    CLOSE(1)

    WRITE(filename,"('/N_accepted_pairs_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) N_accepted_pairs_total(1:(it_restart+1))
    CLOSE(1)

    WRITE(filename,"('/N_added_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) N_added_total(1:(it_restart+1))
    CLOSE(1)

    ! WRITE(filename,"('/ncp_remainder.txt')")
    ! filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    ! READ(1,*) ncp_remainder
    ! CLOSE(1)




    ! for now, dir_cur = "Output/data"
    ! WRITE(filename,"('Output/data/x_',I7.7,'.txt')") (it_restart)
    WRITE(filename,"('/x_',I7.7,'.txt')") (it_restart)
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) x_vec(1:N_simulated,:)
    CLOSE(1)

    ! WRITE(filename,"('Output/data/v_',I7.7,'.txt')") (it_restart)
    WRITE(filename,"('/v_',I7.7,'.txt')") (it_restart)
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) v_vec(1:N_simulated,:)
    CLOSE(1)

    ! WRITE(filename,"('Output/data/i_',I7.7,'.txt')") (it_restart)
    WRITE(filename,"('/i_',I7.7,'.txt')") (it_restart)
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) i_cell_vec(1:N_simulated,:)
    CLOSE(1)

    ! WRITE(filename,"('Output/data/Npc_',I7.7,'.txt')") (it_restart)
    WRITE(filename,"('/Npc_',I7.7,'.txt')") (it_restart)
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,STATUS='old')! ,access='direct',recl=4,iostat=ok)
    READ(1,*) Npc_slice
    CLOSE(1)




END SUBROUTINE RESTART_PARAMETERS_READIN







