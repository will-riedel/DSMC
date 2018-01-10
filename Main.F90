MODULE CONTAIN
    IMPLICIT NONE
    SAVE

!-----------------------------------------------------------------------
!*******************INPUT PARAMETERS*************************
!-----------------------------------------------------------------------
    REAL(8)::n,ns,Fn, n_b,ns_b, x_split
    INTEGER::n_cells_x,n_cells_y
    INTEGER::dt_to_save
    REAL(8)::tmax,dt
    REAL(8)::dx_0,dy_0,dx_factor,dy_factor,dx_inlet, x_inlet
    LOGICAL::include_source,close_inlet,close_outlet,include_gun_boundaries,use_homogenous_grid,restart_simulation
    LOGICAL::include_two_beams
    CHARACTER(80)::dir_cur,dir_temp,x_grid_type,y_grid_type,source_type,initial_distribution
    INTEGER::it_restart

!-----------------------------------------------------------------------
!*******************DYNAMIC ARRAYS*************************
!-----------------------------------------------------------------------
    REAL(8), ALLOCATABLE, DIMENSION(:):: x_cells_vec, y_cells_vec, y_cells_half, dx_cells_vec, dy_cells_vec
    REAL(8), ALLOCATABLE, DIMENSION(:):: y_cells_vec_b, dy_cells_vec_b
    REAL(8), ALLOCATABLE, DIMENSION(:):: counter_vec, i_range_x,i_range_y, t_vec
    INTEGER, ALLOCATABLE, DIMENSION(:):: N_total,N_candidate_pairs_total,N_accepted_pairs_total,N_collisions_total,N_added_total
    INTEGER, ALLOCATABLE, DIMENSION(:):: flux_upstream_total,flux_downstream_total
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: Vc, x_walls, vr_max, x_vec,v_vec!,ncp_remainder
    INTEGER, ALLOCATABLE, DIMENSION(:,:):: i_cell_vec, i_cell_vec_prev, Npc_slice
    INTEGER, ALLOCATABLE, DIMENSION(:,:):: starting_index, final_index, Npc_added, i_cell_lim_x, i_cell_lim_y
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: x_vec_prev,v_vec_prev, xs_vec,vs_vec,xs_vec_prev,vs_vec_prev
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: x_vec_unsorted,v_vec_unsorted, i_cell_vec_unsorted
    LOGICAL,ALLOCATABLE, DIMENSION(:):: reflected_in, reflected_out, in_column, in_cell, removed_from_sim, entered_sim
    LOGICAL,ALLOCATABLE, DIMENSION(:):: removed_from_sim_unsorted
    INTEGER,ALLOCATABLE, DIMENSION(:):: i_counting, i_column, i_cur

    REAL(8),ALLOCATABLE,DIMENSION(:,:):: xr_vec,xr_vec_prev, xr_vec_new
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: vr_vec,vr_vec_prev, vr_vec_new
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: xr_walls,xy0,xyt
    LOGICAL,ALLOCATABLE,DIMENSION(:,:):: collision_occured
    REAL(8),ALLOCATABLE,DIMENSION(:,:):: collision_dt
    REAL(8),ALLOCATABLE,DIMENSION(:):: x0,y0,xt,yt,m,b,xc,yc, dt_cross,min_collision_dt, wall_angle_vec, rn_vec
    LOGICAL,ALLOCATABLE,DIMENSION(:):: first_collision,crossed,collision_occured_any
    INTEGER,ALLOCATABLE,DIMENSION(:):: i_cross, i_first


!-----------------------------------------------------------------------
!*******************MISC. VARIBLES*************************
!-----------------------------------------------------------------------
    REAL(8):: m_g,d_g,vth,c_s,vr_max_0,v_avg, xmin,xmax,ymin,ymax,ymid, n_inf, V_total, v_beam
    REAL(8):: ws,ts,hs,Vs,xs_min,xs_max,xs2_min,xs2_max,ys_min,ys_max,t, N_candidate_pairs_real
    REAL(8):: Nc0,Nc_sim,m_r,collision_ratio, Num_s_exact, Num_s_frac, b_source_A, b_source_B, b_source_barrier,Theta_source
    REAL(8):: Num_s_exact_b, Num_s_frac_b
    REAL(8):: alpha_x,alpha_y,neg_offset,pos_offset, accommodation 
    REAL (8):: t0_total,t0_BC,t0_collisions,t0_loop,t_temp,t_total,t_BC,t_collisions,t_loop, t0_index,t_index
    ! REAL (8):: t0_BC1,t0_BC2,t0_BC3,t0_BC4,t0_BC5,t_BC1,t_BC2,t_BC3,t_BC4,t_BC5
    REAL(8)::t0_source,t_source,t0_init,t_init,t0_save,t_save
    INTEGER:: N_all,N_simulated, nt, n_saved, nw, N_expected, N_array, Num_s, Num_s_b, N_entered, cx,cy,Npc_max, ii
    INTEGER:: N_init_a,N_init_b
    INTEGER:: nmax, nmax_left,nmax_right, nx, ny, ny_b, n_cells, cx_min_collisions, cy_min_collisions
    INTEGER:: N_candidate_pairs,N_accepted_pairs,Npc_cur, num_walls, N_collisions, N_added, N_removed, N_specular, N_diffuse
    REAL(8), DIMENSION(2,2):: x_lim, Rotation_mat_neg, Rotation_mat_pos, x_source_corners
    REAL(8), DIMENSION(2):: y_inlet,xy_w
    REAL(8), DIMENSION(3):: alpha_vec
    INTEGER, DIMENSION(2):: n_cells_vec
    LOGICAL:: file_exists, finding_wall_cells
    CHARACTER(80)::filename

    CHARACTER(16):: string_in
    INTEGER:: Num_r,counter, dir_cur_length,dir_temp_length, cell_lim_buffer
    REAL(8):: xw1,xw2,yw1,yw2,xw1_0,yw1_0,xw2_0,yw2_0,m_w,b_w,Theta



!-----------------------------------------------------------------------
!*******************TEMPORARY/SCRATCH VARIBLES*************************
!-----------------------------------------------------------------------
    REAL(8):: a1,b1,c1
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: a_mat,b_mat,c_mat
    REAL(8), ALLOCATABLE, DIMENSION(:):: a_vec,b_vec,c_vec
    INTEGER, ALLOCATABLE, DIMENSION(:):: i_temp_vec
    LOGICAL, ALLOCATABLE, DIMENSION(:):: a_log_vec
    
END MODULE CONTAIN




MODULE PROPERTIES
    IMPLICIT NONE
    SAVE

! fundamental constants
    REAL(8), PARAMETER:: Pi = 4.d0*atan(1.d0)
    REAL(8), PARAMETER:: k_b = 1.3806503d-23
    REAL(8), PARAMETER:: m_p= 1.67262178d-27

! gas properties
    REAL(8), PARAMETER:: T_g = 293.d0
    REAL(8), PARAMETER:: dH2 = 2.89d-10
    REAL(8), PARAMETER:: dN2 = 3.64d-10

! gun geometry dimensions
    REAL(8), PARAMETER:: inlet_height = 1.d0
    ! REAL(8), PARAMETER:: inlet_height = 0.01d0
    ! REAL(8), PARAMETER:: inlet_height = 0.5
    ! REAL(8), PARAMETER:: inlet_height = 0.025
    REAL(8), PARAMETER:: inlet_length = .1d0
    REAL(8), PARAMETER:: gun_length = .26d0
    REAL(8), PARAMETER:: gun_height = .05d0
    REAL(8), PARAMETER:: outlet_height = 1.d0
    ! REAL(8), PARAMETER:: outlet_height = 0.05
    REAL(8), PARAMETER:: outlet_length = .25d0

! number of dimensions (keep at 2 right now, even for 1-D simulation)
    INTEGER, PARAMETER:: ndim = 2

END MODULE PROPERTIES





! ##############################################################################################
! ##############################################################################################
! ##############################################################################################
! ##############################################################################################
! ##############################################################################################

PROGRAM MAIN
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER:: i
    ! CHARACTER(80)::filename

    WRITE(*,*) "Starting computation..."

    ! N_array = 1.d8
    ! ALLOCATE(x_vec(N_array,2))
    ! ALLOCATE(first_collision(N_array))
    ! ! 1 --------------------------------
    ! CALL CPU_TIME(t0_test)
    ! x_vec(1:N_array,1) = 0.
    ! x_vec(1:N_array,2) = 0.
    ! first_collision(1:N_array) = .true.
    ! ! WHERE(x_vec(1:N_array,1) == 0.)
    ! WHERE(first_collision(1:N_array) .eqv. .true.)
    !     x_vec(1:N_array,2) = 1.
    ! ELSEWHERE
    !     x_vec(1:N_array,2) = 2.
    ! ENDWHERE

    ! CALL CPU_TIME(t_temp)
    ! WRITE(*,*) "t1 = ", (t_temp-t0_test)

    ! ! 2 --------------------------------
    ! CALL CPU_TIME(t0_test)
    ! x_vec(1:N_array,1) = 0.
    ! x_vec(1:N_array,2) = 0.
    ! first_collision(1:N_array) = .true.
    ! DO i = 1,N_array
    !     ! IF (x_vec(i,1) == 0.) THEN
    !     IF (first_collision(i) .eqv. .true.) THEN
    !         x_vec(i,2) = 1.
    !     ELSE
    !         x_vec(i,2) = 2.
    !     END IF
    ! END DO

    ! CALL CPU_TIME(t_temp)
    ! WRITE(*,*) "t2 = ", (t_temp-t0_test)



!-----------------------------------------------------------------------
!*******************READ INPUT FILE************************
!-----------------------------------------------------------------------
! WRITE(*,*) "GH 1"
    CALL INPUT_PARAMETERS_READIN
!-----------------------------------------------------------------------
!*******************SET INITIAL CONDITIONS*****************
!-----------------------------------------------------------------------
! WRITE(*,*) "GH 2"
    CALL INITIALIZE
! !-----------------------------------------------------------------------
! !*******************MAIN LOOP******************************
! !-----------------------------------------------------------------------

    ! WRITE(*,*) "ns,n,N_expected,N_array=",ns,n,N_expected
    ! WRITE(*,*) "nx,ny=",nx,ny
    ! WRITE(*,*) "GH 3"
    DO ii = it_restart+2,nt
        t = t_vec(ii)
        IF (MOD(ii-1,dt_to_save) == 0) THEN
            WRITE(*,*) ii-1,", t=",t,", N=",N_simulated,' ------------------------'
            ! WRITE(*,*) "t_collisions,t_BC,t_loop,t_test=",t_collisions,t_BC,t_loop,t_test
        END IF

        IF ( (t>ts) .and. (include_source .EQV. .true.) ) THEN
            close_inlet = .true.
            include_source = .false.
        END IF

        x_vec_prev = x_vec
        v_vec_prev = v_vec

        ! collisionless motion --------------------------------------------------------
        ! x_vec = x_vec + dt*v_vec(:,1:2)
        IF (N_simulated > 0) THEN
            x_vec(1:N_simulated,:) = x_vec(1:N_simulated,:) + dt*v_vec(1:N_simulated,1:2)
        END IF

        ! WRITE(*,*) "GH 4"
        ! Boundary Condition implementation -------------------------------------------
        ! (removing exiting particles should be absorbed into BC implementation)
        IF (N_simulated > 0) THEN
            Num_r = N_simulated
            counter = 0
            CALL COMPUTE_REFLECTION
            ! CALL SPECULAR_REFLECTION_1D
        END IF

        ! WRITE(*,*) "GH 5"
        ! Input from Source/Reservoir -------------------------------------------------
        CALL INITIALIZE_SOURCE

        ! WRITE(*,*) "GH 6"
        IF (N_simulated > 0) THEN

            ! Update Cell Index -----------------------------------------------------------
            CALL UPDATE_CELL_INDEX
            ! CALL UPDATE_CELL_INDEX_TEMP

            ! WRITE(*,*) "GH 7"
            ! Collisions ------------------------------------------------------------------
            CALL RUN_COLLISIONS

        END IF
        
        ! update and save current data ---------------------------------------------
        N_total(ii) = N_simulated

        ! NOTE: saving with the label (ii-1)
        ! WRITE(*,*) "GH 8"
        IF ( ( MOD(ii-1,dt_to_save) == 0 ) .or. ( ii == nt ) ) THEN
            CALL SAVE_DATA
        END IF

    END DO



    ! Final processing/printing

    N_collisions = SUM(N_collisions_total)
    N_candidate_pairs = SUM(N_candidate_pairs_total)
    N_accepted_pairs = SUM(N_accepted_pairs_total)
    N_added = SUM(N_added_total)
    call CPU_TIME(t_temp)
    ! t_total = t_temp-t0_total
    t_total = t_total + (t_temp-t0_total)

    ! calculated collision rate compared to analytical solution (at equilbrium)
    m_r = m_g/2
    Nc0 = n**2/2*d_g**2*(8*Pi*k_b*T_g/m_r)**.5*(V_total)
    Nc_sim = N_collisions*Fn/(nt*dt)
    collision_ratio = Nc_sim/Nc0


    WRITE(*,*) "N_collisions = ",N_collisions
    WRITE(*,*) "N_candidate_pairs =",N_candidate_pairs
    WRITE(*,*) "N_accepted_pairs =",N_accepted_pairs
    WRITE(*,*) "collision acceptance rate = ",N_collisions*1.0/N_candidate_pairs
    WRITE(*,*) "Computed Collision Rate / Analytic Collision rate = ", collision_ratio
    WRITE(*,*) "N_diffuse, N_specular (not updated), Effective Accomm. Coeff. = ", & 
    N_diffuse, N_specular, (N_diffuse+0.)/(N_diffuse+N_specular)
    WRITE(*,*) "N_expected =",N_expected
    WRITE(*,*) "N_added = ",N_added
    WRITE(*,*) "peak N = ",MAXVAL(N_total)
    WRITE(*,*) "(nx,ny) = (",((/ nx, ny /)),")"
    WRITE(*,*) "nt = ",nt
    WRITE(*,*) "n_saved = ",n_saved
    WRITE(*,*) "computation time (total) = ",t_total
    WRITE(*,*) "computation time (BC's) = ",t_BC
    WRITE(*,*) "computation time (collisions) = ",t_collisions
    WRITE(*,*) "computation time (looping in collisions)=",t_loop
    WRITE(*,*) "computation time (indexing)=",t_index
    WRITE(*,*) "computation time (source)=",t_source
    WRITE(*,*) "computation time (initialization)=",t_init
    WRITE(*,*) "computation time (saving)=",t_save

    WRITE(*,*) "flux down, up = ",SUM(flux_downstream_total),SUM(flux_upstream_total)

    ! WRITE(*,*) "t0_BC1 = ",t_BC1
    ! WRITE(*,*) "t0_BC2 = ",t_BC2
    ! WRITE(*,*) "t0_BC3 = ",t_BC3
    ! WRITE(*,*) "t0_BC4 = ",t_BC4
    ! WRITE(*,*) "t0_BC5 = ",t_BC5



    WRITE(*,*) '------ done ------'

    STOP
END PROGRAM MAIN





