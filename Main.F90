MODULE CONTAIN
    IMPLICIT NONE
    SAVE

!-----------------------------------------------------------------------
!*******************INPUT PARAMETERS*************************
!-----------------------------------------------------------------------
    REAL(8)::n,ns,Fn
    INTEGER::n_cells_x,n_cells_y
    INTEGER::dt_to_save
    REAL(8)::tmax,dt
    REAL(8)::dx_0,dx_factor,dy_factor
    LOGICAL::include_source,close_inlet,close_outlet,include_gun_boundaries,use_homogenous_grid,restart_simulation
    CHARACTER(80)::dir_cur
    INTEGER::it_restart

!-----------------------------------------------------------------------
!*******************DYNAMIC ARRAYS*************************
!-----------------------------------------------------------------------
    REAL(8), ALLOCATABLE, DIMENSION(:):: x_cells_vec, y_cells_vec, dx_cells_vec, dy_cells_vec, i_range, t_vec
    INTEGER, ALLOCATABLE, DIMENSION(:):: N_total,N_candidate_pairs_total,N_accepted_pairs_total,N_collisions_total,N_added_total
    REAL(8), ALLOCATABLE, DIMENSION(:):: counter_vec
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: Vc, x_walls, vr_max, x_vec,v_vec, ncp_remainder
    INTEGER, ALLOCATABLE, DIMENSION(:,:):: i_cell_vec, i_cell_vec_prev, Npc_slice, starting_index, Npc_added
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: x_vec_prev,v_vec_prev, xs_vec,vs_vec,xs_vec_prev,vs_vec_prev
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: x_vec_unsorted,v_vec_unsorted, i_cell_vec_unsorted
    LOGICAL,ALLOCATABLE, DIMENSION(:):: reflected_in, reflected_out, in_column, in_cell, removed_from_sim, entered_sim
    LOGICAL,ALLOCATABLE, DIMENSION(:):: removed_from_sim_unsorted
    INTEGER,ALLOCATABLE, DIMENSION(:):: i_counting, i_column, i_cur

!-----------------------------------------------------------------------
!*******************MISC. VARIBLES*************************
!-----------------------------------------------------------------------
    REAL(8):: m_g,d_g,vth,c_s,vr_max_0,v_avg, xmin,xmax,ymin,ymax,ymid, n_inf, V_total
    REAL(8):: ws,ts,hs,Vs,xs_min,xs_max,ys_min,ys_max,t, N_candidate_pairs_real
    REAL(8):: Nc0,Nc_sim,m_r,collision_ratio
    REAL(8):: alpha_x,alpha_y,neg_offset,pos_offset
    REAL (8):: t0,t0_BC,t0_collisions,t0_loop,t_temp,t_final,t_BC,t_collisions,t_loop, t0_test,t_test
    INTEGER:: nmax, nx, ny, n_cells, N_all,N_simulated, nt, n_saved, nw, N_expected, N_array, Num_s, N_entered, cx,cy,Npc_max, ii
    INTEGER:: N_candidate_pairs,N_accepted_pairs,Npc_cur, num_walls, N_collisions, N_added, N_removed
    REAL(8), DIMENSION(2,2):: x_lim
    INTEGER, DIMENSION(2):: n_cells_vec
    LOGICAL:: file_exists
    CHARACTER(80)::filename


!-----------------------------------------------------------------------
!*******************TEMPORARY/SCRATCH VARIBLES*************************
!-----------------------------------------------------------------------
    REAL(8):: a1,b1,c1
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: a_mat,b_mat,c_mat
    REAL(8), ALLOCATABLE, DIMENSION(:):: a_vec,b_vec,c_vec
    INTEGER, ALLOCATABLE, DIMENSION(:):: i_temp_vec
    
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

! gun geometry dimensions
    REAL(8), PARAMETER:: inlet_height = 1.d0
    ! REAL(8), PARAMETER:: inlet_height = 0.5
    REAL(8), PARAMETER:: inlet_length = .1d0
    REAL(8), PARAMETER:: gun_length = .26d0
    REAL(8), PARAMETER:: gun_height = .05d0
    REAL(8), PARAMETER:: outlet_height = 1.d0
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
    CHARACTER(10)::str_file_num, line

    ! ALLOCATE(a_mat(6,2))
    ! ALLOCATE(b_mat(6,2))
    ! a_mat(:,1) = 1
    ! a_mat(:,2) = 2
    ! WRITE(*,*) "a_mat = ",a_mat

    ! ! save position data
    ! WRITE(filename,"('Output/data/a.txt')")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! WRITE(1) a_mat
    ! ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! ! WRITE(1,"(E12.5)") x_vec
    ! CLOSE(1)


    ! OPEN(UNIT=1,FILE=filename,STATUS='old', FORM='unformatted')! ,access='direct',recl=4,iostat=ok)
    ! READ(1) b_mat
    ! CLOSE(1)
    ! WRITE(*,*) "b_mat = ",b_mat
    ! WRITE(*,*) "shape(b_mat) = ",SHAPE(b_mat)

    ! WRITE(filename,"('Output/data/data.txt')")
    ! OPEN(UNIT=100,FILE=filename)
    ! READ(100,*) line
    ! CLOSE(100)

    ! WRITE(*,*)"line=",line
    ! WRITE(*,*)"line(3:4)=",line(3:4)




!-----------------------------------------------------------------------
!*******************READ INPUT FILE************************
!-----------------------------------------------------------------------
    CALL INPUT_PARAMETERS_READIN
!-----------------------------------------------------------------------
!*******************SET INITIAL CONDITIONS*****************
!-----------------------------------------------------------------------
    CALL INITIALIZE
! !-----------------------------------------------------------------------
! !*******************MAIN LOOP******************************
! !-----------------------------------------------------------------------
         
    DO ii = it_restart+2,nt
        t = t_vec(ii)
        IF (MOD(ii-1,dt_to_save) == 0) THEN
            WRITE(*,*) ii-1,", t=",t,", N=",N_simulated,' ------------------------'
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
        ! IF (MOD(ii-1,dt_to_save) == 0) THEN
        !     ! WRITE(*,*) "x=",x_vec(1:10,1)
        !     ! WRITE(*,*) "xp=",x_vec_prev(1:10,1)
        !     ! WRITE(*,*) "v=",v_vec(1:10,1)
        !     a1 = COUNT(v_vec(:,1) > 1.d-5)
        !     b1 = COUNT(v_vec(:,1) < -1.d-5)
        !     c1 = COUNT(ABS(v_vec(:,1)) <= 1.d-5)
        !     WRITE(*,*) "a1=",a1
        !     WRITE(*,*) "b1=",b1
        !     WRITE(*,*) "c1=",c1
        ! END IF



        ! Boundary Condition implementation -------------------------------------------
        ! (removing exiting particles should be absorbed into BC implementation)
        IF (N_simulated > 0) THEN
            CALL SPECULAR_REFLECTION('SIM_PARTICLES___',N_array,0)
        END IF


        ! Input from Source/Reservoir -------------------------------------------------
        IF (include_source .EQV. .true.) THEN
            CALL INITIALIZE_SOURCE
            N_added_total(ii) = N_entered            
        END IF


        IF (N_simulated > 0) THEN
            ! Update Cell Index -----------------------------------------------------------
            CALL UPDATE_CELL_INDEX


            ! Collisions ------------------------------------------------------------------
            CALL RUN_COLLISIONS
        END IF
        

        ! update and save current data ---------------------------------------------
        N_total(ii) = N_simulated


        ! NOTE: saving with the label (ii-1)
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
    t_final = t_temp-t0

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
    WRITE(*,*) "N_expected =",N_expected
    WRITE(*,*) "N_added = ",N_added
    WRITE(*,*) "peak N = ",MAXVAL(N_total)
    WRITE(*,*) "(nx,ny) = (",((/ nx, ny /)),")"
    WRITE(*,*) "nt = ",nt
    WRITE(*,*) "n_saved = ",n_saved
    WRITE(*,*) "computation time (total) = ",t_final
    WRITE(*,*) "computation time (BC's) = ",t_BC
    WRITE(*,*) "computation time (collisions) = ",t_collisions
    WRITE(*,*) "computation time (looping in collisions)=",t_loop
    WRITE(*,*) "computation time (test)=",t_test

    WRITE(*,*) '------ done ------'

    STOP
END PROGRAM MAIN





