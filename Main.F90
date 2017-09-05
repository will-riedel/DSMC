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
    CHARACTER(80)::filename
    CHARACTER(10)::str_file_num

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
    ! DO i = it_restart+2,50
        ! print current step

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

    ! filename = "Output/data/mat.txt"
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! WRITE(1) v_vec(:,1)
    ! CLOSE(1)


    ! WRITE(filename,"('mat_',I7.7,'.txt')") 10
    ! ! WRITE(str_file_num,'(I0)') 10
    ! ! filename = "Output/data/mat_" // str_file_num //".txt"
    ! WRITE(*,*) filename

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
    WRITE(*,*) "computation time (total) = ",t_final
    WRITE(*,*) "computation time (BC's) = ",t_BC
    WRITE(*,*) "computation time (collisions) = ",t_collisions
    WRITE(*,*) "computation time (looping in collisions)=",t_loop
    WRITE(*,*) "computation time (test)=",t_test

    WRITE(*,*) '------ done ------'

    STOP
END PROGRAM MAIN



SUBROUTINE INITIALIZE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    REAL(8)::dx,dy
    INTEGER::i,j,current_sum
    ! This suboutine initializes various matrices/parameters in the simulation

     ! WRITE(*,*) "ii=",ii    
     ! WRITE(*,*) "it_restart=",it_restart    
     ! WRITE(*,*) "nt=",nt

    IF (restart_simulation .EQV. .false.) THEN
        it_restart = -1
    END IF


    ! simulation geometry and gas properties -----------------------------------------    
    m_g = 2*m_p
    d_g = dH2
    vth = SQRT(k_b*T_g/m_g)
    c_s = Pi*d_g**2
    vr_max_0 = 5*vth
    ! x_lim(1,:) = (/ 0.d0,inlet_length+gun_length+outlet_length /)
    x_lim = transpose(reshape( (/ 0.d0,inlet_length+gun_length+outlet_length,   0.d0,outlet_height /) , (/2,2/) ))
    xmin = x_lim(1,1)
    xmax = x_lim(1,2)
    ymin = x_lim(2,1)
    ymax = x_lim(2,2)    

    ! nx = n_cells_vec(1)
    ! ny = n_cells_vec(2)
    ny = 1
    cy  = 1


    IF (use_homogenous_grid .EQV. .true.) THEN
        ! set up equally spaced grid points
        ALLOCATE( x_cells_vec(nx) )
        ALLOCATE( y_cells_vec(ny) )
        ! dx = ( x_lim(1,2) - x_lim(1,1) ) / (nx)
        dx = ( xmax - xmin ) / (nx)
        DO i = 1 , nx
            ! x_cells_vec(i) = x_lim(1,1) + dx*(i-1)
            x_cells_vec(i) = xmin + dx*(i-1)
        END DO
        ! dy = ( x_lim(2,2) - x_lim(2,1) ) / (ny)
        dy = ( ymax - ymin ) / (ny)
        DO i = 1 , ny
            ! y_cells_vec(i) = x_lim(2,1) + dy*(i-1)
            y_cells_vec(i) = ymin + dy*(i-1)
        END DO
    ELSE
        ! xmax = x_lim(1,2)
        alpha_x = -LOG(1./dx_factor)/xmax
        n_inf = 1/(dx_0*alpha_x)
        nmax = INT(n_inf*(1-EXP(-alpha_x*xmax)))
        nx = nmax+1
        ! ny = 1

        ALLOCATE( i_range(nx) )
        ALLOCATE( x_cells_vec(nx) )
        ALLOCATE( y_cells_vec(ny) )


        i_range = (/ (i,i=0,nmax) /)
        x_cells_vec = -1/alpha_x*LOG(1-i_range/n_inf)
        y_cells_vec = 0
        ! n_cells_vec(1) = nmax+1
        ! n_cells_vec(2) = 1

    END IF



    n_cells=nx*ny
    ALLOCATE(dx_cells_vec(nx))
    ALLOCATE(dy_cells_vec(ny))
    dx_cells_vec(1:(nx-1)) = x_cells_vec(2:nx) - x_cells_vec(1:(nx-1))
    dx_cells_vec(nx) = x_lim(1,2) - x_cells_vec(nx)
    dy_cells_vec = 1
    ALLOCATE(Vc(nx,ny))
    DO i = 1,nx
        DO j = 1,ny
            Vc(i,j) = dx_cells_vec(i)*dy_cells_vec(j)
        END DO
    END DO

    V_total = ( x_lim(1,2)-x_lim(1,1) ) * ( x_lim(2,2)-x_lim(2,1) )                 ! total simulation volume
    
    ! IF (include_source .EQV. .true.) THEN
    !     N_all = 0
    !     N_simulated = 0
    ! ELSE
    !     N_all  = INT(V_total*n/Fn)
    !     N_simulated = N_all
    ! ENDIF

    ! time step parameters ---------------------------------------------------------
    nt = INT(tmax/dt)+1
    ALLOCATE(t_vec(nt))
    DO i = 1,nt
        t_vec(i) = (i-1)*dt
    END DO

    CALL CPU_TIME(t0)
    t_final = 0
    t_test=0
    t_collisions = 0
    t_BC = 0
    t_loop = 0
    n_saved = 0

    ! source/reservoir parameters --------------------------------------------------
    ws = 10*vth*dt                                                                  ! width of source cell
    ts = 5.d-4                                                                      ! pulse width of source valve opening
    
    hs = inlet_height
    Vs = ws*hs
    Num_s = INT(ns*Vs/Fn)                                                           ! Number of source particles in the reservoir

    ! Set up wall geometry ---------------------------------------------------------
    ! draw boundaries with vertical and horizontal lines (each row is endpoints of wall: (x1,y1,x2,y2))
    ! include vertical walls at inlet/outlet as first two rows
    ALLOCATE(x_walls(4,2))
    num_walls = 4
    x_walls(:,1) = (/ xmin,ymin,xmin,ymax /)                                        ! left vertical wall
    x_walls(:,2) = (/ xmax,ymin,xmax,ymax /)                                        ! right vertical wall

    ! bottom/top walls, gun geometry?

    nw = 2


    ! set up vectors/IC's ----------------------------------------------------------

    ! calculate expected total number of particles in simulation ()
    IF (include_source .EQV. .true.) THEN
        N_all = 0
        N_simulated = N_all
        v_avg = SQRT(8*k_b*T_g/(Pi*m_g))
        N_expected = INT( ns*v_avg*MIN(tmax,ts)/(4*Fn) )
        N_array = INT(N_expected*1.25)
    ELSE
        N_all  = INT(V_total*n/Fn)
        N_simulated = N_all
        N_expected = N_all
        N_array = N_expected
    END IF
    Npc_max = N_array


    ALLOCATE(x_vec(N_array,ndim))
    ALLOCATE(v_vec(N_array,3))
    ALLOCATE(i_cell_vec(N_array,2))
    ALLOCATE(i_cell_vec_prev(N_array,2))
    
    ALLOCATE(x_vec_prev(N_array,ndim))
    ALLOCATE(v_vec_prev(N_array,3))

    ALLOCATE(x_vec_unsorted(N_array,ndim))
    ALLOCATE(v_vec_unsorted(N_array,3))
    ALLOCATE(i_cell_vec_unsorted(N_array,2))

    ALLOCATE(reflected_out(N_array))
    ALLOCATE(reflected_in(N_array))
    ALLOCATE(removed_from_sim(N_array))
    ALLOCATE(removed_from_sim_unsorted(N_array))

    ALLOCATE(in_column(N_array))
    ALLOCATE(in_cell(N_array))
    ALLOCATE(i_counting(N_array))
    ALLOCATE(i_cur(N_array))

    ALLOCATE(xs_vec(Num_s,ndim))
    ALLOCATE(xs_vec_prev(Num_s,ndim))
    ALLOCATE(vs_vec(Num_s,3))
    ALLOCATE(vs_vec_prev(Num_s,3))
    ALLOCATE(counter_vec(Num_s))
    ALLOCATE(entered_sim(Num_s))

    x_vec(:,:) = 0
    v_vec(:,:) = 0
    i_cell_vec(:,:) = 0
    i_cell_vec_prev(:,:) = 0

    x_vec_prev(:,:) = 0
    v_vec_prev(:,:) = 0

    reflected_out(:) = .false.
    reflected_in(:) = .false.
    removed_from_sim(:) = .false.
    removed_from_sim_unsorted(:) = .false.

    i_counting = (/ (i,i=1,N_array) /)
    in_column(:) = .false.
    in_cell(:) = .false.
    i_cur(:) = 0

    xs_vec(:,:) = 0
    vs_vec(:,:) = 0
    x_vec_prev(:,:) = 0
    v_vec_prev(:,:) = 0

    counter_vec(:) = 0
    entered_sim(:) = .false.



    

    IF (restart_simulation .EQV. .false.) THEN
        IF (N_all >= 0) THEN
            CALL RANDOM_NUMBER(x_vec)
            x_vec(:,1) = x_vec(:,1)*(xmax-xmin)+xmin
            x_vec(:,2) = x_vec(:,2)*(ymax-ymin)+ymin
            CALL RANDN(v_vec(:,1),N_all)
            CALL RANDN(v_vec(:,2),N_all)
            CALL RANDN(v_vec(:,3),N_all)
            v_vec = v_vec*vth
            ! x_vec(:,2) = 0.5
        END IF

        ALLOCATE(N_total(nt))
        ALLOCATE(N_candidate_pairs_total(nt))
        ALLOCATE(N_accepted_pairs_total(nt))
        ALLOCATE(N_collisions_total(nt))
        ALLOCATE(N_added_total(nt))
        ALLOCATE(Npc_slice(nx,ny))
        ALLOCATE(starting_index(nx,ny))     ! this will break if you go back to 2D

        ALLOCATE(Npc_added(nx,ny))
        ALLOCATE(ncp_remainder(nx,ny))
        ALLOCATE(vr_max(nx,ny))
        
        N_total(:) = 0
        N_candidate_pairs_total(:) = 0
        N_accepted_pairs_total(:) = 0
        N_collisions_total(:) = 0
        N_added_total(:) = 0
        Npc_slice(:,:) = 0
        starting_index(:,:) = 0

        Npc_added(:,:) = 0
        ncp_remainder(:,:) = 0
        vr_max(:,:) = vr_max_0
    ELSE
        ! Add stuff for restarting here (need to be able to load data)
    END IF

    CALL UPDATE_CELL_INDEX

END SUBROUTINE INITIALIZE




