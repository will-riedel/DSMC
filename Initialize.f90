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
        ! need to set dir_cur differently, it's currently set to the restart_directory input
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

    xs_min = xmin - ws
    xs_max = xmin
    ymid = (ymin+ymax)/2
    ys_min = ymid-hs/2
    ys_max = ymid+hs/2

    ! Set up wall geometry ---------------------------------------------------------
    ! draw boundaries with vertical and horizontal lines (each row is endpoints of wall: (x1,y1,x2,y2))
    ! include vertical walls at inlet/outlet as first two rows
    ALLOCATE(x_walls(4,2))
    ! num_walls = 4
    num_walls = 2
    x_walls(:,1) = (/ xmin,ymin,xmin,ymax /)                                        ! left vertical wall
    x_walls(:,2) = (/ xmax,ymin,xmax,ymax /)                                        ! right vertical wall
    ! x_walls(1,1) = xmin
    ! x_walls(2,1) = ymin
    ! x_walls(3,1) = xmin
    ! x_walls(4,1) = ymax

    ! x_walls(1,2) = xmax
    ! x_walls(2,2) = ymin
    ! x_walls(3,2) = xmax
    ! x_walls(4,2) = ymax

    ! bottom/top walls, gun geometry?

    nw = 2


    ! set up vectors/IC's ----------------------------------------------------------

    ! calculate expected total number of particles in simulation ()
    IF (include_source .EQV. .true.) THEN
        N_all = 0
        N_simulated = 0
        v_avg = SQRT(8*k_b*T_g/(Pi*m_g))
        N_expected = INT( ns*v_avg*MIN(tmax,ts)*(hs*1)/(4*Fn) )     ! hs*1 = cross-sectional area of inlet
        N_array = INT(N_expected*1.25)
    ELSE
        N_all  = INT(V_total*n/Fn)
        N_simulated = INT(V_total*n/Fn)
        N_expected = N_simulated
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
    ALLOCATE(vs_vec(Num_s,3))
    ALLOCATE(xs_vec_prev(Num_s,ndim))
    ALLOCATE(vs_vec_prev(Num_s,3))
    ALLOCATE(counter_vec(Num_s))
    ALLOCATE(entered_sim(Num_s))

    x_vec(:,:) = 0
    v_vec(:,:) = 0
    i_cell_vec(:,:) = 0
    i_cell_vec_prev(:,:) = 0

    x_vec_prev(:,:) = 0
    v_vec_prev(:,:) = 0

    x_vec_unsorted(:,:) = 0
    v_vec_unsorted(:,:) = 0
    i_cell_vec_unsorted(:,:) = 0

    reflected_out(:) = .false.
    reflected_in(:) = .false.
    removed_from_sim(:) = .false.
    removed_from_sim_unsorted(:) = .false.

    in_column(:) = .false.
    in_cell(:) = .false.
    i_counting = (/ (i,i=1,N_array) /)
    i_cur(:) = 0

    xs_vec(:,:) = 0
    vs_vec(:,:) = 0
    x_vec_prev(:,:) = 0
    v_vec_prev(:,:) = 0
    counter_vec(:) = 0
    entered_sim(:) = .false.



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


    

    IF (restart_simulation .EQV. .false.) THEN
        IF (N_simulated >= 0) THEN
            CALL RANDOM_NUMBER(x_vec)
            x_vec(:,1) = x_vec(:,1)*(xmax-xmin)+xmin
            x_vec(:,2) = x_vec(:,2)*(ymax-ymin)+ymin
            CALL RANDN(v_vec(:,1),N_simulated)
            CALL RANDN(v_vec(:,2),N_simulated)
            CALL RANDN(v_vec(:,3),N_simulated)
            v_vec = v_vec*vth
            ! x_vec(:,2) = 0.5
        END IF

        
    ELSE
        ! Add stuff for restarting here (need to be able to load data)

        ! it_restart = it_restart+1


        CALL RESTART_PARAMETERS_READIN
        ! WRITE(*,*) "N_total) = ",N_total
        ! WRITE(*,*) "N_candidate_pairs_total = ",N_candidate_pairs_total
        ! ! WRITE(*,*) "N_candidate_pairs_total = ",N_accepted_pairs_total
        ! WRITE(*,*) "ncp_remainder = ",ncp_remainder


        ! WRITE(*,*) "N_simulated = ",N_simulated
        ! ##########################################################################################
        ! ##########################################################################################
        ! ##########################################################################################
        ! ##########################################################################################
        ! Sometimes when this runs, N=1002, and it gets a backtrace error. Other times N=0 and it's fine.
        ! Not sure why either can happen

        ! ##########################################################################################
        ! ##########################################################################################
        ! ##########################################################################################
        ! ##########################################################################################
        ! ##########################################################################################
        ! ##########################################################################################



        ! for now, dir_cur = "Output/data"
        ! WRITE(filename,"('Output/data/a.txt')")
        WRITE(filename,"('Output/data/x_',I7.7,'.txt')") (it_restart)
        OPEN(UNIT=1,FILE=filename,STATUS='old', FORM='unformatted')! ,access='direct',recl=4,iostat=ok)
        READ(1) x_vec(1:N_simulated,:)
        CLOSE(1)

        WRITE(filename,"('Output/data/v_',I7.7,'.txt')") (it_restart)
        OPEN(UNIT=1,FILE=filename,STATUS='old', FORM='unformatted')! ,access='direct',recl=4,iostat=ok)
        READ(1) v_vec(1:N_simulated,:)
        CLOSE(1)

        WRITE(filename,"('Output/data/i_',I7.7,'.txt')") (it_restart)
        OPEN(UNIT=1,FILE=filename,STATUS='old', FORM='unformatted')! ,access='direct',recl=4,iostat=ok)
        READ(1) i_cell_vec(1:N_simulated,:)
        CLOSE(1)

        WRITE(filename,"('Output/data/Npc_',I7.7,'.txt')") (it_restart)
        OPEN(UNIT=1,FILE=filename,STATUS='old', FORM='unformatted')! ,access='direct',recl=4,iostat=ok)
        READ(1) Npc_slice
        CLOSE(1)



        


        

    END IF

    CALL UPDATE_CELL_INDEX

END SUBROUTINE INITIALIZE







