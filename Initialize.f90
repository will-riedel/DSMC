SUBROUTINE INITIALIZE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    REAL(8)::dx,dy
    INTEGER::i,j,current_sum
    ! This suboutine initializes various matrices/parameters in the simulation

    CALL CPU_TIME(t0_init)



    IF (restart_simulation .EQV. .false.) THEN
        it_restart = -1
        ! it_restart = 0
        ! need to set dir_cur differently, it's currently set to the restart_directory input
    END IF


    ! simulation geometry and gas properties -----------------------------------------    

    IF (gas_type(1:8) == "NITROGEN") THEN
        m_g = 28*m_p
        d_g = dN2
        WRITE(*,*) "Note: using Nitrogen"
    ELSE
        m_g = 2*m_p
        d_g = dH2
        WRITE(*,*) "Note: using Hydrogen"
    END IF

    ! m_g = 2*m_p
    ! d_g = dH2
    ! WRITE(*,*) "Note: using Hydrogen"
    ! m_g = 28*m_p
    ! d_g = dN2
    ! WRITE(*,*) "Note: using Nitrogen"

    vth = SQRT(k_b*T_g/m_g)
    c_s = Pi*d_g**2
    vr_max_0 = 5*vth
    ! x_lim = transpose(reshape( (/ 0.d0,inlet_length+gun_length+outlet_length,   0.d0,outlet_height /) , (/2,2/) ))
    ! x_lim = transpose(reshape( (/ 0.d0,0.61d0,   0.d0,.05d0 /) , (/2,2/) ))  ! full-gun geometry, I think
    ! x_lim = transpose(reshape( (/ 0.d0,0.61d0,   0.d0,.025d0 /) , (/2,2/) ))  ! half-gun geometry, I think
    ! x_lim = transpose(reshape( (/ 0.d0,0.8d0,   0.d0,0.2d0 /) , (/2,2/) ))      ! split problem
    
    !x_lim = transpose(reshape( (/ 0.d0,0.2761d0,   0.d0,.025d0 /) , (/2,2/) ))  ! updated gun-geometry
    x_lim = transpose(reshape( (/ 0.d0,0.2761d0+.045,   0.d0,.025d0 /) , (/2,2/) ))  ! updated gun-geometry extended to target

    xmin = x_lim(1,1)
    xmax = x_lim(1,2)
    ymin = x_lim(2,1)
    ymax = x_lim(2,2)    
    ymid = (ymin+ymax)/2

    IF (geometry_type(1:11) == "CYLINDRICAL") THEN
        RWF = RWF_input
    ELSE
        RWF = 1
    END IF

    ! source_type = "DOWNSTREAM"  ! "NORMAL", "2BEAM", "DOWNSTREAM"
    ! x_grid_type = "EXP_1"  ! "HOMOG", "EXP_1"
    ! y_grid_type = "EXP_1"  ! "HOMOG", "EXP_1", "EXP_2" , "SPLIT"
    ! iniital_distribution = "HOMOG"  ! "EMPTY", "HOMOG", "SPLIT"


    ! Set up x-grid -------------------------------------------------------------------

    ! x_split = 1.d10
    x_split = x_inlet

    IF (x_grid_type(1:5) == "HOMOG") THEN
        ! set up equally spaced grid points
        ALLOCATE( x_cells_vec(nx) )
        dx = ( xmax - xmin ) / (nx)
        DO i = 1 , nx
            x_cells_vec(i) = xmin + dx*(i-1)
        END DO

        x_split = 1.d10

    ELSE IF (x_grid_type(1:5) == "EXP_1") THEN
        ! ! find x_cells
        ! alpha_x = -LOG(1./dx_factor)/xmax
        ! n_inf = 1/(dx_0*alpha_x)
        ! nmax = INT(n_inf*(1-EXP(-alpha_x*xmax)))
        ! nx = nmax+1

        ! ALLOCATE( i_range_x(nx) )
        ! ALLOCATE( x_cells_vec(nx) )

        ! i_range_x = (/ (i,i=0,nmax) /)
        ! x_cells_vec = -1/alpha_x*LOG(1-i_range_x/n_inf)


        ! this includes an extra homogeneous portion for the high-density bottlececk region
        ! if you set that region to a large dx, then it disappears and should result in a normal 1-directional exponential
        alpha_x = -LOG(1./dx_factor)/(xmax-x_inlet)
        n_inf = 1./(dx_0*alpha_x)

        nmax_left = INT(FLOOR( (x_inlet-xmin)/dx_inlet ))
        nmax_right = INT( n_inf*(1.-EXP(-alpha_x*(xmax-x_inlet))) )
        nmax = INT( nmax_right + nmax_left)
        nx = nmax+1

        ALLOCATE( i_range_x(nx) )
        ALLOCATE( x_cells_vec(nx) )

        i_range_x = (/ (i,i=0,nmax) /)
        DO j = 1,nmax_left+1
            x_cells_vec(j) = xmin + dx_inlet*(j-1)
        END DO
        DO j = nmax_left+2,nx
            x_cells_vec(j) = -1/alpha_x*LOG(1-(i_range_x(j)-nmax_left)/n_inf) + x_inlet
        END DO

        ! WRITE(*,*) "x_cells_vec=",x_cells_vec

        x_split = 1.d10

    ELSE IF (x_grid_type(1:5) == "EXP_2") THEN
        WRITE(*,*) "--- symmetric x-grid not implemented yet ---"
        STOP
    ELSE
        WRITE(*,*) "--- Can't identify x-grid type ---"
        STOP
    END IF

    ! Set up y-grid -------------------------------------------------------------------
    ny_b = INT(ny*(ns_b/ns))
    
    IF (y_grid_type(1:5) == "HOMOG") THEN
        ! set up equally spaced grid points
        ALLOCATE( y_cells_vec(ny) )
        ALLOCATE( y_cells_vec_b(ny_b) )
        dy = ( ymax - ymin ) / (ny)
        DO i = 1 , ny
            y_cells_vec(i) = ymin + dy*(i-1)
        END DO


    ELSE IF (y_grid_type(1:5) == "SPLIT") THEN
        ! set up two sets of equally spaced grid points, one to use on each side of the split

        ALLOCATE( y_cells_vec(ny) )
        dy = ( ymax - ymin ) / (ny)
        DO i = 1 , ny
            y_cells_vec(i) = ymin + dy*(i-1)
        END DO

        ALLOCATE( y_cells_vec_b(ny_b) )
        dy = ( ymax - ymin ) / (ny_b)
        DO i = 1 , ny_b
            y_cells_vec_b(i) = ymin + dy*(i-1)
        END DO

        ! x_split = (xmax-xmin)/2 + xmin
        ! x_split = x_inlet

    ELSE IF (y_grid_type(1:5) == "EXP_1") THEN
        ! find y_cells 
        alpha_y = -LOG(1./dy_factor)/(ymax-ymin)
        n_inf = 1/(dy_0*alpha_y)
        nmax = INT(n_inf*(1-EXP(-alpha_y*(ymax-ymin) )))
        ny = nmax+1

        ALLOCATE( i_range_y(ny) )
        ALLOCATE( y_cells_vec(ny) )
        ALLOCATE( y_cells_vec_b(ny_b) )


        i_range_y = (/ (i,i=0,nmax) /)
        y_cells_vec = -1/alpha_y*LOG(1-i_range_y/n_inf)

    ELSE IF (y_grid_type(1:5) == "EXP_2") THEN
        ! find y_cells 
        IF (dy_0 < (ymax-ymin) ) THEN
            alpha_y = -LOG(1./dy_factor)/(ymax-ymid)
            n_inf = 1./(dy_0*alpha_y)
            nmax = INT(n_inf*( 1.-EXP( -alpha_y*(ymax-ymid) ) ))
            ny = 2*(nmax+1)

            ALLOCATE( i_range_y(nmax+1) )
            ALLOCATE( y_cells_half(nmax+1) )
            ALLOCATE( y_cells_vec(ny) )
            ALLOCATE( y_cells_vec_b(ny_b) )


            i_range_y = (/ (i,i=0,nmax) /)
            y_cells_half = -1/alpha_y*LOG(1.-i_range_y/n_inf)

            y_cells_vec(1) = ymin
            DO i = 2,(nmax+1)
                y_cells_vec(i) = ymid - y_cells_half( (nmax+1)-(i-2) )
            END DO
            DO i = (nmax+2),ny
                y_cells_vec(i) = ymid + y_cells_half( i-(nmax+2)+1 )
            END DO
        
        ELSE
            ny = 1
            ALLOCATE( y_cells_vec(ny) )
            y_cells_vec = 0
            ALLOCATE( y_cells_vec_b(ny_b) )

        END IF
    ELSE
        WRITE(*,*) "--- Can't identify y-grid type ---"
        STOP
    END IF




    ! set the minimum cells to start collisions
    cx_min_collisions = 1
    cy_min_collisions = 1


    n_cells=nx*ny
    ALLOCATE(dx_cells_vec(nx))
    ALLOCATE(dy_cells_vec(ny))
    ALLOCATE(dy_cells_vec_b(ny_b))
    

    dx_cells_vec(1:(nx-1)) = x_cells_vec(2:nx) - x_cells_vec(1:(nx-1))
    dx_cells_vec(nx) = xmax - x_cells_vec(nx)
    IF (ny == 1) THEN
        ! dy_cells_vec = 1
        dy_cells_vec = (ymax-ymin)
        dy_cells_vec_b = (ymax-ymin)
    ELSE
        dy_cells_vec(1:(ny-1)) = y_cells_vec(2:ny) - y_cells_vec(1:(ny-1))
        dy_cells_vec(ny) = ymax - y_cells_vec(ny)
        dy_cells_vec_b(1:(ny_b-1)) = y_cells_vec_b(2:ny_b) - y_cells_vec_b(1:(ny_b-1))
        dy_cells_vec_b(ny_b) = ymax - y_cells_vec_b(ny_b)
    END IF
    ALLOCATE(Vc(nx,ny))
    Vc(:,:) = 1

    DO cx = 1,nx
        IF (y_grid_type(1:5) == "SPLIT") THEN
            IF (x_cells_vec(cx) < x_split) THEN
                DO cy = 1,ny
                    Vc(cx,cy) = dx_cells_vec(cx)*dy_cells_vec(cy)
                END DO
            ELSE 
                DO cy = 1,ny_b
                    Vc(cx,cy) = dx_cells_vec(cx)*dy_cells_vec_b(cy)
                END DO
            END IF

        ELSE
            DO cy = 1,ny
                IF (geometry_type(1:11) == "CYLINDRICAL") THEN   
                    Vc(cx,cy) = dx_cells_vec(cx)*dy_cells_vec(cy)*(y_cells_vec(cy)+dy/2)*2*Pi       ! assumes dtheta = 2*pi (so r*d_theta = 2*pi*r)
                ELSE
                    Vc(cx,cy) = dx_cells_vec(cx)*dy_cells_vec(cy)       ! assumes dz = 1
                END if
            END DO
        END IF
    END DO

    ! V_total = ( xmax-xmin ) * ( ymax-ymin )                 ! total simulation volume
    V_total = SUM(Vc)      

    IF (geometry_type(1:11) == "CYLINDRICAL") THEN
        ALLOCATE(WF_cell_vec(nx,ny))
        DO i = 1,nx
            DO j = 1,ny
                WF_cell_vec(i,j) = 1 + RWF*(y_cells_vec(j)+dy/2)/ymax
            END DO
        END DO
    END IF                                 ! total simulation volume

    ! time step parameters ---------------------------------------------------------
    nt = INT(tmax/dt)+1
    ALLOCATE(t_vec(nt))
    DO i = 1,nt
        t_vec(i) = (i-1)*dt
    END DO

    CALL CPU_TIME(t0_total)
    t_total = 0
    t_index = 0
    t_collisions = 0
    t_BC = 0
    t_loop = 0
    n_saved = 0

    ! source/reservoir parameters --------------------------------------------------
    
    ! ts = 5.d-4                                                                      ! pulse width of source valve opening
    ! ts = 1.d10                                                                    ! pulse width of source valve opening


    ! IF (include_two_beams .EQV. .true.) THEN
    IF (source_type(1:5) == "2BEAM") THEN
        ! ws = v_beam*dt
        ws = MAX(v_beam*dt,10*vth*dt)
    ELSE
        ws = 10*vth*dt                                                                  ! width of source cell
    END IF

    IF ( (source_type(1:6) == "NORMAL") .OR. (source_type(1:5) == "2BEAM") ) THEN
        hs = y_inlet(2)-y_inlet(1)

        xs_min = xmin - ws
        xs_max = xmin
        xs2_min = xmax
        xs2_max = xmax + ws
        ! ys_min = ymid-hs/2
        ! ys_max = ymid+hs/2
        ys_min = y_inlet(1)
        ys_max = y_inlet(2)

    ELSE IF (source_type(1:10) == "DOWNSTREAM") THEN
        ! angled boundary

        x_source_corners(1,1) = .0235
        x_source_corners(1,2) = .0315-.025
        x_source_corners(2,1) = .0216
        x_source_corners(2,2) = .0334-.025
        hs = SQRT( (x_source_corners(2,1)-x_source_corners(1,1))**2 + &
                        (x_source_corners(2,2)-x_source_corners(1,2))**2 )

        xs_min = x_source_corners(1,1) - ws/SQRT(2.)
        xs_max = x_source_corners(1,1)

        b_source_A = x_source_corners(1,2) - 1*x_source_corners(1,1)
        b_source_B = x_source_corners(2,2) - 1*x_source_corners(2,1)
        b_source_barrier = x_source_corners(1,2) - (-1)*x_source_corners(1,1)

        Theta_source = Pi/4.    

        ! WRITE(*,*) "ws,hs,Vs,Num_s_exact,Num_s=",ws,hs,Vs,Num_s_exact,Num_s
        ! WRITE(*,*) "xs_min,xs_max=",xs_min,xs_max
        ! WRITE(*,*) "bA,bB,b_barrier=",b_source_A,b_source_B,b_source_barrier
        ! WRITE(*,*) "x_source_corners=",x_source_corners

    END IF

    IF (geometry_type(1:9) == "CARTESIAN") THEN
        Vs = ws*hs
        Num_s_exact = ns*Vs/Fn
        Num_s_exact_b = ns_b*Vs/Fn
    ELSE IF (geometry_type(1:11) == "CYLINDRICAL") THEN
        ALLOCATE(Vs_cell_vec(nx,ny))
        ! Vs_cell_vec = dy*2*Pi*(y_cells_vec+dy/2)*ws
        ! Num_s_exact = SUM(ns*Vs_cell_vec/(Fn*WF_cell_vec))
        ! Num_s_exact_b = SUM(ns_b*Vs_cell_vec/(Fn*WF_cell_vec))
        WF_avg = 1 + RWF*(.5*(ys_min + ys_max)/ymax)
        Vs = Pi*(ys_max**2 - ys_min**2)*ws
        Num_s_exact = ns*Vs/(Fn*WF_avg)
        Num_s_exact_b = ns_b*Vs/(Fn*WF_avg)
    ELSE
        WRITE(*,*) "Error: can't parse geometry type"
    END IF

    Num_s = CEILING(Num_s_exact)
    Num_s_frac = Num_s_exact - FLOOR(Num_s_exact)
    Num_s_b = CEILING(Num_s_exact_b)
    Num_s_frac_b = Num_s_exact_b - FLOOR(Num_s_exact_b)



    ! Set up wall geometry ---------------------------------------------------------
    ! draw boundaries with vertical and horizontal lines (each row is endpoints of wall: (x1,y1,x2,y2))
    ! include vertical walls at inlet/outlet as first two rows
    
    x_walls(:,1) = (/ xmin,ymin,xmin,ymax /)                                        ! left vertical wall
    x_walls(:,2) = (/ xmax,ymin,xmax,ymax /)                                        ! right vertical wall
    x_walls(:,3) = (/ xmin,ymin,xmax,ymin /)                                        ! bottom horizontal wall
    x_walls(:,4) = (/ xmin,ymax,xmax,ymax /)                                        ! top horizontal wall

    ! x_walls(:,5) = (/ 0.,0.,0.3,0.5 /)                                        ! top horizontal wall
    ! x_walls(:,5) = (/ 0.,0.,0.5,0.5 /)                                        ! top horizontal wall
    ! x_walls(:,5) = (/ .15,0.,.15,1. /)                                        ! top horizontal wall
    ! x_walls(:,5) = (/ .1,0.,0.1,.05 /)  
    ! x_walls(:,6) = (/ .025,.05,.05,.025 /)  

    ! num_walls = 6

    ALLOCATE(wall_angle_vec(num_walls))
    DO i = 1,num_walls
        xw1 = x_walls(1,i)
        yw1 = x_walls(2,i)
        xw2 = x_walls(3,i)
        yw2 = x_walls(4,i)
        IF (xw2==xw1) THEN
            wall_angle_vec(i) = Pi/2.
        ELSE
            wall_angle_vec(i) = ATAN((yw2-yw1)/(xw2-xw1))
            IF (xw2 < xw1) THEN
                wall_angle_vec(i) = wall_angle_vec(i) + Pi
            END IF
        END IF
    END DO

    ALLOCATE(i_cell_lim_x(num_walls,2))
    ALLOCATE(i_cell_lim_y(num_walls,2))
    i_cell_lim_x(:,1) = 1
    i_cell_lim_x(:,2) = nx
    i_cell_lim_y(:,1) = 1
    i_cell_lim_y(:,2) = ny

    ! find cell range for each wall

    CALL FIND_WALL_CELLS


    N_specular = 0
    N_diffuse = 0

    ! set up vectors/IC's ----------------------------------------------------------

    ! calculate expected total number of particles in simulation ()

    
    

    IF (initial_distribution(1:5) == "EMPTY") THEN
        N_simulated = 0
    ELSE IF (initial_distribution(1:5) == "HOMOG") THEN
        IF (geometry_type(1:9) == "CARTESIAN") THEN
            N_init_a = INT(V_total*n/Fn)
        ELSE IF (geometry_type(1:11) == "CYLINDRICAL") THEN
            ! N_init_a = INT(V_total*n/(Fn*RWF))
            N_init_a = INT(SUM(Vc*n/(Fn*WF_cell_vec)))

        ELSE
            WRITE(*,*) "Error interpreting geometry type"
            STOP
        END IF

        N_simulated =  N_init_a
    ELSE IF (initial_distribution(1:5) == "SPLIT") THEN
        IF (geometry_type(1:9) == "CARTESIAN") THEN
            N_init_a = INT(.5*V_total*n/Fn)
            N_init_b = INT(.5*V_total*n_b/Fn)
        ELSE IF (geometry_type(1:11) == "CYLINDRICAL") THEN
            N_init_a = INT(.5*V_total*n/(Fn*RWF))
            N_init_b = INT(.5*V_total*n_b/(Fn*RWF))
        ELSE
            WRITE(*,*) "Error interpreting geometry type"
            STOP
        END IF
        
        N_simulated =  N_init_a + N_init_b
    ELSE
        WRITE(*,*) "Error interpreting initial distribution input"
        STOP
    END IF


    IF (include_source .EQV. .true.) THEN

        IF (source_type(1:5) == "2BEAM") THEN

            ! ! v_avg = v_beam*2 ! the factor of 2 is for both sides
            ! ! N_expected = 2*Num_s*nt
            ! ! N_expected = CEILING(2*Num_s_exact*nt)
            ! N_expected = CEILING((Num_s_exact+Num_s_exact_b)*nt)    ! this was used for when using beams, not reservoirs
            v_avg = SQRT(8*k_b*T_g/(Pi*m_g))
            N_expected = INT( ns*v_avg*MIN(tmax,ts)*(hs*1)/(4*Fn) ) + INT( ns_b*v_avg*MIN(tmax,ts)*(hs*1)/(4*Fn) )     ! hs*1 = cross-sectional area of inlet

        ELSE
            v_avg = SQRT(8*k_b*T_g/(Pi*m_g))
            IF (geometry_type(1:9) == "CARTESIAN") THEN    
                N_expected = INT( ns*v_avg*MIN(tmax,ts)*(hs*1)/(4*Fn) )     ! (nc/4)*(ts)*(A) , hs*1 = cross-sectional area of inlet
            ELSE IF (geometry_type(1:11) == "CYLINDRICAL") THEN    
                !N_expected = INT( ns*v_avg*MIN(tmax,ts)*(hs*1)/(4*Fn*RWF) )     ! hs*1 = cross-sectional area of inlet
                N_expected = INT( ns*v_avg*MIN(tmax,ts)*(Pi*(ys_max**2 - ys_min**2))/(4*Fn*WF_avg) )     ! (nc/4)*(ts)*(A)
                

                N_expected = N_expected / 2. ! getting a lot less
            END IF

        END IF

        N_array = INT(N_expected*1.5) + N_simulated
        !N_array = INT(N_expected*2.5) + N_simulated
        N_expected = N_expected + N_simulated

    ELSE
        N_expected = N_simulated
        N_array = N_expected+1

        IF (geometry_type(1:11) == "CYLINDRICAL") THEN    
            N_array = INT(N_array*1.5)
            ! N_array = INT(N_array*2)
        END IF

    END IF
    !IF (geometry_type(1:11) == "CYLINDRICAL") THEN    
    !    N_array = N_array*1.5    
    !END IF    
    Npc_max = N_array


    ! WRITE(*,*) "ns,Fn,ys_min,ys_max = ",ns,Fn,ys_min,ys_max
    ! WRITE(*,*) "Num_s = ",Num_s
    ! WRITE(*,*) "WF_avg, RWF = ",WF_avg, RWF
    ! WRITE(*,*) "v_avg,tmax,ts,ws = ",v_avg,tmax,ts,ws




    WRITE(*,*) "N_init_a,N_simulated,N_expected,N_array=",N_init_a,N_simulated,N_expected,N_array










    ! WRITE(*,*) "---",ws,hs,Num_s,Num_s_exact,Num_s_frac
    ! WRITE(*,*) y_inlet
    ! WRITE(*,*) source_type
    ! WRITE(*,*) v_avg,N_expected
    ! STOP

    ! WRITE(*,*) "N_expected,ws,dt = ",N_expected,ws,dt

    ALLOCATE(x_vec(N_array,ndim))
    ALLOCATE(v_vec(N_array,3))
    ALLOCATE(i_cell_vec(N_array,2))
    ALLOCATE(i_cell_vec_prev(N_array,2))
    
    ALLOCATE(x_vec_prev(N_array,ndim))
    ALLOCATE(v_vec_prev(N_array,3))

    ALLOCATE(x_vec_unsorted(N_array,ndim))
    ALLOCATE(v_vec_unsorted(N_array,3))
    ALLOCATE(i_cell_vec_unsorted(N_array,2))

    ALLOCATE(weight_factor_vec(N_array))
    ALLOCATE(weight_factor_vec_prev(N_array))
    ALLOCATE(weight_factor_vec_unsorted(N_array))
    ALLOCATE(reflected_out(N_array))
    ALLOCATE(reflected_in(N_array))
    ALLOCATE(removed_from_sim(N_array))
    ALLOCATE(removed_from_sim_unsorted(N_array))
    ALLOCATE(rn_vec(N_array))

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
    rn_vec(:) = 0

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
    ALLOCATE(NWF_escaped_total(nt))
    ALLOCATE(flux_upstream_total(nt))
    ALLOCATE(flux_downstream_total(nt))
    ALLOCATE(Npc_slice(nx,ny))
    ALLOCATE(starting_index(nx,ny))
    ALLOCATE(final_index(nx,ny))

    ALLOCATE(Npc_added(nx,ny))
    ! ALLOCATE(ncp_remainder(nx,ny))
    ALLOCATE(vr_max(nx,ny))
    
    N_total(:) = 0
    N_candidate_pairs_total(:) = 0
    N_accepted_pairs_total(:) = 0
    N_collisions_total(:) = 0
    N_added_total(:) = 0
    NWF_escaped_total(:) = 0
    flux_upstream_total(:) = 0
    flux_downstream_total(:) = 0
    Npc_slice(:,:) = 0
    starting_index(:,:) = 0
    final_index(:,:) = 0


    Npc_added(:,:) = 0
    ! ncp_remainder(:,:) = 0
    vr_max(:,:) = vr_max_0



    ALLOCATE(xr_vec(N_array,ndim))
    ALLOCATE(xr_vec_prev(N_array,ndim))
    ALLOCATE(xr_vec_new(N_array,ndim))
    ALLOCATE(vr_vec(N_array,3))
    ALLOCATE(vr_vec_prev(N_array,3))
    ALLOCATE(vr_vec_new(N_array,3))
    ALLOCATE(xr_walls(4,num_walls))
    ALLOCATE(collision_occured(N_array,num_walls))
    ALLOCATE(collision_occured_any(N_array))
    ALLOCATE(collision_dt(N_array,num_walls))
    ALLOCATE(min_collision_dt(N_array))
    ALLOCATE(x0(N_array))
    ALLOCATE(y0(N_array))
    ALLOCATE(xt(N_array))
    ALLOCATE(yt(N_array))
    ALLOCATE(xy0(N_array,ndim))
    ALLOCATE(xyt(N_array,ndim))
    ALLOCATE(m(N_array))
    ALLOCATE(b(N_array))
    ALLOCATE(xc(N_array))
    ALLOCATE(yc(N_array))
    ALLOCATE(dt_cross(N_array))
    ALLOCATE(first_collision(N_array))
    ALLOCATE(crossed(N_array))
    ALLOCATE(i_cross(N_array))
    ALLOCATE(i_first(N_array))

    xr_vec(:,:) = 0
    xr_vec_prev(:,:) = 0
    xr_vec_new(:,:) = 0
    vr_vec(:,:) = 0
    vr_vec_prev(:,:) = 0
    vr_vec_new(:,:) = 0
    !vn(:,:) = 0
    xr_walls(:,:) = 0
    collision_occured(:,:) = .false.
    collision_occured_any(:) = .false.
    collision_dt(:,:) = 0
    x0(:) = 0
    y0(:) = 0
    xt(:) = 0
    yt(:) = 0
    xy0(:,:) = 0
    xyt(:,:) = 0
    m(:) = 0
    b(:) = 0
    xc(:) = 0
    yc(:) = 0
    dt_cross(:) = 0
    first_collision(:) = .false.
    crossed(:) = .false.
    i_cross(:) = 0
    i_first(:) = 0

    N_entered = 0
    N_good_prob = 0
    N_bad_prob = 0

    IF (restart_simulation .EQV. .false.) THEN

        IF (initial_distribution(1:5) == "EMPTY") THEN
            ! no initialization
        ELSE IF (initial_distribution(1:5) == "HOMOG") THEN
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

        ELSE IF (initial_distribution(1:5) == "SPLIT") THEN
            IF (N_simulated > 2) THEN
                CALL RANDOM_NUMBER(x_vec)
                x_vec(1:N_init_a,1) = x_vec(1:N_init_a,1)*(x_split-xmin)+xmin
                x_vec((N_init_a+1):(N_init_a+N_init_b),1) = x_vec((N_init_a+1):(N_init_a+N_init_b),1)*(xmax-x_split)+x_split
                x_vec(:,2) = x_vec(:,2)*(ymax-ymin)+ymin
                CALL RANDN(v_vec(:,1),N_simulated)
                CALL RANDN(v_vec(:,2),N_simulated)
                CALL RANDN(v_vec(:,3),N_simulated)
                v_vec = v_vec*vth
                ! x_vec(:,2) = 0.5
            END IF

        ELSE 
            WRITE(*,*) "Error interpreting initial distribution input"
            STOP
        END IF


        ! clear out the current directory
        CALL SYSTEM( "mkdir " // dir_cur(1:dir_cur_length) )
        ! WRITE(*,*) "rm " // dir_cur(1:dir_cur_length) // "/*.txt"
        CALL SYSTEM( "rm " // dir_cur(1:dir_cur_length) // "/*.txt" )

        ! save input parameter file to data directory -------------------------------------------------
        filename = dir_cur(1:dir_cur_length)
        ! WRITE(*,*) "cp Input/Input_Parameters " // dir_cur(1:dir_cur_length) // "/Input_Parameters" 
        CALL SYSTEM( "cp Input/Input_Parameters " // dir_cur(1:dir_cur_length) // "/Input_Parameters" )
        
    ELSE
        CALL RESTART_PARAMETERS_READIN

        WRITE(*,*) "N_simulated=",N_simulated
        ! WRITE(*,*) "N_collisions_total=",N_simulated
        WRITE(*,*) "t_total=",t_total
        WRITE(*,*) "t_BC=",t_BC
        WRITE(*,*) "t_collisions=",t_collisions
        WRITE(*,*) "t_loop=",t_loop

        WRITE(*,*) "N_sim=",N_simulated
        ! WRITE(*,*) "N_total=",N_total(1:n_saved)

    END IF


    CALL CPU_TIME(t_temp)
    t_init = t_init + (t_temp-t0_init)


    IF (geometry_type(1:11) == "CYLINDRICAL") THEN    
        ! CALL UPDATE_WEIGHTS
        weight_factor_vec(:) = 1 + RWF*x_vec(:,2)/ymax
        weight_factor_vec_prev(:) = 1 + RWF*x_vec(:,2)/ymax



        

    ELSE
        weight_factor_vec(:) = 1
        weight_factor_vec_prev(:) = 1
    END IF

    CALL UPDATE_CELL_INDEX

END SUBROUTINE INITIALIZE







