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
        ! it_restart = 0
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
    ymid = (ymin+ymax)/2

    ! nx = n_cells_vec(1)
    ! ny = n_cells_vec(2)
    ! ny = 1
    ! cy  = 1


    IF (use_homogenous_grid .EQV. .true.) THEN
        ! set up equally spaced grid points
        ALLOCATE( x_cells_vec(nx) )
        ALLOCATE( y_cells_vec(ny) )
        dx = ( xmax - xmin ) / (nx)
        DO i = 1 , nx
            x_cells_vec(i) = xmin + dx*(i-1)
        END DO
        dy = ( ymax - ymin ) / (ny)
        DO i = 1 , ny
            y_cells_vec(i) = ymin + dy*(i-1)
        END DO
    ELSE
        ! ! find x_cells
        ! alpha_x = -LOG(1./dx_factor)/xmax
        ! n_inf = 1/(dx_0*alpha_x)
        ! nmax = INT(n_inf*(1-EXP(-alpha_x*xmax)))
        ! nx = nmax+1

        ! ALLOCATE( i_range_x(nx) )
        ! ALLOCATE( x_cells_vec(nx) )

        ! i_range_x = (/ (i,i=0,nmax) /)
        ! x_cells_vec = -1/alpha_x*LOG(1-i_range_x/n_inf)


        alpha_x = -LOG(1./dx_factor)/(xmax-x_inlet)
        n_inf = 1./(dx_0*alpha_x)

        nmax_left = INT(FLOOR( (x_inlet-xmin)/dx_inlet ))
        nmax_right = INT( n_inf*(1.-EXP(-alpha_x*(xmax-x_inlet))) )
        nmax = INT( nmax_right + nmax_left)
        nx = nmax+1
        ! WRITE(*,*) "nmax_left,nmax_right,nmax,nx=",nmax_left,nmax_right,nmax,nx
        ! WRITE(*,*) "x_inlet,xmin,dx_inlet=",x_inlet,xmin,dx_inlet
        ! WRITE(*,*) "product = ",(x_inlet-xmin)/dx_inlet

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


        ! find y_cells
        IF (dy_0 < (ymax-ymin) ) THEN
            
            alpha_y = -LOG(1./dy_factor)/(ymax-ymid)
            n_inf = 1./(dy_0*alpha_y)
            nmax = INT(n_inf*( 1.-EXP( -alpha_y*(ymax-ymid) ) ))
            ny = 2*(nmax+1)

            ALLOCATE( i_range_y(nmax+1) )
            ALLOCATE( y_cells_half(nmax+1) )
            ALLOCATE( y_cells_vec(ny) )

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

        END IF

        ! WRITE(*,*) "ny,nmax,n_inf,dy_factor=",ny,nmax,n_inf,dy_factor
        ! WRITE(*,*) "i_range_y=",i_range_y
        ! WRITE(*,*) "y_cells_half=",y_cells_half
        ! WRITE(*,*) "ymin,ymax,ymid=",ymin,ymax,ymid
        ! WRITE(*,*) "ycv=",y_cells_vec

    END IF



    n_cells=nx*ny
    ALLOCATE(dx_cells_vec(nx))
    ALLOCATE(dy_cells_vec(ny))
    dx_cells_vec(1:(nx-1)) = x_cells_vec(2:nx) - x_cells_vec(1:(nx-1))
    dx_cells_vec(nx) = xmax - x_cells_vec(nx)
    IF (ny == 1) THEN
        ! dy_cells_vec = 1
        dy_cells_vec = (ymax-ymin)
    ELSE
        dy_cells_vec(1:(ny-1)) = y_cells_vec(2:ny) - y_cells_vec(1:(ny-1))
        dy_cells_vec(ny) = ymax - y_cells_vec(ny)
    END IF
    ALLOCATE(Vc(nx,ny))
    DO cx = 1,nx
        DO cy = 1,ny
            Vc(cx,cy) = dx_cells_vec(cx)*dy_cells_vec(cy)
        END DO
    END DO

    V_total = ( xmax-xmin ) * ( ymax-ymin )                 ! total simulation volume
    

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
    t_BC1 = 0
    t_BC2 = 0
    t_BC3 = 0
    t_BC4 = 0
    t_BC5 = 0
    t_loop = 0
    n_saved = 0

    ! source/reservoir parameters --------------------------------------------------
    
    ts = 5.d-4                                                                      ! pulse width of source valve opening

    IF (include_two_beams .EQV. .true.) THEN
        ws = v_beam*dt
    ELSE
        ws = 10*vth*dt                                                                  ! width of source cell
    END IF

    ! hs = inlet_height
    hs = y_inlet(2)-y_inlet(1)
    Vs = ws*hs
    ! Num_s = INT(ns*Vs/Fn)                                                           ! Number of source particles in the reservoir
    Num_s_exact = ns*Vs/Fn
    Num_s = CEILING(Num_s_exact)                                                           ! Number of source particles in the reservoir

    xs_min = xmin - ws
    xs_max = xmin
    xs2_min = xmax
    xs2_max = xmax + ws
    ! ys_min = ymid-hs/2
    ! ys_max = ymid+hs/2
    ys_min = y_inlet(1)
    ys_max = y_inlet(2)

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
    
    ! DO i = 1,num_walls
    !     WRITE(*,*) "i,cxmin,cxmax,cymin,cymax=", &
    !                 i,i_cell_lim_x(i,1),i_cell_lim_x(i,2),i_cell_lim_y(i,1),i_cell_lim_y(i,2)
    ! END DO



    N_specular = 0
    N_diffuse = 0

    ! set up vectors/IC's ----------------------------------------------------------


    if (include_two_beams .EQV. .true.) THEN
        include_source = .true.
    END IF

    ! calculate expected total number of particles in simulation ()
    IF (include_source .EQV. .true.) THEN
        ! N_all = 0
        N_simulated = 0
        IF (include_two_beams .EQV. .true.) THEN
            v_avg = v_beam*2 ! the factor of 2 is for both sides
            ! N_expected = 2*Num_s*nt
            N_expected = CEILING(2*Num_s_exact*nt)
        ELSE
            v_avg = SQRT(8*k_b*T_g/(Pi*m_g))
            N_expected = INT( ns*v_avg*MIN(tmax,ts)*(hs*1)/(4*Fn) )     ! hs*1 = cross-sectional area of inlet
        END IF
        N_array = INT(N_expected*1.25)
    ELSE
        ! N_all  = INT(V_total*n/Fn)
        N_simulated = INT(V_total*n/Fn)
        N_expected = N_simulated
        N_array = N_expected+1
    END IF
    Npc_max = N_array

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
    ALLOCATE(Npc_slice(nx,ny))
    ALLOCATE(starting_index(nx,ny))
    ALLOCATE(final_index(nx,ny))

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
    final_index(:,:) = 0

    Npc_added(:,:) = 0
    ncp_remainder(:,:) = 0
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


        ! clear out the current directory
        CALL SYSTEM( "mkdir " // dir_cur(1:dir_cur_length) )
        ! WRITE(*,*) "rm " // dir_cur(1:dir_cur_length) // "/*.txt"
        CALL SYSTEM( "rm " // dir_cur(1:dir_cur_length) // "/*.txt" )

        ! WRITE(filename,"('Output/data/x_',I7.7,'.txt')") (it_restart)
        ! WRITE(filename,"(dir_cur(1:dir_cur_length),'/x_',I7.7,'.txt')") (it_restart)
        ! WRITE(filename,"('/x_',I7.7,'.txt')") (25)
        ! WRITE(*,*) filename, "$$$"
        ! filename = dir_cur(1:dir_cur_length) // filename
        ! WRITE(*,*) filename, "$$$"

        ! WRITE(*,*) "got here"
        ! STOP
        
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

    CALL UPDATE_CELL_INDEX

END SUBROUTINE INITIALIZE







