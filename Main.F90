MODULE CONTAIN
    IMPLICIT NONE
    SAVE
!-----------------------------------------------------------------------
!*******************PARALLEL VARIBLES**********************
!-----------------------------------------------------------------------
    ! Parallel Assembly/Solver Variables
    ! INTEGER:: Nparts, Neqn_Global
!-----------------------------------------------------------------------
!*******************DYNAMIC ARRAYS*************************
!-----------------------------------------------------------------------
    ! Relevant Allocatable Arrays
    ! REAL(8), ALLOCATABLE, DIMENSION(:):: x, y
    ! REAL(8), ALLOCATABLE, DIMENSION(:):: U, Material_Label, &
    !      Valboundary
    ! REAL(8), ALLOCATABLE, DIMENSION(:):: Valinitial
    ! INTEGER, ALLOCATABLE, DIMENSION(:):: nodeboundary_x, &
    !      nodeboundary_y, Dofboundary, Flux_Label
    ! INTEGER, ALLOCATABLE, DIMENSION(:):: Nodeinitial, Dofinitial
    ! INTEGER, ALLOCATABLE, DIMENSION(:, :):: Node_Type_Part_X, Node_Type_Part_Y
    ! INTEGER, ALLOCATABLE, DIMENSION(:, :):: Node_Type_Phi_X, Node_Type_Phi_Y
    ! INTEGER, ALLOCATABLE, DIMENSION(:, :):: Globaldof, Localdof
    ! REAL(8), ALLOCATABLE, DIMENSION(:, :):: n_i_plus, n_i_minus, n_i_minus_orig
    ! REAL(8), ALLOCATABLE, DIMENSION(:, :):: n_e_plus, n_e_minus, n_e_minus_orig
    ! REAL(8), ALLOCATABLE, DIMENSION(:, :):: Error_Offset_n_i, Error_Offset_n_e
    ! REAL(8), ALLOCATABLE, DIMENSION(:, :):: phi_plus
    ! REAL(8), ALLOCATABLE, DIMENSION(:, :):: Sigma_Charge
    ! REAL(8), ALLOCATABLE, DIMENSION(:, :, :):: K_i, K_e
!-----------------------------------------------------------------------
!*******************MISC. VARIBLES*************************
!-----------------------------------------------------------------------
    ! LOGICAL:: Petsc_Logical, Assembly_Check
    ! INTEGER:: Myrank, Mysize, Ierr, Iter, Fiter, Rf, Nonzeros_Per_Row
    ! REAL(8):: Clock_Start, Clock_Finish
    ! INTEGER:: Nnodes, Nnodes_x, Nnodes_y, Ndimensions, Nboundary, Ninitial
    ! INTEGER:: NFlux_Part, NFlux_Phi, Size
    ! INTEGER:: Ndof, Neqn, Print, Jacobian_Index, Max_Dof_of_Node, Runge_Kutta_Order
    ! REAL(8):: Initial_Time, Final_Time, Error_Offset, Scaling_Factor, Time_Error_Threshold
    ! REAL(8):: Delta_Time
    ! INTEGER:: Stencil_Width, Jacobian_Size




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
REAL(8), ALLOCATABLE, DIMENSION(:):: N_total,n_candidate_pairs_total,n_accepted_pairs_total,n_collisions_total,N_added_total
REAL(8), ALLOCATABLE, DIMENSION(:):: counter_vec
REAL(8), ALLOCATABLE, DIMENSION(:,:):: Vc, x_walls, vr_max, x_vec,v_vec, i_cell_vec, Npc_slice,ncp_remainder
REAL(8), ALLOCATABLE, DIMENSION(:,:):: x_vec_prev,v_vec_prev, xs_vec,vs_vec,xs_vec_prev,vs_vec_prev

!-----------------------------------------------------------------------
!*******************MISC. VARIBLES*************************
!-----------------------------------------------------------------------
    REAL(8):: m_g,d_g,vth,c_s,vr_max_0,v_avg, xmin,xmax,ymin,ymax,ymid, alpha_x,n_inf, V_total, t_collisions,t_BC,t_loop
    REAL(8):: ws,ts,hs,Vs,xs_min,xs_max,ys_min,ys_max,t, n_candidate_pairs_real
    INTEGER:: nmax, nx, ny, N_all, nt, n_saved, nw, N_expected, N_array, Num_s, N_entered, cx,cy
    INTEGER:: n_candidate_pairs,n_accepted_pairs,Npc_cur, num_walls
    REAL(8), DIMENSION(2,2):: x_lim
    INTEGER, DIMENSION(2):: n_cells_vec

    REAL(8):: a,b,c
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: a_mat,b_mat,c_mat
    REAL(8), ALLOCATABLE, DIMENSION(:):: a_vec,b_vec,c_vec
    
    ! INTEGER, ALLOCATABLE, DIMENSION(:):: i_range

    ! REAL(8)::

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













PROGRAM MAIN
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER:: Timestep,i
    REAL(8):: Tot_Time


    ! ALLOCATE(a_vec(10))
    ! ALLOCATE(b_vec(10))
    ! DO i = 1,10 
    !     a_vec(i) = i
    ! END DO

    ! ! ! b_vec(a_vec>5) = 1
    ! ! FORALL(i=1:10,a_vec(i)>5)
    ! !     b_vec(i)=1
    ! ! END FORALL

    ! WHERE (a_vec>5)
    !     b_vec = 1
    ! ELSEWHERE
    !     b_vec = 0
    ! END WHERE

    ! WRITE(*,*) "a_vec=",INT(a_vec)
    ! WRITE(*,*) "b_vec=",INT(b_vec)

!-----------------------------------------------------------------------
!*******************READ INPUT FILE************************
!-----------------------------------------------------------------------
    ! CALL INPUT_READIN
    CALL INPUT_PARAMETERS_READIN
!-----------------------------------------------------------------------
!*******************SET INITIAL CONDITIONS*****************
!-----------------------------------------------------------------------
    ! tot_time = initial_time
    ! Timestep = 0
    ! assembly_check = .true.
    CALL INITIALIZE
! !-----------------------------------------------------------------------
! !*******************MAIN LOOP******************************
! !-----------------------------------------------------------------------
!     DO
!        ! Increment the total time
!        tot_time = tot_time + delta_time
!        timestep = timestep + 1
!        ! Solve/Form Matrices Using PETSc
!        CALL PETSC
    DO i = it_restart+2,nt
    ! DO i = it_restart+2,50
        ! print current step
        t = t_vec(i)
        IF (MOD(i,dt_to_save) == 0) THEN
            WRITE(*,*) i,", t=",t,", N=",N_all,' ------------------------'
        END IF

        IF ( (t>ts) .and. (include_source .EQV. .true.) ) THEN
            close_inlet = .true.
            include_source = .false.
        END IF

        x_vec_prev = x_vec
        v_vec_prev = v_vec

        ! collisionless motion --------------------------------------------------------
        x_vec = x_vec +dt*v_vec(:,1:2)

        ! Boundary Condition implementation -------------------------------------------
        ! (removing exiting particles should be absorbed into BC implementation)
        IF (N_all > 0) THEN
            CALL SPECULAR_REFLECTION(x_vec,x_vec_prev,v_vec,v_vec_prev,N_all,x_walls,num_walls,xmin,xmax,0)
            ! CALL SPECULAR_REFLECTION('SIM_PARTICLES___')
        END IF

        ! Input from Source/Reservoir -------------------------------------------------
        IF (include_source .EQV. .true.) THEN
            CALL INITIALIZE_SOURCE
            N_added_total(i) = N_entered            
        END IF

        ! Update Cell Index -----------------------------------------------------------
        CALL UPDATE_CELL_INDEX

        ! Collisions ------------------------------------------------------------------
        ! Just do separate collisions subroutine here
        CALL RUN_COLLISIONS
        


        ! update and save current data ---------------------------------------------



    END DO

    WRITE(*,*) '------ done ------'

    STOP
END PROGRAM MAIN





SUBROUTINE linspace(z, l, k, n)
    IMPLICIT NONE
    !// Argument declarations
    REAL(8), DIMENSION(n)::z
    REAL(8)::l,k
    INTEGER::n
    INTEGER::i
    REAL(8)::d
    d = (k-l)/n
    z(1) = l
    DO i = 2, n-1
        z(i) = z(i-1) + d
    END DO
    z(1) = l
    z(n) = k
    RETURN
END SUBROUTINE linspace

SUBROUTINE LOG2INT(int_vec,log_vec,n)
    IMPLICIT NONE
    INTEGER,DIMENSION(n)::int_vec
    LOGICAL,DIMENSION(n)::log_vec
    INTEGER::n,i

    WHERE(log_vec .eqv. .true.)
        int_vec = 1
    ELSEWHERE
        int_vec = 0
    END WHERE
END SUBROUTINE LOG2INT


! generate random numbers with normal distribution (only in 1 for now)
SUBROUTINE RANDN(y_out,n)
! SUBROUTINE RANDN(y_out)
    IMPLICIT NONE
    ! uses polar form of Box-Muller Transformaton (described here: http://www.design.caltech.edu/erik/Misc/Gaussian.html)
    ! REAL (8), DIMENSION(n):: x1,x2,w,y1,y2,y_out,i
    REAL(8), DIMENSION(n):: w,y_out
    REAL(8):: x1,x2,y1
    INTEGER:: n,i
    ! REAL(8):: x1,x2,w,y_out

    DO i=1,n
        w(i)=1.d10
        DO WHILE (w(i) >= 1.0)
            CALL RANDOM_NUMBER(x1)
            CALL RANDOM_NUMBER(x2)
            x1 = 2*x1 - 1
            x2 = 2*x2 - 1
            w(i) = x1*x1 + x2*x2
        END DO

        w(i) = SQRT( (-2*LOG(w(i))) / w(i) )
        y_out(i) = x1*w(i)        ! could also have y_out2 = x2*w
    END DO


END SUBROUTINE RANDN





SUBROUTINE OUTPUT_TO_SCREEN
    USE CONTAIN
    IMPLICIT NONE
    ! This Suboutine Prints Out Relevant Domain Info
    WRITE(*,*)
    WRITE(*,*) ' ********   ANALYSIS INFORMATION ********'
    ! WRITE(*,*) ' processor rank =', myrank
    ! WRITE(*,*) ' number of nodes in x =', nnodes_x
    ! WRITE(*,*) ' number of nodes in y=', nnodes_y
    ! WRITE(*,*) ' total number of nodes =', nnodes
    ! WRITE(*,*) ' number of bcs   =', nboundary
    ! WRITE(*,*) ' number of ics   =', ninitial
    ! WRITE(*,*) ' number of flux bcs   =', nflux_part+nflux_phi
    ! WRITE(*,*) ' number of equations =', neqn
    ! WRITE(*,*) ' number of dofs =', ndof


    WRITE(*,*) 'n=',n
    WRITE(*,*) 'tmax=',tmax
    WRITE(*,*) 'include_source=',include_source
    WRITE(*,*) 'use_homogenous_grid=',use_homogenous_grid
    WRITE(*,*) 'dir_cur = ',dir_cur

END SUBROUTINE OUTPUT_TO_SCREEN







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
        ELSE IF (line == '*N_CELLLS_Y____') THEN
          ! number of cells in the y direction (if using homogeneous grid)
          READ(100,*) n_cells_y
          n_cells_vec(2) = n_cells_y
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
        ELSE IF (line == '*CLOSE_INLET___') THEN
          ! whether to put boundary at inlet (after any source is turned off)
          READ(100,*) close_inlet
        ELSE IF (line == '*CLOSE_OUTLET__') THEN
          ! whether to put boundary at the outlet
          READ(100,*) close_outlet
        ELSE IF (line == '*INCLUDE_GUN_BC') THEN
          ! include the boundaries of the gun geometry
          READ(100,*) include_gun_boundaries
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







SUBROUTINE INITIALIZE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    REAL(8)::dx,dy
    INTEGER::i,j
    ! This suboutine initializes various matrices/parameters in the simulation

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
    

    ! ALLOCATE( x_cells_vec(1) )
    ! ALLOCATE( y_cells_vec(1) )

    IF (use_homogenous_grid .EQV. .true.) THEN
        ! set up equally spaced grid points
        ALLOCATE( x_cells_vec(n_cells_vec(1)) )
        ALLOCATE( y_cells_vec(n_cells_vec(2)) )
        dx = ( x_lim(1,2) - x_lim(1,1) ) / n_cells_vec(1)
        DO i = 1 , n_cells_vec(1)
            x_cells_vec(i) = x_lim(1,1) + dx*(i-1)
        END DO
        dy = ( x_lim(2,2) - x_lim(2,1) ) / n_cells_vec(2)
        DO i = 1 , n_cells_vec(2)
            y_cells_vec(i) = x_lim(2,1) + dy*(i-1)
        END DO
    ELSE
        xmax = x_lim(1,2)
        alpha_x = -LOG(1/dx_factor)/xmax
        n_inf = 1/(dx_0*alpha_x)
        nmax = INT(n_inf*(1-EXP(-alpha_x*xmax)))
        ALLOCATE(i_range(nmax+1))
        CALL linspace(i_range,0.d0,REAL(nmax,8),nmax+1)
        ALLOCATE( x_cells_vec(nmax+1) )
        ALLOCATE( y_cells_vec(1) )
        x_cells_vec = -1/alpha_x*LOG(1-i_range/n_inf)
        ! WRITE(*,*) SHAPE(x_cells_vec)
        y_cells_vec = 0

        n_cells_vec(1) = nmax+1
        n_cells_vec(2) = 1
    END IF

    nx = n_cells_vec(1)
    ny = n_cells_vec(2)
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
    N_all  = INT(V_total*n/Fn)


    ! time step parameters ---------------------------------------------------------
    ! WRITE(*,*) INT(tmax/dt)+1
    nt = INT(tmax/dt)+1
    ALLOCATE(t_vec(nt))
    DO i = 1,nt
        t_vec(i) = i*dt
    END DO

    ! t0 = time.time() ##########
    t_collisions = 0
    t_BC = 0
    t_loop = 0
    n_saved = 0

    ! source/reservoir parameters --------------------------------------------------
    ws = 10*vth*dt                                                                  ! width of source cell
    ts = 5.d-4                                                                      ! pulse width of source valve opening
    
    hs = inlet_height
    Vs = ws*hs
    Num_s = INT(ns*Vs/Fn)
    WRITE(*,*) "Num_s =",Num_s 

    ! Set up wall geometry ---------------------------------------------------------
    ! draw boundaries with vertical and horizontal lines (each row is endpoints of wall: (x1,y1,x2,y2))
    ! include vertical walls at inlet/outlet as first two rows
    ALLOCATE(x_walls(4,2))
    num_walls = 4
    xmin = x_lim(1,1)
    xmax = x_lim(1,2)
    ymin = x_lim(2,1)
    ymax = x_lim(2,2)
    x_walls(:,1) = (/ xmin,ymin,xmin,ymax /)                                        ! left vertical wall
    x_walls(:,2) = (/ xmax,ymin,xmax,ymax /)                                        ! right vertical wall

    ! bottom/top walls, gun geometry?

    nw = 2


    ! set up vectors/IC's ----------------------------------------------------------

    ! calculate expected total number of particles in simulation ()
    IF (include_source .EQV. .true.) THEN
        N_all = 0
        v_avg = SQRT(8*k_b*T_g/(Pi*m_g))
        N_expected = INT( ns*v_avg*MIN(tmax,ts)/(4*Fn) )
        N_array = INT(N_expected*1.25)
    ELSE
        N_expected = N_all
        N_array = N_expected
    END IF

    ALLOCATE(x_vec(N_array,ndim))
    ALLOCATE(v_vec(N_array,3))
    ALLOCATE(i_cell_vec(N_array,2))
    ALLOCATE(x_vec_prev(N_array,ndim))
    ALLOCATE(v_vec_prev(N_array,3))

    ALLOCATE(xs_vec(Num_s,ndim))
    ALLOCATE(xs_vec_prev(Num_s,ndim))
    ALLOCATE(vs_vec(Num_s,3))
    ALLOCATE(vs_vec_prev(Num_s,3))
    ALLOCATE(counter_vec(Num_s))


    ALLOCATE(vr_max(nx,ny))
    vr_max(:,:) = vr_max_0

    IF (restart_simulation .EQV. .false.) THEN
        IF (N_all > -1) THEN
            ! need to allocate x_vec if there's nothing there yet?
            ! DO j = 1,ndim
            !     xmin = x_lim
            ! x_vec(:,1) = np.random.rand(N_all,)*(xmax-xmin)+xmin
            ! y_vec(:,2) = np.random.rand(N_all,)*(ymax-ymin)+ymin
            ! DO i = 1,N_all
            !     DO j = 1,3
            !         v_vec(i,j) = vth*RAND()
            CALL RANDOM_NUMBER(x_vec)
            x_vec(:,1) = x_vec(:,1)*(xmax-xmin)+xmin
            x_vec(:,2) = x_vec(:,2)*(ymax-ymin)+ymin
            ! CALL RANDOM_NUMBER(v_vec)   ! need to convert this to normal distribution #########
            CALL RANDN(v_vec(:,1),N_all)
            CALL RANDN(v_vec(:,2),N_all)
            CALL RANDN(v_vec(:,3),N_all)
            v_vec = v_vec*vth
            x_vec(:,2) = 0.5
        END IF

        ALLOCATE(N_total(nt))
        ALLOCATE(n_candidate_pairs_total(nt))
        ALLOCATE(n_accepted_pairs_total(nt))
        ALLOCATE(n_collisions_total(nt))
        ALLOCATE(N_added_total(nt))
        ! ALLOCATE(Npc_slice(nx,ny,1))
        ALLOCATE(Npc_slice(nx,ny))
        ALLOCATE(ncp_remainder(nx,ny))
    ELSE
        ! Add stuff for restarting here (need to be able to load data)
    END IF

    CALL UPDATE_CELL_INDEX




END SUBROUTINE INITIALIZE





SUBROUTINE UPDATE_CELL_INDEX
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    ! note: if (x,y) == (0,0), then just set index to -1000 or something (or 0, since indices start with 1 here)

    
    ! WRITE(*,*) "--- have not included UPDATE_CELL_INDEX function yet ---"

END SUBROUTINE UPDATE_CELL_INDEX





! SUBROUTINE SPECULAR_REFLECTION(string_in)
SUBROUTINE SPECULAR_REFLECTION(xr_vec,xr_vec_prev,vr_vec,vr_vec_prev,N_all,xr_walls,num_walls,xmin,xmax,counter)
    ! USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    ! CHARACTER(16):: string_in
    REAL(8),DIMENSION(N_all,ndim):: xr_vec,xr_vec_prev, x_coll, xr_vec_new
    REAL(8),DIMENSION(N_all,3):: vr_vec,vr_vec_prev, vr_vec_new
    REAL(8),DIMENSION(num_walls,4)::xr_walls
    LOGICAL,DIMENSION(N_all,num_walls):: collision_occured
    REAL(8),DIMENSION(N_all,num_walls):: collision_dt
    REAL(8),DIMENSION(N_all):: x0,y0,xt,yt,m,b,xc,yc, dt_cross
    LOGICAL,DIMENSION(N_all):: reflected_in, reflected_out, first_collision
    INTEGER,DIMENSION(N_all):: i_cross, i_first
    REAL(8):: xw1,xw2,yw1,yw2,xmin,xmax, temp
    INTEGER:: N_all,num_walls,counter, i

    ! IF (string_in == 'SOURCE_PARTICLES') THEN
    !     xr_vec = xs_vec
    !     xr_vec_prev = xs_vec_prev
    !     vr_vec = vs_vec
    !     vr_vec_prev = vs_vec_prev
    !     Num_r = Num_s
    !     xr_walls = x_walls(5:)
    ! ELSE IF (string_in == 'SIM_PARTICLES___') THEN
    !     xr_vec = x_vec
    !     xr_vec_prev = x_vec_prev
    !     vr_vec = v_vec
    !     vr_vec_prev = v_vec_prev
    !     Num_r = N_all
    !     xr_walls = x_walls
    ! ELSE
    !     WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    ! END IF

    xr_vec_new = xr_vec
    vr_vec_new = vr_vec
    ! i_refl_out/in?
    collision_occured(:,:) = .false.

    collision_dt = collision_dt + 1e10

    DO i = 1,num_walls
        xw1 = xr_walls(i,1)
        yw1 = xr_walls(i,2)
        xw2 = xr_walls(i,3)
        yw2 = xr_walls(i,4)

        x0 = xr_vec_prev(:,1)
        y0 = xr_vec_prev(:,2)
        xt = xr_vec(:,1)
        yt = xr_vec(:,2)
        m = vr_vec_prev(:,2)/vr_vec_prev(:,1)
        b = y0-m*x0

        ! if (xw1=xw2)
        ! vertical boundary
        IF (yw2 > yw1) THEN
            temp = yw1
            yw1 = yw2
            yw2 = temp
        END IF

        yc = m*xw1+b
        xc(:) = xw1

        ! i_cross = INT( ((x0<xw1) .and. (xw1<xt)) .or. ((x0>xw1).and.(xw1>xt)) )
        CALL LOG2INT( i_cross, ((x0<xw1).and.(xw1<xt)) .or. ((x0>xw1).and.(xw1>xt)) , N_all )
        dt_cross = (xc-x0)/vr_vec_prev(:,1)

        collision_occured(i_cross,i) = .true.
        collision_dt(i_cross,i) = dt_cross(i_cross)

        ! CALL LOG2INT( i_first, (collision_occured(:,i) .eqv. .true.) .and. (collision_dt(:,i) == MINVAL(collision_dt,2)) , N_all)
        ! xr_vec_new(i_first,1) = 2*xw1 - xr_vec(i_first,1)
        ! vr_vec_new(i_first,1) = -vr_vec(i_first,1)
        ! x_coll(i_first,1) = xc(i_first)
        ! x_coll(i_first,2) = yc(i_first)

        ! IF ((xw1==xmin).and.(xw2==xmin)) THEN
        !     reflected_in(i_first) = .true.
        ! ELSE IF (xw1==xmax).and.(xw2==xmax) THEN
        !     reflected_out(i_first) = .true.
        ! END IF


        first_collision = (collision_occured(:,i) .eqv. .true.) .and. (collision_dt(:,i) == MINVAL(collision_dt,2))
        WHERE (first_collision) 
            xr_vec_new(:,1) = 2*xw1 - xr_vec(:,1)
            vr_vec_new(:,1) = -vr_vec(:,1)
            x_coll(:,1) = xc
            x_coll(:,2) = yc
        ELSEWHERE
        END WHERE

        IF ((xw1==xmin).and.(xw2==xmin)) THEN
            WHERE(first_collision)
                reflected_in = .true.
            ELSEWHERE
            END WHERE
        ELSE IF ((xw1==xmax).and.(xw2==xmax)) THEN
            WHERE(first_collision)
                reflected_out = .true.
            ELSEWHERE
            END WHERE
        END IF

        ! other angles/wall-types go here

        xr_vec = xr_vec_new
        vr_vec = vr_vec_new

    END DO


    
    ! WRITE(*,*) "--- have not included SPECULAR_REFLECTION function yet ---"
END SUBROUTINE SPECULAR_REFLECTION





SUBROUTINE INITIALIZE_SOURCE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    xs_min = xmin - ws
    xs_max = xmin
    ymid = (ymin+ymax)/2
    ys_min = ymid-hs/2
    ys_max = ymid+hs/2
    ! CALL RAND(xs_vec)
    CALL RANDOM_NUMBER(xs_vec(:,1))
    xs_vec(:,1) = xs_vec(:,1)*(xs_max-xs_min) + xs_min
    xs_vec(:,2) = 0.5

    ! CALL RANDOM_NUMBER(vs_vec)   ! NEED TO IMPLEMENT NORMAL DISTRIUBUTION #########################
    CALL RANDN(vs_vec(:,1),Num_s) 
    CALL RANDN(vs_vec(:,2),Num_s) 
    CALL RANDN(vs_vec(:,3),Num_s) 
    vs_vec = vs_vec*vth


    ! keep incoming source particles within inlet
    xs_vec_prev = xs_vec
    vs_vec_prev = vs_vec
    xs_vec = xs_vec + dt*vs_vec(:,1:2)


    ! (not needed right now)
    ! CALL SPECULAR_REFLECTION('SOURCE_PARTICLES') ###############

    N_entered = COUNT(xs_vec(:,1) > 0)


    IF (N_entered > 0) THEN
        x_vec( (N_all+1):(N_all+1+N_entered) , : ) = xs_vec
        v_vec( (N_all+1):(N_all+1+N_entered) , : ) = vs_vec
        N_all = N_all + N_entered
    END IF

END SUBROUTINE INITIALIZE_SOURCE



SUBROUTINE RUN_COLLISIONS
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    INTEGER::i,j,k


    DO cx = 1,n_cells_vec(1)

        cy = 1

        ! collision processing
        n_candidate_pairs = 0
        n_accepted_pairs = 0
        Npc_cur = COUNT(i_cell_vec(:,1) == cx)
        Npc_slice(cx,cy) = Npc_cur

        IF (Npc_cur >= 2) THEN
            n_candidate_pairs_real = .5*Npc_cur*(Npc_cur-1)*Fn*c_s*vr_max(cx,cy)*dt/Vc(cx,cy)   ! number of candidate collision pairs
            n_candidate_pairs_real = n_candidate_pairs_real + ncp_remainder(cx,cy)
            ncp_remainder(cx,cy) = MOD(n_candidate_pairs_real,1.0)
            n_candidate_pairs = INT(n_candidate_pairs_real)


            FORALL(i = 1:N_all , i_cell_vec(i,1) == cx)
                    ! ######## this FORALL/WHERE stuff isn't really working out. Can only have assignemnt statements inside ##########

                ! DO k = 1,n_candidate_pairs


                ! END DO

            END FORALL


        END IF
    END DO


END SUBROUTINE RUN_COLLISIONS






