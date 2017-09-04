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
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: Vc, x_walls, vr_max, x_vec,v_vec, i_cell_vec, i_cell_vec_prev, particles_in_cell,Npc_slice,ncp_remainder
    REAL(8), ALLOCATABLE, DIMENSION(:,:):: x_vec_prev,v_vec_prev, xs_vec,vs_vec,xs_vec_prev,vs_vec_prev
    LOGICAL,ALLOCATABLE, DIMENSION(:):: reflected_in, reflected_out, in_column, in_cell, removed_from_sim, entered_sim
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













PROGRAM MAIN
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER:: i
    CHARACTER(80)::filename
    CHARACTER(10)::str_file_num


    ! WRITE(*,*) "got here"

    ! ALLOCATE(a_vec(8))
    ! ALLOCATE(b_vec(8))
    ! ALLOCATE(c_vec(8))
    ! ALLOCATE(i_temp_vec(8))
    ! ALLOCATE(b_mat(8,2))
    ! DO i = 1,8
    !     a_vec(i) = i
    !     b_vec(i) = i
    ! END DO

    ! c_vec = a_vec*b_vec

    ! ! n_vec = (/ (i,i=1,10) /)

    ! ! ! b_vec(a_vec>5) = 1
    ! ! FORALL(i=1:10,a_vec(i)>5)
    ! !     b_vec(i)=1
    ! ! END FORALL

    ! ! WHERE (a_vec>5)
    ! !     b_vec = 1
    ! !     b_mat(:,1) = 1
    ! ! ELSEWHERE
    ! !     b_vec = 0
    ! !     b_mat(:,1) = 5
    ! ! END WHERE

    ! ! b_vec = a_vec((/i,2/))

    ! WRITE(*,*) "a_vec=",INT(a_vec)
    ! WRITE(*,*) "b_vec=",INT(b_vec)
    ! ! WRITE(*,*) "i_temp_vec=",INT(i_temp_vec)
    ! WRITE(*,*) "c_vec=",INT(c_vec)

    ! ! CALL RANDOM_NUMBER(a1)
    ! ! WRITE(*,*) a1
    ! ! ALLOCATE(a_vec(10000))
    ! ! CALL RANDN(a_vec,10000)

    ! ! INQUIRE(FILE="Output/data/mat.txt",EXIST=file_exists)
    ! ! IF (file_exists .EQV. .false.) THEN
    ! OPEN(UNIT=1,FILE="Output/data/mat.txt",FORM="UNFORMATTED")
    ! WRITE(1) a_vec
    ! CLOSE(1)
    ! ! ELSE
    ! !     WRITE(*,*) "file already exists!"
    ! ! END IF
    ! ! DEALLOCATE(a_vec)

!-----------------------------------------------------------------------
!*******************READ INPUT FILE************************
!-----------------------------------------------------------------------
    ! CALL INPUT_READIN
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
        x_vec = x_vec +dt*v_vec(:,1:2)
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
            ! CALL SPECULAR_REFLECTION(x_vec,x_vec_prev,v_vec,v_vec_prev,N_all,x_walls,num_walls,xmin,xmax,0)
            ! WRITE(*,*) "--- specular_reflection ---"
            CALL SPECULAR_REFLECTION('SIM_PARTICLES___',N_array,0)
        END IF

        ! Input from Source/Reservoir -------------------------------------------------
        IF (include_source .EQV. .true.) THEN
            ! WRITE(*,*) "--- initialize_source ---"

            CALL INITIALIZE_SOURCE
            N_added_total(ii) = N_entered            
        END IF

        IF (N_simulated > 0) THEN
            ! Update Cell Index -----------------------------------------------------------
            ! WRITE(*,*) "--- update_cell_index ---"
            CALL UPDATE_CELL_INDEX

            ! Collisions ------------------------------------------------------------------
            ! Just do separate collisions subroutine here
            ! WRITE(*,*) "--- run_collisions ---"
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
    ! t_final = time.time()-t0
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
    ! WRITE(*,*) "n_cells_vec = ",n_cells_vec
    WRITE(*,*) "(nx,ny) = (",((/ nx, ny /)),")"
    WRITE(*,*) "nt = ",nt
    WRITE(*,*) "computation time (total) = ",t_final
    WRITE(*,*) "computation time (BC's) = ",t_BC
    WRITE(*,*) "computation time (collisions) = ",t_collisions
    WRITE(*,*) "computation time (looping in collisions)=",t_loop
    WRITE(*,*) "computation time (test)=",t_test

    ! WRITE(*,*) "i_cell_vec = ",i_cell_vec
    ! WRITE(*,*) "Npc_slice = ",Npc_slice


    WRITE(*,*) '------ done ------'

    STOP
END PROGRAM MAIN


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


    ! ALLOCATE( x_cells_vec(1) )
    ! ALLOCATE( y_cells_vec(1) )
    nx = n_cells_vec(1)
    ! ny = n_cells_vec(2)
    ny = 1
    cy  = 1
    Npc_max = 250

    IF (use_homogenous_grid .EQV. .true.) THEN
        ! set up equally spaced grid points
        ALLOCATE( x_cells_vec(nx) )
        ALLOCATE( y_cells_vec(ny) )
        dx = ( x_lim(1,2) - x_lim(1,1) ) / nx
        DO i = 1 , nx
            x_cells_vec(i) = x_lim(1,1) + dx*(i-1)
        END DO
        dy = ( x_lim(2,2) - x_lim(2,1) ) / ny
        DO i = 1 , ny
            y_cells_vec(i) = x_lim(2,1) + dy*(i-1)
        END DO
    ELSE
        xmax = x_lim(1,2)
        alpha_x = -LOG(1./dx_factor)/xmax
        n_inf = 1/(dx_0*alpha_x)
        nmax = INT(n_inf*(1-EXP(-alpha_x*xmax)))
        ALLOCATE(i_range(nmax+1))
        ! CALL linspace(i_range,0.d0,REAL(nmax,8),nmax+1)
        i_range = (/ (i,i=0,nmax) /)

        ALLOCATE( x_cells_vec(nmax+1) )
        ALLOCATE( y_cells_vec(1) )
        x_cells_vec = -1/alpha_x*LOG(1-i_range/n_inf)
        ! WRITE(*,*) SHAPE(x_cells_vec)
        y_cells_vec = 0

        n_cells_vec(1) = nmax+1
        n_cells_vec(2) = 1
        nx = nmax+1
        ! ny = 1
    END IF

    ! nx = n_cells_vec(1)
    ! ny = n_cells_vec(2)
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
    N_all  = INT(V_total*n/Fn)
    N_simulated = N_all

    ! time step parameters ---------------------------------------------------------
    ! WRITE(*,*) INT(tmax/dt)+1
    nt = INT(tmax/dt)+1
    ALLOCATE(t_vec(nt))
    DO i = 1,nt
        t_vec(i) = (i-1)*dt
    END DO

    ! t0 = time.time() ##########
    CALL CPU_TIME(t0)
    t_final = 0
    ! ##############
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
    Num_s = INT(ns*Vs/Fn)
    ! WRITE(*,*) "Num_s =",Num_s 

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
        N_expected = N_all
        N_array = N_expected
    END IF

    ALLOCATE(x_vec(N_array,ndim))
    ALLOCATE(v_vec(N_array,3))
    ALLOCATE(i_cell_vec(N_array,2))
    ALLOCATE(i_cell_vec_prev(N_array,2))
    ALLOCATE(x_vec_prev(N_array,ndim))
    ALLOCATE(v_vec_prev(N_array,3))
    ALLOCATE(particles_in_cell(nx,Npc_max))

    ALLOCATE(reflected_out(N_array))
    ALLOCATE(reflected_in(N_array))
    ALLOCATE(removed_from_sim(N_array))
    removed_from_sim(:) = .false.

    ALLOCATE(in_column(N_array))
    ALLOCATE(in_cell(N_array))
    ALLOCATE(i_counting(N_array))
    i_counting = (/ (i,i=1,N_array) /)
    ALLOCATE(i_cur(N_array))

    ALLOCATE(xs_vec(Num_s,ndim))
    ALLOCATE(xs_vec_prev(Num_s,ndim))
    ALLOCATE(vs_vec(Num_s,3))
    ALLOCATE(vs_vec_prev(Num_s,3))
    ALLOCATE(counter_vec(Num_s))
    ALLOCATE(entered_Sim(Num_s))




    ALLOCATE(vr_max(nx,ny))
    vr_max(:,:) = vr_max_0

    IF (restart_simulation .EQV. .false.) THEN
        IF (N_all >= 0) THEN
            CALL RANDOM_NUMBER(x_vec)
            x_vec(:,1) = x_vec(:,1)*(xmax-xmin)+xmin
            x_vec(:,2) = x_vec(:,2)*(ymax-ymin)+ymin
            ! CALL RANDOM_NUMBER(v_vec)   ! need to convert this to normal distribution #########
            CALL RANDN(v_vec(:,1),N_all)
            CALL RANDN(v_vec(:,2),N_all)
            CALL RANDN(v_vec(:,3),N_all)
            ! CALL RANDN(v_vec(:,1),N_simulated)
            ! CALL RANDN(v_vec(:,2),N_simulated)
            ! CALL RANDN(v_vec(:,3),N_simulated)
            v_vec = v_vec*vth
            x_vec(:,2) = 0.5
        END IF

        ALLOCATE(N_total(nt))
        ALLOCATE(N_candidate_pairs_total(nt))
        ALLOCATE(N_accepted_pairs_total(nt))
        ALLOCATE(N_collisions_total(nt))
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
    IF (N_all > 0) THEN

        i_cell_vec_prev = i_cell_vec

        IF (use_homogenous_grid .EQV. .true.) THEN
            i_cell_vec(1:N_all,1) = CEILING( (x_vec(1:N_all,1)-xmin)/(xmax-xmin)*nx )
            ! i_cell_vec(:,1) = CEILING( (x_vec(:,1)-xmin)/(xmax-xmin)*nx )
            IF (ny>1) THEN
                i_cell_vec(1:N_all,2) = CEILING( (x_vec(1:N_all,2)-xmin)/(xmax-xmin)*ny )
                ! i_cell_vec(:,2) = CEILING( (x_vec(:,2)-xmin)/(xmax-xmin)*ny )
            ELSE
                i_cell_vec(1:N_all,2) = 1
                ! i_cell_vec(:,2) = 1
            END IF


        ELSE
            alpha_x  = -LOG(1/dx_factor)/xmax
            n_inf = 1/(dx_0*alpha_x)
            i_cell_vec(1:N_all,1) = FLOOR(n_inf*(1-EXP( -alpha_x*x_vec(1:N_all,1) )))+1
            ! i_cell_vec(:,1) = CEILING(n_inf*(1-EXP( -alpha_x*x_vec(:,1) )))


            ! WRITE(*,*) "alpha_x = ",alpha_x
            ! WRITE(*,*) "1/dx_factor = ",1/dx_factor
            ! WRITE(*,*) "xmax = ",xmax
            ! WRITE(*,*) "dx_0 = ",dx_0
            ! WRITE(*,*) "n_inf = ",n_inf

            ! WRITE(*,*) "val_1 = ",-alpha_x*x_vec(1:20,1)
            ! WRITE(*,*) "EXP(val_1) = ",EXP( -alpha_x*x_vec(1:20,1) )
            ! WRITE(*,*) "x_vec= ",x_vec(1:20,1)



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
                ! WHERE (x_vec(:,2) < ymid)
                !     i_cell_vec(:,2) = FLOOR( neg_offset - FLOOR( n_inf*(1-EXP(-alpha_y*(ymid - x_vec(:,2))))) )
                ! ELSEWHERE
                !     i_cell_vec(:,2) = CEILING( pos_offset - FLOOR( n_inf*(1-EXP(-alpha_y*(x_vec(:,2) - ymid)))) )
                ! END WHERE

                i_cell_vec(1:N_all,2) = i_cell_vec(1:N_all,2) + 1
                ! i_cell_vec(:,2) = i_cell_vec(:,2) + 1
            ELSE
                i_cell_vec(1:N_all,2) = 1
                ! i_cell_vec(:,2) = 1
            ENDIF

        END IF

        WHERE (removed_from_sim .EQV. .true.)
            i_cell_vec(:,1) = 0
            i_cell_vec(:,2) = 0
        ELSEWHERE
        END WHERE


        

    ENDIF
END SUBROUTINE UPDATE_CELL_INDEX





SUBROUTINE SPECULAR_REFLECTION(string_in,Num_r,counter)
! SUBROUTINE SPECULAR_REFLECTION(xr_vec,xr_vec_prev,vr_vec,vr_vec_prev,N_all,xr_walls,num_walls,xmin,xmax,counter)
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    ! REAL(8),DIMENSION(Num_r,ndim):: xr_vec,xr_vec_prev, x_coll, xr_vec_new
    ! REAL(8),DIMENSION(Num_r,3):: vr_vec,vr_vec_prev, vr_vec_new
    ! REAL(8),DIMENSION(num_walls,4)::xr_walls
    ! LOGICAL,DIMENSION(Num_r,num_walls):: collision_occured
    ! REAL(8),DIMENSION(Num_r,num_walls):: collision_dt
    ! REAL(8),DIMENSION(Num_r):: x0,y0,xt,yt,m,b,xc,yc, dt_cross
    ! LOGICAL,DIMENSION(Num_r):: reflected_in, reflected_out, first_collision
    ! INTEGER,DIMENSION(Num_r):: i_cross, i_first
    ! REAL(8):: xw1,xw2,yw1,yw2,xmin,xmax, temp
    ! INTEGER:: Num_r,num_walls,counter, i

    CHARACTER(16):: string_in
    REAL(8),DIMENSION(Num_r,ndim):: xr_vec,xr_vec_prev, x_coll, xr_vec_new
    REAL(8),DIMENSION(Num_r,3):: vr_vec,vr_vec_prev, vr_vec_new
    REAL(8),DIMENSION(num_walls,4)::xr_walls
    LOGICAL,DIMENSION(Num_r,num_walls):: collision_occured
    REAL(8),DIMENSION(Num_r,num_walls):: collision_dt
    REAL(8),DIMENSION(Num_r):: x0,y0,xt,yt,m,b,xc,yc, dt_cross
    LOGICAL,DIMENSION(Num_r):: first_collision,crossed
    INTEGER,DIMENSION(Num_r):: i_cross, i_first
    REAL(8):: xw1,xw2,yw1,yw2, temp
    INTEGER:: Num_r,counter, i

    CALL CPU_TIME(t0_BC)

    IF (string_in == 'SOURCE_PARTICLES') THEN
        xr_vec = xs_vec
        xr_vec_prev = xs_vec_prev
        vr_vec = vs_vec
        vr_vec_prev = vs_vec_prev
        ! xr_walls = x_walls(5:)
    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        xr_vec = x_vec
        xr_vec_prev = x_vec_prev
        vr_vec = v_vec
        vr_vec_prev = v_vec_prev
        xr_walls = x_walls
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF

    xr_vec_new = xr_vec
    vr_vec_new = vr_vec
    ! i_refl_out/in?
    collision_occured(:,:) = .false.
    reflected_in(:) = .false.
    reflected_out(:) = .false.

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
        ! CALL LOG2INT( i_cross, ((x0<xw1).and.(xw1<xt)) .or. ((x0>xw1).and.(xw1>xt)) , Num_r )
        crossed = ((x0<xw1).and.(xw1<xt)) .or. ((x0>xw1).and.(xw1>xt))
        dt_cross = (xc-x0)/vr_vec_prev(:,1)


        ! collision_occured(i_cross,i) = .true.
        ! collision_dt(i_cross,i) = dt_cross(i_cross)
        WHERE (crossed)
            collision_occured(:,i) = .true.
            collision_dt(:,i) = dt_cross
        ELSEWHERE
        END WHERE

        ! CALL LOG2INT( i_first, (collision_occured(:,i) .eqv. .true.) .and. (collision_dt(:,i) == MINVAL(collision_dt,2)) , Num_r)
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

        IF (string_in == 'SIM_PARTICLES___') THEN
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
        END IF

        ! other angles/wall-types go here



    END DO

    ! xr_vec = xr_vec_new
    ! vr_vec = vr_vec_new



    IF (string_in == 'SOURCE_PARTICLES') THEN
        xs_vec = xr_vec_new
        vs_vec = vr_vec_new
    ELSE IF (string_in == 'SIM_PARTICLES___') THEN
        x_vec = xr_vec_new
        v_vec = vr_vec_new
    ELSE
        WRITE(*,*) "Error identifying source/sim particles for specular reflection"
    END IF



    N_removed = 0
    ! remove(ignore) exiting particles from simulation
    IF (close_outlet .EQV. .false.) THEN
        WHERE( reflected_out .EQV. .true.)
            removed_from_sim = .true.
        ELSEWHERE
        END WHERE
        N_removed = N_removed + COUNT(reflected_out)
    END IF
    IF (close_inlet .EQV. .false.) THEN
        WHERE( reflected_in .EQV. .true.)
            removed_from_sim = .true.
        ELSEWHERE
        END WHERE
        N_removed = N_removed + COUNT(reflected_in)
    END IF
    N_simulated = N_simulated - N_removed
    ! WRITE(*,*) "N_simulated=",N_simulated
    ! WRITE(*,*) "N_removed=",N_removed

    CALL CPU_TIME(t_temp)
    t_BC = t_BC + (t_temp-t0_BC)
    
END SUBROUTINE SPECULAR_REFLECTION





SUBROUTINE INITIALIZE_SOURCE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    ! INTEGER,ALLOCATABLE, DIMENSION(:):: i_cur


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

    entered_sim = (xs_vec(:,1) > 0)
    N_entered = COUNT(entered_sim)
    ! ALLOCATE(i_cur(N_entered))
    ! i_cur(:) = 0
    i_cur(1:N_entered) = PACK(i_counting , entered_sim)

    ! WRITE(*,*) SIZE(i_cur),N_entered,N_all,N_simulated

    IF (N_entered > 0) THEN
        x_vec( (N_all+1):(N_all+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
        v_vec( (N_all+1):(N_all+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)
        N_all = N_all + N_entered
        N_simulated = N_simulated + N_entered
    END IF
    ! DEALLOCATE(i_cur)



END SUBROUTINE INITIALIZE_SOURCE



SUBROUTINE RUN_COLLISIONS
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    REAL(8):: rn, vr, cosT,sinT,alpha,alpha2
    INTEGER:: i,j,k
    REAL(8),DIMENSION(3)::normal_vec,v_temp !,v0,v1
    ! INTEGER,ALLOCATABLE, DIMENSION(:):: i_cur

    CALL CPU_TIME(t0_collisions)

    ! cy = 1

    ! DO cx = 1,n_cells_vec(1)
    DO cx = 1,nx
        CALL CPU_TIME(t0_test)

        ! CALL LOG2INT( i_column, i_cell_vec(:,1) == cx , N_all)
        in_column = (i_cell_vec(:,1) == cx)
        in_cell  = in_column !1D, cy=1


        ! collision processing
        N_candidate_pairs = 0
        N_accepted_pairs = 0
        Npc_cur = COUNT(in_column)
        Npc_slice(cx,cy) = Npc_cur

        ! ALLOCATE(i_cur(Npc_cur))
        i_cur(1:Npc_cur) = PACK(i_counting , in_column)
        CALL CPU_TIME(t_temp)
        t_test = t_test + (t_temp-t0_test)

        IF (Npc_cur >= 2) THEN
            N_candidate_pairs_real = .5*Npc_cur*(Npc_cur-1)*Fn*c_s*vr_max(cx,cy)*dt/Vc(cx,cy)   ! number of candidate collision pairs
            N_candidate_pairs_real = N_candidate_pairs_real + ncp_remainder(cx,cy)
            ncp_remainder(cx,cy) = MOD(N_candidate_pairs_real,1.0)
            N_candidate_pairs = FLOOR(N_candidate_pairs_real)
            ! IF (N_candidate_pairs > 1) THEN
            !     WRITE(*,*) "---"
            !     WRITE(*,*) "Npc_cur=",Npc_cur
            !     WRITE(*,*) "Fn=",Fn
            !     WRITE(*,*) "c_s=",c_s
            !     WRITE(*,*) "vr_max=",vr_max(cx,cy)
            !     WRITE(*,*) "dt=",dt
            !     WRITE(*,*) "Vc=",Vc(cx,cy)
            !     WRITE(*,*) "i_cur = ",i_cur
            !     WRITE(*,*) "N_candidate_pairs_real=",N_candidate_pairs_real
            !     WRITE(*,*) "N_candidate_pairs=",N_candidate_pairs
            !     WRITE(*,*) "ncp_remainder=",ncp_remainder(cx,cy)
            ! END IF

            CALL CPU_TIME(t0_loop)
            DO k=1,N_candidate_pairs

                CALL RANDOM_NUMBER(rn)
                i = FLOOR(rn*Npc_cur)+1
                j = i
                DO WHILE (i==j)
                    CALL RANDOM_NUMBER(rn)
                    j = FLOOR(rn*Npc_cur)+1
                END DO
                ! i_pair = i_cur((/i,j/)) 
                ! i_pair = i_counting( i_cur((/i,j/)) )

                ! WRITE(*,*) "--- got here, k= ", k, ", Npc_cur=",Npc_cur
                ! WRITE(*,*) "(i,j) = ",i,j
                ! WRITE(*,*) "(i_cur(i),i_cur(j) = ",i_cur(i),i_cur(j)
                ! WRITE(*,*) "size(i_counting) = ",SIZE(i_counting)
                ! WRITE(*,*) "N_all,N_array = ",N_all,N_array

                i = i_counting(i_cur(i))
                j = i_counting(i_cur(j))

                ! WRITE(*,*) "--- got here ---"


                vr = SQRT( (v_vec(i,1)-v_vec(j,1))**2 + (v_vec(i,2)-v_vec(j,2))**2 + (v_vec(i,3)-v_vec(j,3))**2 )

                IF (vr > vr_max(cx,cy)) THEN
                    vr_max(cx,cy) = vr
                ENDIF
                CALL RANDOM_NUMBER(rn)
                IF (vr/vr_max(cx,cy) > rn) THEN
                    ! process collision

                    CALL RANDOM_NUMBER(alpha)
                    cosT = 1-2*alpha
                    sinT = SQRT(1-cosT**2)
                    CALL RANDOM_NUMBER(alpha2)
                    normal_vec(1) = cosT
                    normal_vec(2) = sinT*COS(2*Pi*alpha2)
                    normal_vec(3) = sinT*SIN(2*Pi*alpha2)

                    v_temp = SUM( (v_vec(j,:)-v_vec(i,:))*normal_vec ) * normal_vec
                    v_vec(i,:) = v_vec(i,:) + v_temp
                    v_vec(j,:) = v_vec(j,:)- v_temp

                    N_collisions_total(ii) = N_collisions_total(ii) + 1
                    N_accepted_pairs = N_accepted_pairs + 1
                ENDIF

            END DO
            CALL CPU_TIME(t_temp)
            t_loop = t_loop + (t_temp-t0_loop)
        END IF

        N_candidate_pairs_total(ii) = N_candidate_pairs_total(ii) + N_candidate_pairs
        N_accepted_pairs_total(ii) = N_accepted_pairs_total(ii) + N_accepted_pairs


        ! DEALLOCATE(i_cur)

    END DO

    CALL CPU_TIME(t_temp)
    t_collisions = t_collisions + (t_temp-t0_collisions)

END SUBROUTINE RUN_COLLISIONS



SUBROUTINE SAVE_DATA
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    CHARACTER(80)::filename

    ! save position data
    WRITE(filename,"('Output/data/x_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) x_vec
    ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! WRITE(1,"(E12.5)") x_vec
    CLOSE(1)

    ! save velocity data
    WRITE(filename,"('Output/data/v_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) v_vec
    CLOSE(1)

    ! save cell data (this is probably not required)
    WRITE(filename,"('Output/data/i_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) i_cell_vec
    CLOSE(1)

    ! save num_per_cell data? This also seems way unnecessary
    WRITE(filename,"('Output/data/Npc_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) Npc_slice
    CLOSE(1)


    ! save miscellaneous data
    WRITE(filename,"('Output/data/data.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    WRITE(1,"(A)") "*n_collisions"   
    WRITE(1,"(I0)") n_collisions
    WRITE(1,"(A)") "*n_collisions_total"   
    WRITE(1,"(I0)") n_collisions_total
    ! WRITE(1,"(A)") "N_candidate_pairs"   
    ! WRITE(1,"(I0)") N_candidate_pairs
    ! WRITE(1,"(A)") "N_accepted_pairs"   
    ! WRITE(1,"(I0)") N_accepted_pairs
    WRITE(1,"(A)") "*N_candidate_pairs_total"   
    WRITE(1,"(I0)") N_candidate_pairs_total
    WRITE(1,"(A)") "*N_accepted_pairs_total"   
    WRITE(1,"(I0)") N_accepted_pairs_total
    WRITE(1,"(A)") "*N_added_total"   
    WRITE(1,"(I0)") N_added_total
    WRITE(1,"(A)") "*ncp_remainder"   
    WRITE(1,"(E12.5)") ncp_remainder
    WRITE(1,"(A)") "*N_added"   
    WRITE(1,"(I0)") N_added
    WRITE(1,"(A)") "*N_total"   
    WRITE(1,"(I0)") N_total
    WRITE(1,"(A)") "*t_final"   
    WRITE(1,"(E12.5)") t_final
    WRITE(1,"(A)") "*t_BC"   
    WRITE(1,"(E12.5)") t_BC
    WRITE(1,"(A)") "*t_collisions"   
    WRITE(1,"(E12.5)") t_collisions
    WRITE(1,"(A)") "*t_loop"   
    WRITE(1,"(E12.5)") t_loop
    WRITE(1,"(A)") "*n"   
    WRITE(1,"(E12.5)") n
    WRITE(1,"(A)") "*ns"   
    WRITE(1,"(E12.5)") ns
    WRITE(1,"(A)") "*Fn"   
    WRITE(1,"(E12.5)") Fn
    WRITE(1,"(A)") "*nx"   
    WRITE(1,"(I0)") nx
    WRITE(1,"(A)") "*ny"   
    WRITE(1,"(I0)") ny
    WRITE(1,"(A)") "*tmax"   
    WRITE(1,"(E12.5)") tmax
    WRITE(1,"(A)") "*nt"   
    WRITE(1,"(I0)") nt
    WRITE(1,"(A)") "*dt"   
    WRITE(1,"(E12.5)") dt
    WRITE(1,"(A)") "*dt_to_save"   
    WRITE(1,"(I0)") dt_to_save
    WRITE(1,"(A)") "*n_saved"   
    WRITE(1,"(I0)") n_saved
    WRITE(1,"(A)") "*include_source"   
    WRITE(1,"(L)") include_source
    WRITE(1,"(A)") "*close_inlet"   
    WRITE(1,"(L)") close_inlet
    WRITE(1,"(A)") "*include_gun_boundaries"   
    WRITE(1,"(L)") include_gun_boundaries
    WRITE(1,"(A)") "*use_homogenous_grid"   
    WRITE(1,"(L)") use_homogenous_grid
    WRITE(1,"(A)") "*dx_0"   
    WRITE(1,"(E12.5)") dx_0
    WRITE(1,"(A)") "*dx_factor"   
    WRITE(1,"(E12.5)") dx_factor
    WRITE(1,"(A)") "*dy_factor"   
    WRITE(1,"(E12.5)") dy_factor
    WRITE(1,"(A)") "*x_cells_vec"   
    WRITE(1,"(E12.5)") x_cells_vec
    WRITE(1,"(A)") "*y_cells_vec"   
    WRITE(1,"(E12.5)") y_cells_vec
    WRITE(1,"(A)") "*x_lim"   
    WRITE(1,"(E12.5)") x_lim
    WRITE(1,"(A)") "*x_walls"   
    WRITE(1,"(E12.5)") x_walls
    WRITE(1,"(A)") "*it_last"   
    WRITE(1,"(I0)") ii
    CLOSE(1)



    ! WRITE(var_cur,"("N_candidate_pairs_total" )")
    ! WRITE(filename,"('Output/data/',,'.txt')") var_cur
    ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! WRITE(1,"(A)") var_cur
    ! WRITE(1,"(I0)") N_candidate_pairs_total
    ! CLOSE(1)
    ! WRITE(filename,"('Output/data/data.txt')")
    ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! WRITE(1,"(A)") "N_accepted_pairs_total"   
    ! WRITE(1,"(I0)") N_accepted_pairs_total
    ! CLOSE(1)


    n_saved = n_saved + 1

    CALL CPU_TIME(t_temp)
    t_collisions = t_collisions + (t_temp-t0_collisions)

END SUBROUTINE SAVE_DATA


