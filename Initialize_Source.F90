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