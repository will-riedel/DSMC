SUBROUTINE INITIALIZE_SOURCE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    REAL(8)::rn



    IF (include_source .EQV. .true.) THEN

        IF (source_type(1:10) == "DOWNSTREAM") THEN
            CALL INITIALIZE_SOURCE_ONE_STREAM_DOWNSTREAM
        ELSE

            CALL RANDOM_NUMBER(rn)
            ! IF (rn < Num_s_exact) THEN
            IF (rn < Num_s_frac) THEN
                Num_s = CEILING(Num_s_exact)
            ELSE
                Num_s = FLOOR(Num_s_exact)
            END IF

            IF (Num_s > 0) THEN
                IF (source_type(1:5) == "2BEAM") THEN
                    CALL INITIALIZE_SOURCE_TWO_STREAM
                ELSE IF (source_type(1:6) == "NORMAL") THEN
                    CALL INITIALIZE_SOURCE_ONE_STREAM
                ! ELSE IF (source_type(1:10) == "DOWNSTREAM") THEN
                !     CALL INITIALIZE_SOURCE_ONE_STREAM_DOWNSTREAM
                END IF
            END IF

        END IF

    END IF

END SUBROUTINE INITIALIZE_SOURCE




SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    
    CALL RANDOM_NUMBER(xs_vec)
    ! CALL RANDOM_NUMBER(xs_vec(:,1))
    xs_vec(:,1) = xs_vec(:,1)*(xs_max-xs_min) + xs_min
    xs_vec(:,2) = xs_vec(:,2)*(ys_max-ys_min) + ys_min
    ! xs_vec(:,2) = 0.5

    CALL RANDN(vs_vec(:,1),Num_s) 
    CALL RANDN(vs_vec(:,2),Num_s) 
    CALL RANDN(vs_vec(:,3),Num_s) 
    vs_vec = vs_vec*vth

    xs_vec_prev = xs_vec
    vs_vec_prev = vs_vec
    xs_vec = xs_vec + dt*vs_vec(:,1:2)


    CALL SPECULAR_REFLECTION_SOURCE

    entered_sim = (xs_vec(:,1) > xs_max)
    N_entered = COUNT(entered_sim)
    i_cur(1:N_entered) = PACK(i_counting , entered_sim)


    IF (N_entered > 0) THEN
        x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
        v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)
        N_simulated = N_simulated + N_entered
    END IF

    N_added_total(ii) = N_entered            


END SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM



SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM_DOWNSTREAM
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    REAL(8)::rn, Num_s_cur, Num_s_exact_cur, Num_s_frac_cur, scale_factor

    
    ! scale_factor = t/tmax
    scale_factor = 1

    Num_s_exact_cur = Num_s_exact*scale_factor
    Num_s_frac_cur = Num_s_exact_cur - FLOOR(Num_s_exact_cur)

    ! scale up the source density based on the assumed input ramp-up profile
    CALL RANDOM_NUMBER(rn)
    ! IF (rn < Num_s_exact) THEN
    IF (rn < Num_s_frac_cur) THEN
        Num_s = CEILING(Num_s_exact_cur)
    ELSE
        Num_s = FLOOR(Num_s_exact_cur)
    END IF



    IF (Num_s > 0) THEN

        CALL RANDOM_NUMBER(xs_vec_prev)

        xs_vec(:,1) = xs_vec_prev(:,1)*(xs_max-xs_min)  + xs_min
        xs_vec(:,2) = (xs_vec_prev(:,1)*(xs_max-xs_min) + xs_min) + b_source_A

        xs_vec(:,1) = xs_vec(:,1) - xs_vec_prev(:,2)*(hs)/SQRT(2.)
        xs_vec(:,2) = xs_vec(:,2) + xs_vec_prev(:,2)*(hs)/SQRT(2.)


        CALL RANDN(vs_vec(:,1),Num_s) 
        CALL RANDN(vs_vec(:,2),Num_s) 
        CALL RANDN(vs_vec(:,3),Num_s) 
        vs_vec = vs_vec*vth

        xs_vec_prev = xs_vec
        vs_vec_prev = vs_vec
        xs_vec = xs_vec + dt*vs_vec(:,1:2)



        entered_sim = ( (xs_vec(:,2) > ((-1)*xs_vec(:,1) + b_source_barrier) ) .and. & 
                        (xs_vec(:,2) > ( (1)*xs_vec(:,1) + b_source_A)) .and. & 
                        (xs_vec(:,2) < ( (1)*xs_vec(:,1) + b_source_B)) )
        N_entered = COUNT(entered_sim)
        i_cur(1:N_entered) = PACK(i_counting , entered_sim)


        IF (N_entered > 0) THEN
            x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
            v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)
            N_simulated = N_simulated + N_entered
        END IF

        N_added_total(ii) = N_entered            

    END IF

END SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM_DOWNSTREAM







SUBROUTINE INITIALIZE_SOURCE_TWO_STREAM
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE


    ! Initialize beam from near side --------------------------------------------------------------
    CALL RANDOM_NUMBER(xs_vec)
    ! CALL RANDOM_NUMBER(xs_vec(:,1))
    xs_vec(:,1) = xs_vec(:,1)*(xs_max-xs_min) + xs_min
    xs_vec(:,2) = xs_vec(:,2)*(ys_max-ys_min) + ys_min
    ! xs_vec(:,2) = 0.5

    ! vs_vec(:,1) = v_beam
    ! vs_vec(:,2) = 0
    ! vs_vec(:,3) = 0
    ! WRITE(*,*) "Num_s=",Num_s
    ! WRITE(*,*) "shape = ",SHAPE(vs_vec)
    CALL RANDN(vs_vec(:,1),Num_s) 
    CALL RANDN(vs_vec(:,2),Num_s) 
    CALL RANDN(vs_vec(:,3),Num_s) 
    vs_vec = vs_vec*vth/10.
    vs_vec(:,1) = vs_vec(:,1) + v_beam

    ! WRITE(*,*) "got here 5"
    ! WRITE(*,*) "Num_s=",Num_s
    ! WRITE(*,*) "shape = ",SHAPE(vs_vec)

    xs_vec_prev = xs_vec
    vs_vec_prev = vs_vec
    xs_vec = xs_vec + dt*vs_vec(:,1:2)


    ! (not needed right now)
    ! CALL SPECULAR_REFLECTION_SOURCE

    ! entered_sim = (xs_vec(:,1) > 0)
    entered_sim = (xs_vec(:,1) > xs_max)
    N_entered = COUNT(entered_sim)
    i_cur(1:N_entered) = PACK(i_counting , entered_sim)


    IF (N_entered > 0) THEN
        x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
        v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)
        ! N_all = N_all + N_entered
        N_simulated = N_simulated + N_entered
    END IF
    N_added_total(ii) = N_added_total(ii) + N_entered            



    ! Initialize beam from far side --------------------------------------------------------------
    CALL RANDOM_NUMBER(xs_vec)
    ! CALL RANDOM_NUMBER(xs_vec(:,1))
    xs_vec(:,1) = xs_vec(:,1)*(xs2_max-xs2_min) + xs2_min
    xs_vec(:,2) = xs_vec(:,2)*(ys_max-ys_min) + ys_min
    ! xs_vec(:,2) = 0.5

    ! vs_vec(:,1) = -v_beam
    ! vs_vec(:,2) = 0
    ! vs_vec(:,3) = 0
    CALL RANDN(vs_vec(:,1),Num_s) 
    CALL RANDN(vs_vec(:,2),Num_s) 
    CALL RANDN(vs_vec(:,3),Num_s) 
    vs_vec = vs_vec*vth/10.
    vs_vec(:,1) = vs_vec(:,1) - v_beam


    xs_vec_prev = xs_vec
    vs_vec_prev = vs_vec
    xs_vec = xs_vec + dt*vs_vec(:,1:2)


    ! (not needed right now)
    ! CALL SPECULAR_REFLECTION_SOURCE

    ! entered_sim = (xs_vec(:,1) > 0)
    entered_sim = (xs_vec(:,1) < xs2_min)
    N_entered = COUNT(entered_sim)
    i_cur(1:N_entered) = PACK(i_counting , entered_sim)


    IF (N_entered > 0) THEN
        x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
        v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)
        ! N_all = N_all + N_entered
        N_simulated = N_simulated + N_entered
    END IF
    N_added_total(ii) = N_added_total(ii) + N_entered            



END SUBROUTINE INITIALIZE_SOURCE_TWO_STREAM


