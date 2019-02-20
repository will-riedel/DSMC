SUBROUTINE INITIALIZE_SOURCE
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    REAL(8)::rn


    CALL CPU_TIME(t0_source)


    IF (include_source .EQV. .true.) THEN

        IF (source_type(1:10) == "DOWNSTREAM") THEN
            CALL INITIALIZE_SOURCE_ONE_STREAM_DOWNSTREAM
        ELSE

            CALL RANDOM_NUMBER(rn)
            ! IF (rn < Num_s_exact) THEN
            IF (rn < Num_s_frac) THEN
                Num_s = CEILING(Num_s_exact)
                Num_s_b = CEILING(Num_s_exact_b)
            ELSE
                Num_s = FLOOR(Num_s_exact)
                Num_s_b = FLOOR(Num_s_exact_b)
            END IF

            ! WRITE(*,*) "GH 5.2: Num_s_exact, Num_s, source_type= ",Num_s_exact,Num_s,source_type

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

    CALL CPU_TIME(t_temp)
    t_source = t_source + (t_temp-t0_source)

END SUBROUTINE INITIALIZE_SOURCE




SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    INTEGER::j,i_temp

    ! WRITE(*,*) "GH 5.3: shape(xs_vec), shape(vs_vec) = ",SHAPE(xs_vec),SHAPE(vs_vec)
    ! WRITE(*,*) "xs_min,xs_max,ys_min,ys_max=",xs_min,xs_max,ys_min,ys_max
    ! ! WRITE(*,*) "vth,dt=",vth,dt


    CALL RANDOM_NUMBER(xs_vec(1:Num_s,:))
    ! CALL RANDOM_NUMBER(xs_vec(:,1))
    xs_vec(1:Num_s,1) = xs_vec(1:Num_s,1)*(xs_max-xs_min) + xs_min
    xs_vec(1:Num_s,2) = xs_vec(1:Num_s,2)*(ys_max-ys_min) + ys_min
    !xs_vec(1:Num_s,3) = 0
    ! xs_vec(:,2) = 0.5
    ! WRITE(*,*) "GH 5.4"


    CALL RANDN(vs_vec(1:Num_s,1),Num_s) 
    CALL RANDN(vs_vec(1:Num_s,2),Num_s) 
    CALL RANDN(vs_vec(1:Num_s,3),Num_s) 
    vs_vec = vs_vec*vth
    ! WRITE(*,*) "GH 5.5"


    xs_vec_prev = xs_vec
    vs_vec_prev = vs_vec
    xs_vec = xs_vec + dt*vs_vec(:,1:2)
    !xs_vec = xs_vec + dt*vs_vec(:,:)


    IF (geometry_type(1:11) == "CYLINDRICAL") THEN   
        CALL ROTATE_SOURCE_TO_AXIAL_PLANE
    END IF

    CALL SPECULAR_REFLECTION_SOURCE


    entered_sim = (xs_vec(1:Num_s,1) > xs_max)
    ! entered_sim = (vs_vec(1:Num_s,1) > 0)
    N_entered = COUNT(entered_sim)


    ! WRITE(*,*) "GH 5.6: Num_s, N_entered = ",Num_s,N_entered
    IF (N_entered > 0) THEN
        ! i_cur(1:N_entered) = PACK(i_counting , entered_sim)

        ! i_cur(1:N_entered) = PACK(i_counting(1:Num_s) , entered_sim(1:Num_s))
        ! x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
        ! v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)

        i_temp = N_simulated+1
        DO j = 1,Num_s
            IF (entered_sim(j) .eqv. .true.) THEN
                x_vec(i_temp,:) = xs_vec(j,:)
                v_vec(i_temp,:) = vs_vec(j,:)

                weight_factor_vec(i_temp) = &
                            1 + RWF*x_vec(i_temp,2)/ymax


                i_temp = i_temp+1
            END IF
        END DO




        
        N_simulated = N_simulated + N_entered
    END IF

    N_added_total(ii) = N_entered            
    ! N_added_total(ii) = N_added_total(ii) + N_entered            


END SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM



SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM_DOWNSTREAM
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    REAL(8)::rn, Num_s_exact_cur, Num_s_frac_cur, scale_factor
    REAL(8)::alpha1,alpha2,scale_max,t_offset
    ! INTEGER::j

    scale_factor = 1

    ! exponential scaling based on no-collision case simulated on 10/20
    scale_max = 2.45/6.8    ! based on simulation at ns = 6.8d25
    alpha1 = 16900.
    alpha2 = 4550.
    t_offset = 15.d-6
    IF (t < ts+t_offset) THEN
        scale_factor = scale_max*( 1-EXP(-alpha1*t) )
    ELSE
        scale_factor = scale_max*EXP( -alpha2*(t-(ts+t_offset)) )
    END IF

    ! WRITE(*,*) scale_max,alpha1,alpha2,t_offset


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
        CALL RANDOM_NUMBER(xs_vec_prev(1:Num_s,:))

        xs_vec(1:Num_s,1) = xs_vec_prev(1:Num_s,1)*(xs_max-xs_min)  + xs_min
        xs_vec(1:Num_s,2) = (xs_vec_prev(1:Num_s,1)*(xs_max-xs_min) + xs_min) + b_source_A

        xs_vec(1:Num_s,1) = xs_vec(1:Num_s,1) - xs_vec_prev(1:Num_s,2)*(hs)/SQRT(2.)
        xs_vec(1:Num_s,2) = xs_vec(1:Num_s,2) + xs_vec_prev(1:Num_s,2)*(hs)/SQRT(2.)


        CALL RANDN(vs_vec(1:Num_s,1),Num_s) 
        CALL RANDN(vs_vec(1:Num_s,2),Num_s) 
        CALL RANDN(vs_vec(1:Num_s,3),Num_s) 
        vs_vec = vs_vec*vth

        xs_vec_prev = xs_vec
        vs_vec_prev = vs_vec
        xs_vec = xs_vec + dt*vs_vec(:,1:2)

        entered_sim = ( (xs_vec(1:Num_s,2) > ((-1)*xs_vec(1:Num_s,1) + b_source_barrier) ) .and. & 
                        (xs_vec(1:Num_s,2) > ( (1)*xs_vec(1:Num_s,1) + b_source_A)) .and. & 
                        (xs_vec(1:Num_s,2) < ( (1)*xs_vec(1:Num_s,1) + b_source_B)) )
        N_entered = COUNT(entered_sim)



        IF (N_entered > 0) THEN
            i_cur(1:N_entered) = PACK(i_counting(1:Num_s) , entered_sim(1:Num_s))

            x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
            v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)

            weight_factor_vec( (N_simulated+1):(N_simulated+1+N_entered) ) = &
                            1 + RWF*x_vec( (N_simulated+1):(N_simulated+1+N_entered) , 2 )/ymax


            N_simulated = N_simulated + N_entered
            ! WRITE(*,*) "got here s-7"
        END IF

        N_added_total(ii) = N_entered            

    END IF

END SUBROUTINE INITIALIZE_SOURCE_ONE_STREAM_DOWNSTREAM







SUBROUTINE INITIALIZE_SOURCE_TWO_STREAM
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE


    ! Initialize beam from near side --------------------------------------------------------------
    ! CALL RANDOM_NUMBER(xs_vec)
    CALL RANDOM_NUMBER(xs_vec(1:Num_s,:))
    xs_vec(1:Num_s,1) = xs_vec(1:Num_s,1)*(xs_max-xs_min) + xs_min
    xs_vec(1:Num_s,2) = xs_vec(1:Num_s,2)*(ys_max-ys_min) + ys_min
    ! xs_vec(:,2) = 0.5

    ! vs_vec(:,1) = v_beam
    ! vs_vec(:,2) = 0
    ! vs_vec(:,3) = 0
    ! WRITE(*,*) "Num_s=",Num_s
    ! WRITE(*,*) "shape = ",SHAPE(vs_vec)
    CALL RANDN(vs_vec(1:Num_s,1),Num_s) 
    CALL RANDN(vs_vec(1:Num_s,2),Num_s) 
    CALL RANDN(vs_vec(1:Num_s,3),Num_s) 
    vs_vec = vs_vec*vth !/10.
    vs_vec(1:Num_s,1) = vs_vec(1:Num_s,1) + v_beam

    ! WRITE(*,*) "got here 5"
    ! WRITE(*,*) "Num_s=",Num_s
    ! WRITE(*,*) "shape = ",SHAPE(vs_vec)

    ! xs_vec_prev = xs_vec
    ! vs_vec_prev = vs_vec
    xs_vec(1:Num_s,:) = xs_vec(1:Num_s,:) + dt*vs_vec(1:Num_s,1:2)


    ! (not needed right now)
    Num_r = Num_s
    CALL SPECULAR_REFLECTION_SOURCE

    ! entered_sim = (xs_vec(:,1) > 0)
    entered_sim = (xs_vec(1:Num_s,1) > xs_max)
    N_entered = COUNT(entered_sim)


    IF (N_entered > 0) THEN
        ! i_cur(1:N_entered) = PACK(i_counting , entered_sim)
        i_cur(1:N_entered) = PACK(i_counting(1:Num_s) , entered_sim(1:Num_s))
        x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
        v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)

        weight_factor_vec( (N_simulated+1):(N_simulated+1+N_entered) ) = &
                            1 + RWF*x_vec( (N_simulated+1):(N_simulated+1+N_entered) , 2 )/ymax

        ! N_all = N_all + N_entered
        N_simulated = N_simulated + N_entered
    END IF
    N_added_total(ii) = N_added_total(ii) + N_entered            



    ! Initialize beam from far side --------------------------------------------------------------
    CALL RANDOM_NUMBER(xs_vec)
    ! CALL RANDOM_NUMBER(xs_vec(:,1))
    xs_vec(1:Num_s_b,1) = xs_vec(1:Num_s_b,1)*(xs2_max-xs2_min) + xs2_min
    xs_vec(1:Num_s_b,2) = xs_vec(1:Num_s_b,2)*(ys_max-ys_min) + ys_min
    ! xs_vec(:,2) = 0.5

    ! vs_vec(:,1) = -v_beam
    ! vs_vec(:,2) = 0
    ! vs_vec(:,3) = 0
    CALL RANDN(vs_vec(1:Num_s_b,1),Num_s_b) 
    CALL RANDN(vs_vec(1:Num_s_b,2),Num_s_b) 
    CALL RANDN(vs_vec(1:Num_s_b,3),Num_s_b) 
    vs_vec = vs_vec*vth !/10.
    vs_vec(1:Num_s_b,1) = vs_vec(1:Num_s_b,1) - v_beam


    ! xs_vec_prev = xs_vec
    ! vs_vec_prev = vs_vec
    xs_vec(1:Num_s_b,:) = xs_vec(1:Num_s_b,:) + dt*vs_vec(1:Num_s_b,1:2)


    ! (not needed right now)
    Num_r = Num_s_b
    CALL SPECULAR_REFLECTION_SOURCE

    ! entered_sim = (xs_vec(:,1) > 0)
    entered_sim = (xs_vec(1:Num_s_b,1) < xs2_min)
    N_entered = COUNT(entered_sim)


    IF (N_entered > 0) THEN
        i_cur(1:N_entered) = PACK(i_counting , entered_sim)
        i_cur(1:N_entered) = PACK(i_counting(1:Num_s_b) , entered_sim(1:Num_s_b))
        x_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = xs_vec(i_cur(1:N_entered),:)
        v_vec( (N_simulated+1):(N_simulated+1+N_entered) , : ) = vs_vec(i_cur(1:N_entered),:)

        weight_factor_vec( (N_simulated+1):(N_simulated+1+N_entered) ) = &
                            1 + RWF*x_vec( (N_simulated+1):(N_simulated+1+N_entered) , 2 )/ymax

        ! N_all = N_all + N_entered
        N_simulated = N_simulated + N_entered
    END IF
    N_added_total(ii) = N_added_total(ii) + N_entered            



END SUBROUTINE INITIALIZE_SOURCE_TWO_STREAM









