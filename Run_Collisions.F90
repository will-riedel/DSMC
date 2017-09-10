SUBROUTINE RUN_COLLISIONS
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE

    REAL(8):: rn, vr, cosT,sinT,alpha,alpha2
    INTEGER:: i,j,k,u, cell_start,cell_end
    REAL(8),DIMENSION(3)::normal_vec,v_temp !,v0,v1

    CALL CPU_TIME(t0_collisions)

    DO cx = 1,nx

        ! collision processing
        N_candidate_pairs = 0
        N_accepted_pairs = 0
        Npc_cur = Npc_slice(cx,cy)

        IF (Npc_cur >= 2) THEN
            cell_start = starting_index(cx,cy)
            cell_end = cell_start + Npc_cur - 1
            i_cur(1:Npc_cur) = (/ (u,u=cell_start,cell_end) /)


            N_candidate_pairs_real = .5*Npc_cur*(Npc_cur-1)*Fn*c_s*vr_max(cx,cy)*dt/Vc(cx,cy)   ! number of candidate collision pairs
            N_candidate_pairs_real = N_candidate_pairs_real + ncp_remainder(cx,cy)
            ncp_remainder(cx,cy) = MOD(N_candidate_pairs_real,1.0)
            N_candidate_pairs = FLOOR(N_candidate_pairs_real)

            CALL CPU_TIME(t0_loop)
            DO k=1,N_candidate_pairs

                CALL RANDOM_NUMBER(rn)
                i = FLOOR(rn*Npc_cur)+1
                j = i
                DO WHILE (i==j)
                    CALL RANDOM_NUMBER(rn)
                    j = FLOOR(rn*Npc_cur)+1
                END DO

                ! set i,j to the total indices of the colliding particles (originally were the indices among particles in the cell)
                i = i_counting(i_cur(i))
                j = i_counting(i_cur(j))

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
                    v_vec(j,:) = v_vec(j,:) - v_temp

                    N_collisions_total(ii) = N_collisions_total(ii) + 1
                    N_accepted_pairs = N_accepted_pairs + 1
                ENDIF

            END DO
            CALL CPU_TIME(t_temp)
            t_loop = t_loop + (t_temp-t0_loop)
        END IF

        N_candidate_pairs_total(ii) = N_candidate_pairs_total(ii) + N_candidate_pairs
        N_accepted_pairs_total(ii) = N_accepted_pairs_total(ii) + N_accepted_pairs

    END DO

    CALL CPU_TIME(t_temp)
    t_collisions = t_collisions + (t_temp-t0_collisions)

END SUBROUTINE RUN_COLLISIONS







