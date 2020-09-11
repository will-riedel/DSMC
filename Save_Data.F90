SUBROUTINE SAVE_DATA
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    ! CHARACTER(80)::filename

    CALL CPU_TIME(t0_save)

    ! use format = ES11.2E3 for output?

    !WRITE(*,*) "GH 8.1"

    ! save position data
    WRITE(filename,"('/x_',I7.7,'.txt')") (ii-1)
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1,"(E12.5)") x_vec(1:N_simulated,:)
    WRITE(1,"(E12.5E3)") x_vec(1:N_simulated,:)
    CLOSE(1)

    !WRITE(*,*) "GH 8.2"

    ! save velocity data
    WRITE(filename,"('/v_',I7.7,'.txt')") (ii-1)
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1,"(E12.5)") v_vec(1:N_simulated,:)
    WRITE(1,"(E12.5E3)") v_vec(1:N_simulated,:)
    CLOSE(1)

    ! ! save weight factor data
    ! WRITE(filename,"('/wf_',I7.7,'.txt')") (ii-1)
    ! filename = dir_cur(1:dir_cur_length) // filename
    ! ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! ! WRITE(1,"(E12.5)") v_vec(1:N_simulated,:)
    ! WRITE(1,"(E12.5E3)") weight_factor_vec(1:N_simulated)
    ! CLOSE(1)

    ! WRITE(*,*) "GH 8.3"


    ! save cell data (this is probably not required)
    WRITE(filename,"('/i_',I7.7,'.txt')") (ii-1)
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) i_cell_vec(1:N_simulated,:)
    WRITE(1,"(I0)") i_cell_vec(1:N_simulated,:)
    CLOSE(1)

    ! WRITE(*,*) "GH 8.4"


    ! ! save num_per_cell data? This also seems way unnecessary
    WRITE(filename,"('/Npc_',I7.7,'.txt')") (ii-1)
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) nx*ny
    ! WRITE(1) Npc_slice
    WRITE(1,"(I0)") Npc_slice
    CLOSE(1)

    ! WRITE(*,*) "GH 8.5"


    ! save some tracking vectors
    WRITE(filename,"('/n_collisions_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    ! CALL SYSTEM( "mv " // filename // ' ' // dir_cur(1:dir_cur_length) // "_temp.txt" )
    ! WRITE(*,*) ("mv " // filename // ' ' // dir_cur(1:dir_cur_length) // "_temp.txt")
    ! IF (ii > 50) THEN
    !     STOP
    ! END IF
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='direct',recl=nt)
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! WRITE(1) nt
    ! WRITE(1) n_collisions_total
    WRITE(1,"(I0)") n_collisions_total
    CLOSE(1)
    ! CALL SYSTEM( "rm " // dir_cur(1:dir_cur_length) // "_temp.txt" )
    ! WRITE(*,*) ( "rm " // dir_cur(1:dir_cur_length) // "_temp.txt" )

    ! WRITE(*,*) "GH 8.6"


    WRITE(filename,"('/N_candidate_pairs_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) nt
    ! WRITE(1) N_candidate_pairs_total
    WRITE(1,"(I0)") N_candidate_pairs_total
    CLOSE(1)

    ! WRITE(*,*) "GH 8.7"

    WRITE(filename,"('/N_accepted_pairs_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) nt
    ! WRITE(1) N_accepted_pairs_total
    WRITE(1,"(I0)") N_accepted_pairs_total
    CLOSE(1)

    ! WRITE(*,*) "GH 8.8"

    WRITE(filename,"('/N_added_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) nt
    ! WRITE(1) N_added_total
    WRITE(1,"(I0)") N_added_total
    CLOSE(1)
    ! WRITE(filename,"('/ncp_remainder.txt')")
    ! filename = dir_cur(1:dir_cur_length) // filename
    ! ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! ! WRITE(1) nx*ny
    ! ! WRITE(1) ncp_remainder
    ! WRITE(1,"(E12.5)") ncp_remainder
    ! CLOSE(1)

    ! WRITE(*,*) "GH 8.9"

    WRITE(filename,"('/N_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) nt
    ! WRITE(1) N_total
    WRITE(1,"(I0)") N_total
    CLOSE(1)

    WRITE(filename,"('/NWF_escaped_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) nt
    ! WRITE(1) N_total
    WRITE(1,"(E12.5E3)") NWF_escaped_total
    CLOSE(1)

    ! WRITE(*,*) "GH 8.10"

    WRITE(filename,"('/x_walls.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED",access='stream')
    ! WRITE(1) num_walls*4
    ! WRITE(1) x_walls(:,1:num_walls)
    WRITE(1,"(E12.5)") x_walls(:,1:num_walls)
    CLOSE(1)

    ! WRITE(*,*) "GH 8.11"

    WRITE(filename,"('/flux_downstream_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    WRITE(1,"(I0)") flux_downstream_total
    CLOSE(1)
    WRITE(filename,"('/flux_upstream_total.txt')")
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    WRITE(1,"(I0)") flux_upstream_total
    CLOSE(1)

    ! WRITE(*,*) "GH 8.12"


    ! save miscellaneous data

    call CPU_TIME(t_temp)
    t_total = t_total + (t_temp-t0_total)
    CALL CPU_TIME(t0_total)

    ! WRITE(filename,"('/data.txt')")
    WRITE(filename,"('/data_',I7.7,'.txt')") (ii-1)
    filename = dir_cur(1:dir_cur_length) // filename
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    

    WRITE(1,"(A)") "*t_total"   
    WRITE(1,"(E12.5)") t_total
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*t_BC"   
    WRITE(1,"(E12.5)") t_BC
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*t_collisions"   
    WRITE(1,"(E12.5)") t_collisions
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*t_loop"   
    WRITE(1,"(E12.5)") t_loop
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*n"   
    WRITE(1,"(E12.5)") n
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*ns"   
    WRITE(1,"(E12.5)") ns
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*n_b"   
    WRITE(1,"(E12.5)") n_b
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*ns_b"   
    WRITE(1,"(E12.5)") ns_b
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*Fn"   
    WRITE(1,"(E12.5)") Fn
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*RWF"   
    WRITE(1,"(E12.5)") RWF
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*m_g"   
    WRITE(1,"(E12.5)") m_g
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*nx"   
    WRITE(1,"(I0)") nx
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*ny"   
    WRITE(1,"(I0)") ny
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*tmax"   
    WRITE(1,"(E12.5)") tmax
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*nt"   
    WRITE(1,"(I0)") nt
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*dt"   
    WRITE(1,"(E12.5)") dt
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*dt_to_save"   
    WRITE(1,"(I0)") dt_to_save
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*n_saved"   
    WRITE(1,"(I0)") n_saved
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*include_source"   
    WRITE(1,"(L)") include_source
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*close_inlet"   
    WRITE(1,"(L)") close_inlet
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*include_gun_boundaries"   
    WRITE(1,"(L)") include_gun_boundaries
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*use_homogenous_grid"   
    ! WRITE(1,"(L)") use_homogenous_grid
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*include_collisions"   
    WRITE(1,"(L)") include_collisions
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*dx_0"   
    WRITE(1,"(E12.5)") dx_0
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*dx_factor"   
    WRITE(1,"(E12.5)") dx_factor
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*dy_factor"   
    WRITE(1,"(E12.5)") dy_factor
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*n_collisions"   
    WRITE(1,"(I0)") n_collisions
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*N_added"   
    WRITE(1,"(I0)") SUM(N_added_total)
    WRITE(1,"(A)") "*N_array"   
    WRITE(1,"(I0)") N_array
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*N_simulated"   
    WRITE(1,"(I0)") N_simulated
    WRITE(1,"(A)") "*N_duplicated"   
    WRITE(1,"(I0)") N_duplicated
    WRITE(1,"(A)") "*N_deleted"   
    WRITE(1,"(I0)") N_deleted
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*x_lim"   
    WRITE(1,"(E12.5)") x_lim
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*it_last"   
    WRITE(1,"(I0)") ii
    WRITE(1,"(A)") "*END"
    CLOSE(1)



    n_saved = n_saved + 1







    ! ! save source velocity data
    ! WRITE(filename,"('/vs_',I7.7,'.txt')") (ii-1)
    ! filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! WRITE(1) vs_vec(1:Num_s,:)
    ! CLOSE(1)
    
    ! ! save source position data
    ! WRITE(filename,"('/xs_',I7.7,'.txt')") (ii-1)
    ! filename = dir_cur(1:dir_cur_length) // filename
    ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! WRITE(1,"(E12.5)") xs_vec(1:Num_s,:)
    ! CLOSE(1)



    CALL CPU_TIME(t_temp)
    t_save = t_save + (t_temp-t0_save)

END SUBROUTINE SAVE_DATA







