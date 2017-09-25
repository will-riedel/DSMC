SUBROUTINE SAVE_DATA
    USE CONTAIN
    USE PROPERTIES
    IMPLICIT NONE
    ! CHARACTER(80)::filename

    ! save position data
    WRITE(filename,"('Output/data/x_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) x_vec(1:N_simulated,:)
    ! OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    ! WRITE(1,"(E12.5)") x_vec
    CLOSE(1)

    ! save velocity data
    WRITE(filename,"('Output/data/v_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) v_vec(1:N_simulated,:)
    CLOSE(1)

    ! save cell data (this is probably not required)
    WRITE(filename,"('Output/data/i_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) i_cell_vec(1:N_simulated,:)
    CLOSE(1)

    ! save num_per_cell data? This also seems way unnecessary
    WRITE(filename,"('Output/data/Npc_',I7.7,'.txt')") (ii-1)
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) Npc_slice
    CLOSE(1)

    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*n_collisions_total"   
    ! WRITE(1,"(I0)") n_collisions_total(1:ii)
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*N_candidate_pairs_total"   
    ! WRITE(1,"(I0)") N_candidate_pairs_total(1:ii)
    ! WRITE(1,"(A)") "*"   
    ! WRITE(1,"(A)") "*N_accepted_pairs_total"   
    ! WRITE(1,"(I0)") N_accepted_pairs_total(1:ii)
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*N_added_total"   
    ! WRITE(1,"(I0)") N_added_total(1:ii)
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*ncp_remainder"   
    ! WRITE(1,"(E12.5)") ncp_remainder
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*N_total"   
    ! WRITE(1,"(I0)") N_total(1:ii)

    ! save some tracking vectors
    WRITE(filename,"('Output/data/n_collisions_total.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) n_collisions_total
    CLOSE(1)
    WRITE(filename,"('Output/data/N_candidate_pairs_total.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) N_candidate_pairs_total
    CLOSE(1)
    WRITE(filename,"('Output/data/N_accepted_pairs_total.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) N_accepted_pairs_total
    CLOSE(1)
    WRITE(filename,"('Output/data/N_added_total.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) N_added_total
    CLOSE(1)
    WRITE(filename,"('Output/data/ncp_remainder.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) ncp_remainder
    CLOSE(1)
    WRITE(filename,"('Output/data/N_total.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) N_total
    CLOSE(1)
    WRITE(filename,"('Output/data/x_walls.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    WRITE(1) x_walls
    CLOSE(1)

    ! save miscellaneous data
    WRITE(filename,"('Output/data/data.txt')")
    OPEN(UNIT=1,FILE=filename,FORM="FORMATTED")
    

    WRITE(1,"(A)") "*t_final"   
    WRITE(1,"(E12.5)") t_final
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
    WRITE(1,"(A)") "*Fn"   
    WRITE(1,"(E12.5)") Fn
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
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*use_homogenous_grid"   
    WRITE(1,"(L)") use_homogenous_grid
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
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*n_collisions_total"   
    ! WRITE(1,"(I0)") n_collisions_total(1:ii)
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*N_candidate_pairs_total"   
    ! WRITE(1,"(I0)") N_candidate_pairs_total(1:ii)
    ! WRITE(1,"(A)") "*"   
    ! WRITE(1,"(A)") "*N_accepted_pairs_total"   
    ! WRITE(1,"(I0)") N_accepted_pairs_total(1:ii)
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*N_added_total"   
    ! WRITE(1,"(I0)") N_added_total(1:ii)
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*ncp_remainder"   
    ! WRITE(1,"(E12.5)") ncp_remainder
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*N_total"   
    ! WRITE(1,"(I0)") N_total(1:ii)
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*x_cells_vec"   
    ! WRITE(1,"(E12.5)") x_cells_vec
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*y_cells_vec"   
    ! WRITE(1,"(E12.5)") y_cells_vec
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*x_lim"   
    WRITE(1,"(E12.5)") x_lim
    ! WRITE(1,"(A)") "*"
    ! WRITE(1,"(A)") "*x_walls"   
    ! WRITE(1,"(E12.5)") x_walls
    WRITE(1,"(A)") "*"
    WRITE(1,"(A)") "*it_last"   
    WRITE(1,"(I0)") ii
    WRITE(1,"(A)") "*END"
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








    ! ! save source velocity data
    ! WRITE(filename,"('Output/data/vs_',I7.7,'.txt')") (ii-1)
    ! OPEN(UNIT=1,FILE=filename,FORM="UNFORMATTED")
    ! WRITE(1) vs_vec(1:Num_s,:)
    ! CLOSE(1)






END SUBROUTINE SAVE_DATA







