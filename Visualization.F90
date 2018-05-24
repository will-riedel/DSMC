SUBROUTINE TECPLOT(Timestep, Tot_Time)
    USE CONTAIN
    IMPLICIT NONE
    INTEGER, INTENT(IN):: Timestep
    integer i,j,iel
    integer restart_unit, tecplot_unit
    real*8 tot_time
    character*6 dofno(15)
    character*80 restart_file_name, tecplot_file_name

    data dofno / 'ni','ne','phi','sigma','DOF05', &
         'DOF06','DOF07','DOF08','DOF09','DOF10', &
         'DOF11','DOF12','DOF13','DOF14','DOF15' /
  
! Write the tecplot file
    tecplot_unit = 10
    write(tecplot_file_name,200) 'Output/tecplot',1,'.dat'
200 format(a,i0,a)
     
    if(timestep == 0) then
       open(unit=tecplot_unit,file=tecplot_file_name)
    else
       open(unit=tecplot_unit,file=tecplot_file_name,status='old',position='append')
    endif

! Write the header
    if(timestep == 0) then
       write (tecplot_unit,'(''title='',1h",''Plasma'',1h")')
       if(ndimensions==1)then
          write (tecplot_unit,'(''variables= x'',10(A1,2x,A6))') &
               (',',dofno(i), i=1,3)
       elseif(ndimensions==2)then
          write (tecplot_unit,'(''variables= x, y'',10(A1,2x,A6))') &
               (',',dofno(i), i=1,3)
       endif
    endif

! Write the current zone information

    if(ndimensions==1)then
       write (tecplot_unit,'(''zone  t ='',1h",''Plasma'',1h",'' , I='',i7, &
            '' , DATAPACKING=BLOCK'')') nnodes
    elseif(ndimensions==2)then
       write (tecplot_unit,'(''zone  t ='',1h",''Plasma'',1h",'' , I='',i7, &
            '' , J='',i7, '' , DATAPACKING=BLOCK'')') nnodes_x, nnodes_y
    endif

    write (tecplot_unit,'(''SolutionTime = '',e20.10)') tot_time

! Write out the grid
    if(timestep == 0) then
       write (tecplot_unit,'(6e13.5)') ((x(i),i=1,nnodes_x),j=1,nnodes_y)
       write (tecplot_unit,'(6e13.5)') ((y(j),i=1,nnodes_x),j=1,nnodes_y)
    else
       if(ndimensions==1)then
          write (tecplot_unit,*) 'VARSHARELIST = ([1]=1)'
       elseif(ndimensions==2)then
          write (tecplot_unit,*) 'VARSHARELIST = ([1, 2]=1)'
       endif
    endif
     
! Write the solution
    write (tecplot_unit,'(6e13.5)') & 
         ((n_i_plus(i,j),i=1,nnodes_x), j=1, nnodes_y)
    write (tecplot_unit,'(6e13.5)') & 
         ((n_e_plus(i,j),i=1,nnodes_x), j=1, nnodes_y)
    write (tecplot_unit,'(6e13.5)') & 
         ((phi_plus(i,j),i=1,nnodes_x), j=1, nnodes_y)
! Close the file
    close(tecplot_unit)	
END SUBROUTINE TECPLOT
