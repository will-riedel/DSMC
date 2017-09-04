! generate normally distributed random number at each element in 1-D vector y_out (of size n)
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