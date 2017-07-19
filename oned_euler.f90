PROGRAM oned_euler
USE constants
USE variables
IMPLICIT NONE
REAL(p2), DIMENSION(0:grd_pts-1) ::xc
TYPE(cell_data)::cell(3)


dx=(x2-x1)/REAL(grd_pts-1)
dt=dx

DO i=0,grd_pts-1
   xc(i)= REAL(i)*dx
END DO
!DO i=0,grd_pts-1
!WRITE(*,*) xc(i)
!END DO


CALL ic(cell)

OPEN(UNIT=100,FILE='p1ic.DAT', STATUS='REPLACE',ACTION='WRITE',IOSTAT=ier)
DO i=0,itime
WRITE(100,*) (cell(1)%p1(i,j),j=0,grd_pts-1)
END DO


CALL p2qic(cell)

OPEN(UNIT=101,FILE='q1ic.DAT', STATUS='REPLACE',ACTION='WRITE',IOSTAT=ier)
DO i=0,itime
WRITE(101,*) (cell(2)%q1(i,j),j=0,grd_pts-1)
END DO


CALL q2eic(cell)

OPEN(UNIT=122,FILE='e1ic.DAT', STATUS='REPLACE',ACTION='WRITE',IOSTAT=ier)
DO i=0,itime
WRITE(122,*) (cell(3)%e1(i,j),j=0,grd_pts-1)
END DO



CALL ftcs(cell)

OPEN(UNIT=103,FILE='p1fc.DAT', STATUS='REPLACE',ACTION='WRITE',IOSTAT=ier)
DO i=0,itime
WRITE(103,*) (cell(1)%p1(i,j),j=0,grd_pts-1)
END DO





END PROGRAM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ic(cell1)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell1(3)
DO i=0,itime
 DO j=0,grd_pts-1
   IF(i==0) THEN
   cell1(1)%p1(i,j) = initial_rho
   cell1(1)%p2(i,j) = initial_vel
   cell1(1)%p3(i,j) = exit_pr
   cell1(1)%p4(i,j) = exit_temp
 ELSE
   cell1(1)%p1(i,j) = zero
   cell1(1)%p2(i,j)= zero
   cell1(1)%p3(i,j) = zero
   cell1(1)%p4(i,j)= zero
 END IF
 END DO
END DO
END SUBROUTINE


SUBROUTINE p2qic(cell2)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data),DIMENSION(3),INTENT(INOUT):: cell2
DO i=0,itime
 DO j=0,grd_pts-1
 IF(i==0) THEN
  cell2(2)%q1(i,j)=cell2(1)%p1(i,j)
  cell2(2)%q2(i,j)=cell2(1)%p1(i,j)*cell2(1)%p2(i,j)
  cell2(2)%q3(i,j)=cell2(1)%p4(i,j)/(gamma-one)+(half*cell2(1)%p1(i,j)*cell2(1)%p2(i,j)*cell2(1)%p2(i,j))
ELSE
  cell2(2)%q1(i,j)=zero
  cell2(2)%q2(i,j)=zero
  cell2(2)%q3(i,j)=zero
END IF
 END DO
END DO
END SUBROUTINE


SUBROUTINE q2eic(cell3)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell3(3)
DO i=0,itime 
 DO j=0,grd_pts-1 
   IF(i==0) THEN
  cell3(3)%e1(i,j)=cell3(2)%q2(i,j)
  cell3(3)%e2(i,j)=(gamma-one)*cell3(2)%q3(i,j)+(half*(3.0_p2-gamma)*cell3(2)%q2(i,j)*cell3(2)%q2(i,j)/cell3(2)%q1(i,j))
  cell3(3)%e3(i,j)=(gamma*cell3(2)%q3(i,j)-(half*(gamma-one)*cell3(2)%q2(i,j)*cell3(2)%q2(i,j)/cell3(2)%q1(i,j))) &
                                                    *(cell3(2)%q2(i,j)/cell3(2)%q1(i,j))
 ELSE 
  cell3(3)%e1(i,j)=zero
  cell3(3)%e2(i,j)=zero
  cell3(3)%e3(i,j)=zero
  END IF
 END DO
END DO
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE ftcs(cell4)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell4(3)
DO i=1,itime
 DO j=1,grd_pts-2
   cell4(2)%q1(i,j)=cell4(2)%q1(i-1,j)-(half*dt/dx*(cell4(3)%e1(i-1,j+1)-cell4(3)%e1(i-1,j-1)))
   cell4(2)%q2(i,j)=cell4(2)%q2(i-1,j)-(half*dt/dx*(cell4(3)%e2(i-1,j+1)-cell4(3)%e2(i-1,j-1))) 
   cell4(2)%q3(i,j)=cell4(2)%q3(i-1,j)-(half*dt/dx*(cell4(3)%e3(i-1,j+1)-cell4(3)%e3(i-1,j-1)))
!Dissipation term to make it stable
!!&+mu2*((cell4(2)%q1(i-1,j+1))-(2*(cell4(2)%q1(i-1,j)))-cell4(2)%q1(i-1,j-1)) &
!!+mu4*((-4*cell4(2)%q1(i-1,j-1))+(6*cell4(2)%q1(i-1,j))-(4*cell4(2)%q1(i-1,j+1))+cell4(2)%q1(i-1,j-2)+cell4(2)%q1(i-1,j+2))
   
   CALL q2p(cell4,i-1,j)
   CALL q2e(cell4,i-1,j)
 END  DO
   CALL ic1(cell4,i)
   CALL p2q1(cell4,i)
   CALL q2e1(cell4,i)
END DO
  
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


SUBROUTINE q2p(cell5,k,l)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell5(3)
INTEGER, INTENT(IN):: k,l
cell5(1)%p1(k+1,l)=cell5(2)%q1(k+1,l)
cell5(1)%p2(k+1,l)=cell5(2)%q2(k+1,l)/cell5(2)%q1(k+1,l)
cell5(1)%p3(k+1,l)=(gamma-1)*(cell5(2)%q3(k+1,l)-((half*cell5(2)%q2(k+1,l)*cell5(2)%q2(k+1,l))/cell5(2)%q1(k+1,l)))
cell5(1)%p4(k+1,l)=cell5(1)%p3(k+1,l)/(gas_const*cell5(1)%p1(k+1,l))
END SUBROUTINE

SUBROUTINE q2e(cell8,m,n)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell8(3)
INTEGER, INTENT(IN):: m,n
cell8(3)%e1(m+1,n)=cell8(2)%q2(m+1,n)
cell8(3)%e2(m+1,n)=((gamma-one)*cell8(2)%q3(m+1,n))+&
                                ((half*(3.0_p2-gamma)*cell8(2)%q2(m+1,n)*cell8(2)%q2(m+1,n))/cell8(2)%q1(m+1,n))
cell8(3)%e3(m+1,n)=((gamma*cell8(2)%q3(m+1,n))-&
                                ((half*(gamma-one)*cell8(2)%q2(m+1,n)*cell8(2)%q2(m+1,n))/cell8(2)%q1(m+1,n))) &
                                                      *((cell8(2)%q2(m+1,n)/cell8(2)%q1(m+1,n)))
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ic1(cell6,l)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell6(3)
INTEGER, INTENT(IN):: l
  cell6(1)%p2(l,0)=cell6(1)%p2(l,1)
  cell6(1)%p4(l,0)=stag_temp-(((gamma-1)*half*cell6(1)%p2(l,1)*cell6(1)%p2(l,1))/(gamma*gas_const))
  cell6(1)%p3(l,0)=stag_pr*((cell6(1)%p4(l,0)/stag_temp)**(gamma/(gamma-one)))
  cell6(1)%p1(l,0)=cell6(1)%p3(l,0)/(gas_const*cell6(1)%p4(l,0))

  cell6(1)%p2(l,grd_pts-1)=cell6(1)%p2(l,grd_pts-2)
  cell6(1)%p4(l,grd_pts-1)=stag_temp-((gamma-1)*half*cell6(1)%p2(l,grd_pts-2)*cell6(1)%p2(l,grd_pts-2)/(gamma*gas_const))
  cell6(1)%p3(l,grd_pts-1)=stag_pr*((cell6(1)%p4(l,grd_pts-1)/stag_temp)**(gamma/(gamma-one)))
  cell6(1)%p1(l,grd_pts-1)=(cell6(1)%p3(l,grd_pts-1))/(gas_const*cell6(1)%p4(l,grd_pts-1))

END SUBROUTINE

SUBROUTINE p2q1(cell7,r)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell7(3)
INTEGER, INTENT(IN):: r
  cell7(2)%q1(r,0)=cell7(1)%p1(r,0)
  cell7(2)%q2(r,0)=cell7(1)%p1(r,0)*cell7(1)%p2(r,0)
  cell7(2)%q3(r,0)=cell7(1)%p4(r,0)/(gamma-one)+(half*cell7(1)%p1(r,0)*cell7(1)%p2(r,0)*cell7(1)%p2(r,0))
  cell7(2)%q1(r,grd_pts-1)=cell7(1)%p1(r,grd_pts-1)
  cell7(2)%q2(r,grd_pts-1)=cell7(1)%p1(r,grd_pts-1)*cell7(1)%p2(r,grd_pts-1)
  cell7(2)%q3(r,grd_pts-1)=(cell7(1)%p4(r,grd_pts-1)/(gamma-one))+(half*cell7(1)%p1(r,grd_pts-1)&
                                                       *cell7(1)%p2(r,grd_pts-1)*cell7(1)%p2(r,grd_pts-1))
END SUBROUTINE

SUBROUTINE q2e1(cell9,p)
USE variables
USE constants
IMPLICIT NONE
TYPE(cell_data), INTENT(INOUT):: cell9(3)
INTEGER, INTENT(IN):: p
  cell9(3)%e1(p,0)=cell9(2)%q2(p,0)
  cell9(3)%e2(p,0)=(gamma-one)*cell9(2)%q3(p,0)+&
                     (half*(3.0_p2-gamma)*cell9(2)%q2(p,0)*cell9(2)%q2(p,0)/cell9(2)%q1(p,0))
  cell9(3)%e3(p,0)=((gamma*cell9(2)%q3(p,0))-(half*(gamma-one)*cell9(2)%q2(p,0)/cell9(2)%q1(p,0))) &
                                                              *(cell9(2)%q2(p,0)/cell9(2)%q1(p,0))
  cell9(3)%e1(p,grd_pts-1)=cell9(2)%q2(p,grd_pts-1)
  cell9(3)%e2(p,grd_pts-1)=(gamma-one)*cell9(2)%q3(p,grd_pts-1)+&
                     (half*(3.0_p2-gamma)*cell9(2)%q2(p,grd_pts-1)*cell9(2)%q2(p,grd_pts-1)/cell9(2)%q1(p,grd_pts-1))
  cell9(3)%e3(p,grd_pts-1)=((gamma*cell9(2)%q3(p,grd_pts-1))-(half*(gamma-one)*cell9(2)%q2(p,grd_pts-1)/cell9(2)%q1(p,grd_pts-1)))&
                                                              *(cell9(2)%q2(p,grd_pts-1)/cell9(2)%q1(p,grd_pts-1))
END SUBROUTINE
