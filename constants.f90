MODULE constants
IMPLICIT NONE

INTEGER, PARAMETER :: p2=SELECTED_REAL_KIND(13)
INTEGER, PARAMETER :: itime=101,grd_pts=101
REAL(p2), PARAMETER :: x2=4.0,x1=0.0
REAL(p2)            :: dx,dt
REAL(p2), PARAMETER :: zero=0.0_p2,mu2=0.01_p2
REAL(p2), PARAMETER :: one=1.0_p2,mu4=0.001_p2
REAL(p2), PARAMETER :: half=0.5_p2
REAL(p2), PARAMETER :: gamma=1.4_p2
REAL(p2), PARAMETER :: gas_const=286.7_p2
INTEGER             :: i,j
INTEGER             :: ier

REAL(p2),PARAMETER  :: stag_pr=101325.0_p2,stag_temp=300.0_p2, & 
                       exit_pr=84000.0_p2 ,exit_temp=300.0_p2
REAL(p2)  :: initial_rho=exit_pr/(gas_const*exit_temp)
REAL(p2)  :: initial_vel=zero

END MODULE constants
