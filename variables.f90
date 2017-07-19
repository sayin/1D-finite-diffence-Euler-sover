MODULE variables
USE constants

IMPLICIT NONE

TYPE :: cell_data
REAL(p2), DIMENSION(0:itime,0:grd_pts-1) :: p1=zero,p2=zero,p3=zero,p4=zero !prmitive variable= [p1=rho,p2=u,p3=p,p4=T]
REAL(p2), DIMENSION(0:itime,0:grd_pts-1) :: q1=zero,q2=zero,q3=zero !conservative variables=[rho,rho*u,rho*e]
REAL(p2), DIMENSION(0:itime,0:grd_pts-1) :: e1=zero,e2=zero,e3=zero  !flux vector=[rho*u,rho*u^2+p,(rho*e+p)*u]
END TYPE cell_data
!TYPE(cell_data)::cell(3)
END MODULE variables
