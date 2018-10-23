!================================================
SUBROUTINE energy(x,y,z,n,e)
!================================================
IMPLICIT NONE
INTEGER :: n,i,j
REAL*8 :: x(n),y(n),z(n),e,e2,r,dx,dy,dz
REAL*8, PARAMETER :: a=4.08,c=34.408,eps=2.04968D-21 !IN JOULE

e=0.0d0

DO i=1,n

   e2=0.0d0

   DO j=1,n

      IF(i.NE.j)THEN

         dx=x(i)-x(j)
         dy=y(i)-y(j)
         dz=z(i)-z(j)

         r=SQRT(dx**2+dy**2+dz**2)

         e=e+0.5d0*(a/r)**10

         e2=e2+(a/r)**8

      ENDIF
   ENDDO
   e=e-c*SQRT(e2)
ENDDO
! Umrechnen der Energie in Joule
e=e*eps
RETURN
!================================================
END SUBROUTINE energy
!================================================

