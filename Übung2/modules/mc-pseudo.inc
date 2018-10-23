!================================================
PROGRAM mc
!================================================

! nstep MC steps 
DO istep=1,nstep 

! randomly displace atoms
   CALL mcmove(x,y,z,n,nacc,e)

ENDDO

!================================================
END PROGRAM mc
!================================================

!================================================
SUBROUTINE mcmove(x,y,z,n,nacc,e)
!================================================

! calculate energy at old positions
CALL energy(x,y,z,n,e)       
eold=e

! random displacement of atoms: 
DO i=1,n   

   CALL random_number(rand)
   xnew(i)=x(i)+delta*(rand-0.5d0)
   CALL random_number(rand) 
   ynew(i)=y(i)+delta*(rand-0.5d0)
   CALL random_number(rand)	
   znew(i)=z(i)+delta*(rand-0.5d0)

ENDDO

! calculate energy at new positions
CALL energy(xnew,ynew,znew,n,e)

! calculate Boltzmann factor: 
! kb = Boltzmann constant, t = temperature
boltz=EXP(-(e-eold)/(kb*t)) 
                                
! accept or reject new positions
CALL random_number(rand) 
IF (rand.LT.boltz) THEN   

   DO i=1,n

      x(i)=xnew(i)
      y(i)=ynew(i)
      z(i)=znew(i)

   ENDDO
         
! count accepted trial moves
   nacc=nacc+1
      
ELSE
         
   e=eold
         
ENDIF
      
RETURN
!================================================
END SUBROUTINE mcmove
!================================================


