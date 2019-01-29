!Modul zur Berechnung der Autokorrelation
module autocorrelation
contains

!Subroutine zur Berechnung der Autokorrelation
subroutine correlate(V,A)
   implicit none
   !(Atom, Koordinate, Wert), Korrelationsfunktion, Summe zur Berechnung von A
   real :: V(:,:,:), A(:), summe
   !Atomzahl, Schrittzahl, Dimension, Laufindices, Form von V
   integer :: N, steps, d, i, j, t, t0, s(3)
  
   !Festlegen der Dimensionen 
   s=shape(V)
   !Anzahl Atome
   N=s(1)
   !Anzahl der Schritte
   steps=s(3)
   !Dimension (sollte 3 sein)
   d=s(2)

   !Berechnung der Autokorrelationsfunktion für jeden Schritt
   do t=1, steps
      summe=0.
      !Mittelung über alle Atome
      do i=1, N
         !Skalarprodukt
         do j=1, d
            summe=summe+V(i,j,t)*V(i,j,1)
         end do
      end do
      !Teilen durch Atomzahl
      A(t)=1./real(N)*summe
   end do

end subroutine

!Autokorrelation mit Mittelung über die Zeitschritte
subroutine timecorrelate(V,A)
   implicit none
   !(Atom, Koordinate, Wert), Summe zur Berechnung von A
   real :: V(:,:,:), A(:), summe
   !Atomzahl, Schrittzahl, Dimension, Laufindices, Form von V
   integer :: N, steps, d, i, j, t, t0, s(3)
   
   !Festlegen der Dimensionen
   s=shape(V)
   !Anzahl Atome
   N=s(1)
   !Anzahl der Schritte
   steps=s(3)
   !Dimension (sollte 3 sein)
   d=s(2)

   !Berechnung der Autokorrelationsfunktion für jeden Zeitschritt
   do t=1, steps
      summe=0.
      !Ändern des Zeitnullpunkts, +1, da unten -1 und alle steps erreicht werden sollen
      do t0=1, (steps-t+1)
         !Mittelung über alle Atome
         do i=1, N
            !Skalarprodukt
            do j=1, d
               !Verschobene Nullpunkte, -1 damit auch der Schritt ohne Verschiebung mitgenommen wird
               summe=summe+V(i,j,t+t0-1)*V(i,j,t0)
            end do
         end do
      end do
      !Teilen durch Atomzahl und Anzahl der Zeitschritte, über die gemittelt wurde
      A(t)=1./(real(N)*real(steps-t+1))*summe
   end do

end subroutine

end module
