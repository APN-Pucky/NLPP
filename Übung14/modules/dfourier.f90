!Ein Modul zur diskreten Fourier-Transformation
module dfourier

!Pi
real, parameter :: PI=2*acos(0.)

contains

subroutine dft(Y,S)
   implicit none
   !Input: Zeitsignal Y, Ausgabe: Spektrum S
   complex :: Y(:), S(:)
   !Anzahl der Punkte, Laufindices
   integer :: N, i, j

   !Festlegen der Anzahl der Datenpunkte
   N=size(Y)

   !DFT
   do i=1, N
      S(i)=(0.0, 0.0)
      do j=1, N
         S(i)=S(i)+Y(j)*exp(-2*PI*(0.0,1.0)*(i-1)*(j-1)/N)
      end do
   end do

end subroutine

!Inverse Fouriertransformation
subroutine invdft(Y,S)
   implicit none
   !Input: Zeitsignal Y, Ausgabe: Spektrum S
   complex :: Y(:), S(:)
   !Anzahl der Punkte, Laufindices
   integer :: N, i, j

   !Festlegen der Anzahl der Datenpunkte
   N=size(Y)

   !Rücktransformation
   do i=1, N
      S(i)=(0.0, 0.0)
      do j=1, N
         S(i)=S(i)+Y(j)*exp(2*PI*(0.0,1.0)*(i-1)*(j-1)/N)
      end do
      S(i)=1./N*S(i)
   end do

end subroutine

end module