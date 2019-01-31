!Ein Programm, das erst die Geschwindigkeits-Autokorrelation der Molekülschwingung eines H2O-Moleküls berechnet und dann daraus das Schwingungsspektrum bestimmt
program spektrum
   !Modul zur Berechnung der Autokorrelation
   use autocorrelation
   !modul zur Berechnung der Fourier-Transformation
   use dfourier
   implicit none
   !Geschwindigkeit, Ort, Autokorrelationsfunktionen ohne und mit zeitlicher Mittelung
   real, allocatable :: V(:,:,:), X(:,:,:), Auto(:), Tauto(:)
   !Spektren der Korrelationsfunktionen
   complex, allocatable ::  SpAuto(:), SpTauto(:)
   !Speicher für den Schritt, Zeitschritt
   real :: c, deltat
   !Ein-/Ausgabefehler, Laufindex, Atomzahl, Dimension, Schrittzahl
   integer :: ierror, ierror2, ierror3, i, N, d, steps


   !###########################
   !Festlegen der Parameter
      !Anzahl der Atome
      N=3
      !Dimensionen (üblicherweise 3)
      d=3
      !Anzahl der Zeitschritte
      steps=2000
      !Zeitschritt
      deltat=20.
   !###########################

   !Allozieren der Arrays
   allocate( V(N,d,steps), X(N,d,steps), Auto(steps), Tauto(steps), Spauto(steps), SpTauto(steps) )
   !Öffnen der Dateien
   open(unit=11, file='traj-h2o', status='old', iostat=ierror)
   open(unit=12, file='correlation.txt', status='replace', iostat=ierror2)
   open(unit=13, file='spectrum.txt', status='replace', iostat=ierror3)

   !Einlesen der Werte
   do i=1, steps
      read(11,*,iostat=ierror) c, X(1,:,i), V(1,:,i)
      read(11,*,iostat=ierror) c, X(2,:,i), V(2,:,i)
      read(11,*,iostat=ierror) c, X(3,:,i), V(3,:,i)
      if (ierror/=0) then
         write(*,*) 'Einlesefehler!'
         exit
      end if
   end do

   !Berechnen der Geschwindigkeits-Autokorrelationsfunktionen
   call correlate(V,Auto)
   call timecorrelate(V,Tauto)

   !Ausgabe
   do i=1, steps
      write(12,*,iostat=ierror2) deltat*(i-1), Auto(i), Tauto(i)
      if (ierror3/=0) then
         write(*,*) 'Ausgabefehler!'
         exit
      end if
   end do
 
   !Fouriertransformation der Korrelationsfunktionen
   !dft verlangt complex-arrays
   call dft((1.0,0.0)*Auto,SpAuto)
   call dft((1.0,0.0)*Tauto,SpTauto)
 
   !Ausgabe
   do i=1, steps
      write(13,*,iostat=ierror3) 1./(steps*(137*0.52918e-8*deltat))*(i-1), abs(SpAuto(i)), abs(SpTauto(i))
      if (ierror3/=0) then
         write(*,*) 'Ausgabefehler!'
         exit
      end if
   end do

end program