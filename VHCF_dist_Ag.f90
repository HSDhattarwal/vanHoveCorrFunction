Program dist
Implicit None
  Integer, Parameter                   :: natoms =108,nframes =80000
  Real*8, Parameter                    :: Lx=15.2565, Ly=15.2565, Lz=15.2565
  Real*8, allocatable                  :: xA(:,:), yA(:,:), zA(:,:)
  Real*8, allocatable                  :: xB(:,:), yB(:,:), zB(:,:)
  Character*4                          :: dum1,dum2
Integer, Parameter                     :: nmolAg = 54, nmolI = 54
Integer                                :: maxbinr, moltotal, nbinAA, nbinBB, nbinAB
Real*8                                 :: pi, delr, drAA, rAA, dxAA, dyAA, dzAA
Real*8                                 :: drBB, rBB, dxBB, dyBB, dzBB
Real*8                                 :: drAB, rAB, dxAB, dyAB, dzAB
Integer                                :: i, k, ii, kk, jj, ll, binr, bint, delt
Real*8                                 :: numdensity, ntotalr, ru, rl, r, C, rho
Real*8, Allocatable, Dimension (:,:)     :: histrAA, GrAA
Real*8, Allocatable, Dimension (:,:)     :: histrBB, GrBB
Real*8, Allocatable, Dimension (:,:)     :: histrAB, GrAB
Integer, Dimension (21)                 :: t_lag, countA

pi = 3.1415926536D0
delr = 0.1
maxbinr = int (Lx / (1 * delr)) ! 2 if truncating at L/2
rho = nmolAg/Lx/Ly/Lz

write(*,*) maxbinr

!!!!!!!!!!!Allocate***************************************************************************************
allocate (xA(nframes,nmolAg)); allocate (yA(nframes,nmolAg)); allocate (zA(nframes,nmolAg))
allocate (xB(nframes,nmolI)); allocate (yB(nframes,nmolI)); allocate (zB(nframes,nmolI))

allocate (histrAA (maxbinr,21)); allocate (GrAA (maxbinr,21))
allocate (histrBB (maxbinr,21)); allocate (GrBB (maxbinr,21))
allocate (histrAB (maxbinr,21)); allocate (GrAB (maxbinr,21))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Open(50,file = 'VHCF-dist_Ag.dat')

!----------------read coordinates--------------------------------------

Open(10,file = '../../AgI-MD-pos-750K.xyz')

xA(:,:)=0.0; yA(:,:)=0.0; zA(:,:)=0.0
xB(:,:)=0.0; yB(:,:)=0.0; zB(:,:)=0.0

do k = 1, nframes
   Read(10,*)
   Read(10,*)
   do ii = 1, natoms/2, 2
      Read(10,*)  dum1, xA(k,ii),yA(k,ii), zA(k,ii)
      Read(10,*)  dum1, xA(k,ii+1),yA(k,ii+1), zA(k,ii+1)
      Read(10,*)  dum2, xB(k,ii),yB(k,ii), zB(k,ii)
      Read(10,*)  dum2, xB(k,ii+1),yB(k,ii+1), zB(k,ii+1)
    !  write(20,'(a4,3F20.10)')  dum2, xB(k,ii+1),yB(k,ii+1), zB(k,ii+1)
   enddo
enddo
write(*,*) k
close(10)

Write(*,*) "coordinates Read"

!-------------------calculate distances (r)-----------------
histrAA(:,:) = 0; histrBB(:,:) = 0; histrAB(:,:) = 0
countA(:) = 0

t_lag=(/0.,500.,1000.,1500.,2000.,2500.,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500.,&
      &  8000.,  8500.,  9000.,  9500., 10000./)

do i = 1, 21
  write(*,*) i, t_lag(i)
Do kk = 1, nframes-t_lag(i)
  jj = kk + t_lag(i)
  Do ii = 1, nmolAg
  Do ll = 1, nmolAg
  IF (ii .ne. ll) Then
  countA(i) = countA(i) + 1

	dxAA = xA (kk, ii) - xA (jj, ll)
	dyAA = yA (kk, ii) - yA (jj, ll)
	dzAA = zA (kk, ii) - zA (jj, ll)
        dxAA = dxAA - (nint(dxAA / Lz)) * Lz
        dyAA = dyAA - (nint(dyAA / Lz)) * Lz
        dzAA = dzAA - (nint(dzAA / Lz)) * Lz
	drAA = (dxAA * dxAA) + (dyAA * dyAA) + (dzAA * dzAA)
	rAA = sqrt (drAA)
	binr = (int (rAA / delr)) + 1
	If (binr .LE. maxbinr) Then
	histrAA (binr,i) = histrAA (binr,i) + 1
 	End If

  End If
  End Do
  End Do
End Do
End Do

write(*,*) countA

Do binr = 1, maxbinr
rl = real (binr - 1) * delr
ru = rl + delr
C = 4.0 * pi * (ru**3 - rl**3) / 3.0
r = (binr - 0.5) * delr

Do i = 1, 21
	GrAA (binr,i) = nmolAg*(real (histrAA (binr,i)) / (real (countA(i)*C))) / rho
        Write (i*100, *) r, GrAA (binr,i)
END Do

!Write (50, *) r, GrAA (binr,1), GrAA (binr,2), GrAA (binr,3), GrAA (binr,4), GrAA (binr,5), GrAA (binr,6),  &
!               & GrAA (binr,7), GrAA (binr,8), GrAA (binr,9), GrAA (binr,10), GrAA (binr,11), GrAA (binr,12), &
!               & GrAA (binr,13), GrAA (binr,14), GrAA (binr,15), GrAA (binr,16), GrAA (binr,17), GrAA (binr,18), &
!               & GrAA (binr,19), GrAA (binr,20), GrAA (binr,21)

End Do
End
