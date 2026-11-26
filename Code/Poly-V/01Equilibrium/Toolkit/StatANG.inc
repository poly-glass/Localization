!====Get_ANG: compute bond angle
SUBROUTINE Get_ANG(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, NStat, NAngle, I, J, K, P, B, NBin, Bin
	
	DOUBLE PRECISION Rmax, DelR, PAng
	DOUBLE PRECISION COSTheta, Theta, SumCOSTheta, SumTheta, COSAng, Ang
	DOUBLE PRECISION dxij, dyij, dzij, rij2
	DOUBLE PRECISION dxjk, dyjk, dzjk, rjk2
	DOUBLE PRECISION, ALLOCATABLE:: RBin(:), VBin(:), PTheta(:)
	
	SAVE NStat, NAngle, SumCOSTheta, SumTheta, NBin, DelR, RBin, VBin, PTheta
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
!		====distribution of square bond length
		DO P = 1, NPoly
			DO B = 1, NBead - 2
				I = (P - 1)*NBead + B
				J = I + 1
				K = I + 2
				dxij = RX(I) - RX(J)
				dyij = RY(I) - RY(J)
				dzij = RZ(I) - RZ(J)
				dxjk = RX(J) - RX(K)
				dyjk = RY(J) - RY(K)
				dzjk = RZ(J) - RZ(K)
				rij2 = dxij**2 + dyij**2 + dzij**2
				rjk2 = dxjk**2 + dyjk**2 + dzjk**2
				COSTheta = (dxij*dxjk + dyij*dyjk + dzij*dzjk)/(rij2*rjk2)**0.5
				Theta = ACOS(COSTheta)
				SumCOSTheta = SumCOSTheta + COSTheta
				SumTheta = SumTheta + Theta
				Bin = INT(Theta/DelR)
				IF (Bin<=NBin) PTheta(Bin) = PTheta(Bin) + 1
			END DO
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		NStat = 0
		
!		====total angles
		NAngle = (NBead - 2)*NPoly

!		====set basic info
		Rmax = PI
		DelR = Rmax/1.0D2
		
!		====number of bins
		NBin = INT(Rmax/DelR)
		
!		====allocate arrays
		ALLOCATE(RBin(0:NBin), VBin(0:NBin), PTheta(0:NBin))
		
!		====bin position and shell volume
		DO I = 0, NBin
			RBin(I) = (I + 0.5D0)*DelR
			VBin(I) = 4*PI*DelR**3*(DBLE(I+1)**3 - DBLE(I)**3)/3
			PTheta(I) = 0.0
		END DO
		
		SumTheta = 0.0
		SumCOSTheta = 0.0

		WRITE (12, "('  To calculate ANG')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=55, FILE='PANG')
		DO I = 0, NBin - 1
			PAng = 4*PI*RBin(I)**2*PTheta(I)/(VBin(I)*DBLE(NStat)*DBLE(NAngle))
			WRITE (55, "(10E17.8)") RBin(I), PAng
		END DO
		CLOSE (55)
		
		COSAng = SumCOSTheta/(DBLE(NStat)*DBLE(NAngle))
		Ang = SumTheta/(DBLE(NStat)*DBLE(NAngle))
		OPEN (UNIT=56, FILE='ANG')
		WRITE (56, "(10E17.8)") COSAng, Ang
		CLOSE (56)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(RBin, VBin, PTheta)
	END IF

RETURN
END
