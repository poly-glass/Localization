!====Get_LB2: compute mean-square bond length
SUBROUTINE Get_LB2(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, NStat, NBond, I, P, B, NBin, Bin
	
	DOUBLE PRECISION Rmax, DelR, dx, dy, dz, r2
	DOUBLE PRECISION SumR2, LB2, PLB2
	DOUBLE PRECISION, ALLOCATABLE:: RBin(:), VBin(:), PR2(:)
	
	SAVE NStat, NBond, SumR2, NBin, DelR, RBin, VBin, PR2
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
!		====distribution of square bond length
		DO P = 1, NPoly
			DO B = 1, NBead - 1
				I = (P - 1)*NBead + B
				dx = RX(I+1) - RX(I)
				dy = RY(I+1) - RY(I)
				dz = RZ(I+1) - RZ(I)
				r2 = dx**2 + dy**2 + dz**2
				SumR2 = SumR2 + r2
				Bin = INT(r2/DelR)
				IF (Bin<=NBin) PR2(Bin) = PR2(Bin) + 1
			END DO
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		NStat = 0
		
!		====total bonds
		NBond = (NBead - 1)*NPoly

!		====set basic info
		Rmax = 2.0D0
		DelR = Rmax/2.0D2
		
!		====number of bins
		NBin = INT(Rmax/DelR)
		
!		====allocate arrays
		ALLOCATE(RBin(0:NBin), VBin(0:NBin), PR2(0:NBin))
		
!		====bin position and shell volume
		DO I = 0, NBin
			RBin(I) = (I + 0.5D0)*DelR
			VBin(I) = 4*PI*DelR**3*(DBLE(I+1)**3 - DBLE(I)**3)/3
			PR2(I) = 0.0
		END DO
		
		SumR2 = 0.0

		WRITE (12, "('  To calculate LB2')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=45, FILE='PLB2')
		DO I = 0, NBin - 1
			PLB2 = 4*PI*RBin(I)**2*PR2(I)/(VBin(I)*DBLE(NStat)*DBLE(NBond))
			WRITE (45, "(10E17.8)") RBin(I), PLB2
		END DO
		CLOSE (45)
		
		LB2 = SumR2/(DBLE(NStat)*DBLE(NBond))
		OPEN (UNIT=46, FILE='LB2')
		WRITE (46, "(10E17.8)") LB2
		CLOSE (46)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(RBin, VBin, PR2)
	END IF

RETURN
END
