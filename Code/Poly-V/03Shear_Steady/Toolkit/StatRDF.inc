!====Get_RDF: compute Radial Distribution Function
SUBROUTINE Get_RDF(Switch)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER Switch, NStat, I, J, NBin, Bin
	
	DOUBLE PRECISION Rmax, DelR, dx, dy, dz, r2, RDF
	DOUBLE PRECISION, ALLOCATABLE:: RBin(:), VBin(:), Gbb(:)

	SAVE NStat, NBin, DelR, RBin, VBin, Gbb
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
		DO I = 1, NPart
			DO J = I + 1, NPart
				dx = MX(I) - MX(J)
				dy = MY(I) - MY(J)
				dz = MZ(I) - MZ(J)
				dx = dx - DelX*ANINT(dy/Box)
				dx = dx - Box*ANINT(dx/Box)
				dy = dy - Box*ANINT(dy/Box)
				dz = dz - Box*ANINT(dz/Box)
				r2 = dx**2 + dy**2 + dz**2
				Bin = INT(r2**0.5/DelR)
				IF (Bin<=NBin) Gbb(Bin) = Gbb(Bin) + 1
			END DO
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		NStat = 0

!		====set basic info
		Rmax = 0.50D0*Box
		DelR = 0.01D0
		
!		====number of bins
		NBin = INT(Rmax/DelR)
		
!		====allocate arrays
		ALLOCATE(RBin(0:NBin), VBin(0:NBin), Gbb(0:NBin))

!		====bin position and shell volume
		DO I = 0, NBin
			RBin(I) = (I + 0.5D0)*DelR
			VBin(I) = 4*PI*DelR**3*(DBLE(I+1)**3 - DBLE(I)**3)/3
			Gbb(I) = 0.0
		END DO

		WRITE (12, "('  To calculate RDF')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=25, FILE='RDF')
		DO I = 0, NBin - 1
			RDF = 2*Vol*Gbb(I)/(VBin(I)*DBLE(NStat)*NPart**2)
			WRITE (25, "(10E17.8)") RBin(I), RDF
		END DO
		CLOSE (25)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(RBin, VBin, Gbb)
	END IF

RETURN
END
