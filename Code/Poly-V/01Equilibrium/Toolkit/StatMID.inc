!====Get_MID: compute mean-square internal distance
SUBROUTINE Get_MID(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, NStat, I, J, P, B, B1, B2, n
	
	DOUBLE PRECISION dx, dy, dz, r2, MID	
	DOUBLE PRECISION, ALLOCATABLE:: Acc(:), Sum(:), Rn2(:)
	
	SAVE NStat, Acc, Sum, Rn2
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
		DO n = 1, NBead - 1
			Sum(n) = 0
			Acc(n) = 0
		END DO
	
		DO P = 1, NPoly
			DO B1 = 1, NBead
				I = (P - 1)*NBead + B1
				DO B2 = B1 + 1, NBead
					J = (P - 1)*NBead + B2
					dx = RX(I) - RX(J)
					dy = RY(I) - RY(J)
					dz = RZ(I) - RZ(J)
					r2 = dx**2 + dy**2 + dz**2
					n = B2 - B1
					Sum(n) = Sum(n) + r2
					Acc(n) = Acc(n) + 1
				END DO
			END DO
		END DO
		
		DO n = 1, NBead - 1
			Rn2(n) = Rn2(n) + Sum(n)/Acc(n)
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		NStat = 0
		
!		====allocate arrays
		ALLOCATE(Acc(NBead), Sum(NBead), Rn2(NBead))
		
		DO n = 1, NBead
			Rn2(n) = 0.0
		END DO

		WRITE (12, "('  To calculate MID')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=47, FILE='MID')
		DO n = 1, NBead - 1
			MID = Rn2(n)/(DBLE(NStat)*DBLE(n))
			WRITE (47, "(10E17.8)") DBLE(n), MID
		END DO
		CLOSE (47)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(Acc, Sum, Rn2)
	END IF

RETURN
END
