!====Get_PSP: compute Polymer-Specific Properties
SUBROUTINE Get_PSP(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, NStat, I, J, P, B
	
	DOUBLE PRECISION dx, dy, dz, r2, L2, SumRe2, SumRg2
	DOUBLE PRECISION BCF, RE2, RG2
	DOUBLE PRECISION, ALLOCATABLE:: Fbb(:)
	DOUBLE PRECISION, ALLOCATABLE:: LX(:), LY(:), LZ(:)
	DOUBLE PRECISION, ALLOCATABLE:: UX(:), UY(:), UZ(:)
	
	SAVE NStat, Fbb, SumRe2, SumRg2
	SAVE LX, LY, LZ, UX, UY, UZ
!=========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
		DO P = 1, NPoly
!			====set bond vectors and unit bond vectors
			DO B = 1, NBead - 1
				I = (P - 1)*NBead + B
				LX(B) = RX(I+1) - RX(I)
				LY(B) = RY(I+1) - RY(I)
				LZ(B) = RZ(I+1) - RZ(I)
				L2 = LX(B)**2 + LY(B)**2 + LZ(B)**2
				UX(B) = LX(B)/L2**0.5
				UY(B) = LY(B)/L2**0.5
				UZ(B) = LZ(B)/L2**0.5
			END DO
			
!			====end-to-end distance and bond correlation function
			DO B = 1, NBead - 1
				Fbb(B) = Fbb(B) + UX(B)*UX(1) + UY(B)*UY(1) + UZ(B)*UZ(1)
			END DO
			
!			====end-to-end distance
			dx = RX(P*NBead) - RX((P-1)*NBead + 1)
			dy = RY(P*NBead) - RY((P-1)*NBead + 1)
			dz = RZ(P*NBead) - RZ((P-1)*NBead + 1)
			r2 = dx**2 + dy**2 + dz**2
			SumRe2 = SumRe2 + r2
			
!			====radius of gyration
			DO B = 1, NBead
				I = (P - 1)*NBead + B
				dx = RX(I) - PX(P)
				dy = RY(I) - PY(P)
				dz = RZ(I) - PZ(P)
				r2 = dx**2 + dy**2 + dz**2
				SumRg2 = SumRg2 + r2
			END DO
		END DO
		
	ELSE IF (Switch==0) THEN
!		====initialization
		NStat = 0
		
		SumRe2 = 0.0
		SumRg2 = 0.0

!		====allocate arrays
		ALLOCATE(Fbb(NBead))
		
		DO I = 1, NBead
			Fbb(I) = 0.0
		END DO
		
!		====allocate arrays
		ALLOCATE(LX(NBead), LY(NBead), LZ(NBead))
		ALLOCATE(UX(NBead), UY(NBead), UZ(NBead))
		
		WRITE (12, "('  To calculate PSP')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=26, FILE='BCF')
		DO I = 1, NBead - 1
			BCF = Fbb(I)/(DBLE(NStat)*DBLE(NPoly))
			WRITE (26, "(10E17.8)") DBLE(I - 1), BCF
		END DO
		CLOSE (26)
		
		RE2 = SumRe2/(DBLE(NStat)*DBLE(NPoly))
		RG2 = SumRg2/(DBLE(NStat)*DBLE(NPoly)*DBLE(NBead))
		
		OPEN (UNIT=27, FILE='Size')
		WRITE (27, "(10E17.8)") RE2, RG2
		CLOSE (27)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(Fbb)
		DEALLOCATE(LX, LY, LZ)
		DEALLOCATE(UX, UY, UZ)
	END IF

RETURN
END
