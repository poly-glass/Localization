!====Get_MSD: compute Mean Squared Displacement and Non-Gaussian Parameter
SUBROUTINE Get_MSD(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, Nt, Nt0, ttel, t, dt, I
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, dx, dy, dz, r2, r4, a2, Sum1, Sum2
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: FR2(:), FR4(:)
		
	SAVE Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0, FR2, FR4
!=========================================================
	IF (Switch==1) THEN
!		====sample
		Nt = Nt + 1

!		====set new t=0
		IF (MOD(Nt, IT0)==0) THEN
			Nt0 = Nt0 + 1
			ttel = MOD(Nt0 - 1, T0Max) + 1
			ttv0(ttel) = Nt
			
			DO I = 1, NPart
				RX0(I, ttel) = RX(I)
				RY0(I, ttel) = RY(I)
				RZ0(I, ttel) = RZ(I)
			END DO
		END IF

		DO t = 1, MIN(Nt0, T0Max)
			dt = Nt - ttv0(t) + 1
			IF (dt>0 .AND. dt<TMax) THEN
				NDyna(dt) = NDyna(dt) + 1
				
				Sum1 = 0.0
				Sum2 = 0.0
				DO I = 1, NPart
					dx = RX(I) - RX0(I, t)
					dy = RY(I) - RY0(I, t)
					dz = RZ(I) - RZ0(I, t)
					r2 = dx**2 + dy**2 + dz**2
					r4 = r2**2
					Sum1 = Sum1 + r2
					Sum2 = Sum2 + r4
				END DO
				FR2(dt) = FR2(dt) + Sum1
				FR4(dt) = FR4(dt) + Sum2
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0
		
!		===allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(FR2(TMax), FR4(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			FR2(I) = 0.0
			FR4(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate MSD_', A1)") FType
	ELSE IF (Switch==2) THEN
!		====output results to file
		OPEN (UNIT=81, FILE='MSD_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyna(I)>0) THEN
!				====MSD
				r2 = FR2(I)/(DBLE(NPart)*DBLE(NDyna(I)))
!				====NGP
				r4 = FR4(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				a2 = 0.0
				IF (r2>0.0) a2 = 0.6*r4/r2**2 - 1
				WRITE (81, "(10E17.8)") MDTime, r2, a2
			END IF
		END DO
		CLOSE (81)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(FR2, FR4)
	END IF

RETURN
END
