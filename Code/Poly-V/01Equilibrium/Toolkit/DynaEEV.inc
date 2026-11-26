!====Get_EEV: compute End-to-End Vector relaxation
SUBROUTINE Get_EEV(Switch)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER Switch, Nt, Nt0, ttel, t, dt, I
	INTEGER, ALLOCATABLE:: NDyn(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, EProd, EEV0, EEV
	DOUBLE PRECISION, ALLOCATABLE:: EX(:), EY(:), EZ(:)
	DOUBLE PRECISION, ALLOCATABLE:: EX0(:,:), EY0(:,:), EZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: SumEEV(:)
	
	SAVE Nt, Nt0, NDyn, ttv0
	SAVE EX, EY, EZ, EX0, EY0, EZ0, SumEEV
!=========================================================
	IF (Switch==1) THEN
!		====sample
		Nt = Nt + 1
		
!		====end-to-end vectors
		DO I = 1, NPoly
			EX(I) = RX(I*NBead) - RX((I-1)*NBead + 1)
			EY(I) = RY(I*NBead) - RY((I-1)*NBead + 1)
			EZ(I) = RZ(I*NBead) - RZ((I-1)*NBead + 1)
		END DO

!		====set new t=0
		IF (MOD(Nt, IT0)==0) THEN
			Nt0 = Nt0 + 1
			ttel = MOD(Nt0 - 1, T0Max) + 1
			ttv0(ttel) = Nt
			DO I = 1, NPoly
				EX0(I, ttel) = EX(I)
				EY0(I, ttel) = EY(I)
				EZ0(I, ttel) = EZ(I)
			END DO
		END IF

		DO t = 1, MIN(Nt0, T0Max)
			dt = Nt - ttv0(t) + 1
			IF (dt>0 .AND. dt<TMax) THEN
				NDyn(dt) = NDyn(dt) + 1
				DO I = 1, NPoly
					EProd = EX(I)*EX0(I, t) + EY(I)*EY0(I, t) + EZ(I)*EZ0(I, t)
					SumEEV(dt) = SumEEV(dt) + EProd
				END DO
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0
		
!		===allocate arrays
		ALLOCATE(NDyn(TMax), ttv0(T0Max))
		ALLOCATE(EX(NPoly), EY(NPoly), EZ(NPoly))
		ALLOCATE(EX0(NPoly, T0Max), EY0(NPoly, T0Max), EZ0(NPoly, T0Max))
		ALLOCATE(SumEEV(TMax))
		
		DO I = 1, TMax
			NDyn(I) = 0
			SumEEV(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate EEV_', A1)") FType
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=81, FILE='EEV_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyn(I)>0) THEN
				EEV = SumEEV(I)/(DBLE(NDyn(I))*DBLE(NPoly))
				IF (I==1) EEV0 = EEV
				WRITE (81, "(10E17.8)") MDTime, EEV/EEV0
			END IF
		END DO
		CLOSE (81)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyn, ttv0)
		DEALLOCATE(EX, EY, EZ)
		DEALLOCATE(EX0, EY0, EZ0)
		DEALLOCATE(SumEEV)
	END IF

RETURN
END
