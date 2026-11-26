!====Get_X4S: compute self part of Four-Point Density Correlation Function
SUBROUTINE Get_X4S(Switch)

	USE SetCOM
	IMPLICIT NONE
	
	DOUBLE PRECISION, PARAMETER:: a2=0.09
	
	INTEGER Switch, Nt, Nt0, ttel, t, dt, I
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, dx, dy, dz, r2, Qs, SumQs
	DOUBLE PRECISION Qs1i, Qs2i, QSS, X4S
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: Qs1(:), Qs2(:)
		
	SAVE Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0, Qs1, Qs2
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
				RX0(I, ttel) = MX(I)
				RY0(I, ttel) = MY(I)
				RZ0(I, ttel) = MZ(I)
			END DO
		END IF

		DO t = 1, MIN(Nt0, T0Max)
			dt = Nt - ttv0(t) + 1
			IF (dt>0 .AND. dt<TMax) THEN
				NDyna(dt) = NDyna(dt) + 1
				
				SumQs = 0
				DO I = 1, NPart
					dx = MX(I) - RX0(I, t)
					dy = MY(I) - RY0(I, t)
					dz = MZ(I) - RZ0(I, t)
					dx = dx - Rate*(dt-1)*Delt*RY0(I, t)
					r2 = dx**2 + dy**2 + dz**2
					IF (r2<a2) THEN
						Qs = 1
					ELSE
						Qs = 0
					END IF
					SumQs = SumQs + Qs
				END DO
				Qs1(dt) = Qs1(dt) + SumQs
				Qs2(dt) = Qs2(dt) + SumQs**2
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0
		
!		====allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(Qs1(TMax), Qs2(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			Qs1(I) = 0.0
			Qs2(I) = 0.0
		END DO

		WRITE (12, "('  To calculate X4S_', A1)") FType
	ELSE IF (Switch==2) THEN
!		====output results to file
		OPEN (UNIT=84, FILE='X4S_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyna(I)>0) THEN
!				====QSS
				QSS = Qs1(I)/(DBLE(NPart)*DBLE(NDyna(I)))
!				====X4S
				Qs1i = Qs1(I)/DBLE(NDyna(I))
				Qs2i = Qs2(I)/DBLE(NDyna(I))
				X4S = Vol*(Qs2i - Qs1i**2)/DBLE(NPart)**2
				WRITE (84, "(10E17.8)") MDTime, QSS, X4S
			END IF
		END DO
		CLOSE (84)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(Qs1, Qs2)
	END IF

RETURN
END
