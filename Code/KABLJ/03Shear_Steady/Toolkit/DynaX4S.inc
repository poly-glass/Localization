!====Get_X4S: compute self part of Four-Point Density Correlation Function
SUBROUTINE Get_X4S(Switch)

	USE SetCOM
	IMPLICIT NONE
	
	DOUBLE PRECISION, PARAMETER:: a2=0.09
	
	INTEGER Switch, Nt, Nt0, ttel, t, dt, I, Natom(0:2)
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, dx, dy, dz, r2, Qs
	DOUBLE PRECISION SumQs0, Qs1i0, Qs2i0, QSS0, X4S0
	DOUBLE PRECISION SumQs1, Qs1i1, Qs2i1, QSS1, X4S1
	DOUBLE PRECISION SumQs2, Qs1i2, Qs2i2, QSS2, X4S2
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: Qs01(:), Qs02(:)
	DOUBLE PRECISION, ALLOCATABLE:: Qs11(:), Qs12(:)
	DOUBLE PRECISION, ALLOCATABLE:: Qs21(:), Qs22(:)
		
	SAVE Natom, Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0
	SAVE Qs01, Qs02, Qs11, Qs12, Qs21, Qs22
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
				
				SumQs0 = 0
				SumQs1 = 0
				SumQs2 = 0
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
					SumQs0 = SumQs0 + Qs
					IF (PType(I)==1) THEN
						SumQs1 = SumQs1 + Qs
					END IF
					IF (PType(I)==2) THEN
						SumQs2 = SumQs2 + Qs
					END IF
				END DO
				Qs01(dt) = Qs01(dt) + SumQs0
				Qs02(dt) = Qs02(dt) + SumQs0**2
				Qs11(dt) = Qs11(dt) + SumQs1
				Qs12(dt) = Qs12(dt) + SumQs1**2
				Qs21(dt) = Qs21(dt) + SumQs2
				Qs22(dt) = Qs22(dt) + SumQs2**2
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Natom(0) = NPart
		Natom(1) = NPart1
		Natom(2) = NPart2
		
		Nt = 0
		Nt0 = 0
		
!		====allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(Qs01(TMax), Qs02(TMax))
		ALLOCATE(Qs11(TMax), Qs12(TMax))
		ALLOCATE(Qs21(TMax), Qs22(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			Qs01(I) = 0.0
			Qs02(I) = 0.0
			Qs11(I) = 0.0
			Qs12(I) = 0.0
			Qs21(I) = 0.0
			Qs22(I) = 0.0
		END DO

		WRITE (12, "('  To calculate X4S_', A1)") FType
	ELSE IF (Switch==2) THEN
!		====output results to file
		OPEN (UNIT=84, FILE='X4S_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyna(I)>0) THEN
!				====QSS
				QSS0 = Qs01(I)/(DBLE(Natom(0))*DBLE(NDyna(I)))
				QSS1 = Qs11(I)/(DBLE(Natom(1))*DBLE(NDyna(I)))
				
				IF(Natom(2) == 0)THEN
					QSS2 = 0.0
				ELSE
					QSS2 = Qs21(I)/(DBLE(Natom(2))*DBLE(NDyna(I)))				
				END IF	
!				====X4S
				Qs1i0 = Qs01(I)/DBLE(NDyna(I))
				Qs2i0 = Qs02(I)/DBLE(NDyna(I))
				Qs1i1 = Qs11(I)/DBLE(NDyna(I))
				Qs2i1 = Qs12(I)/DBLE(NDyna(I))
				Qs1i2 = Qs21(I)/DBLE(NDyna(I))
				Qs2i2 = Qs22(I)/DBLE(NDyna(I))
				X4S0 = Vol*(Qs2i0 - Qs1i0**2)/DBLE(Natom(0))**2
				X4S1 = Vol*(Qs2i1 - Qs1i1**2)/DBLE(Natom(1))**2
				
				IF(Natom(2) == 0)THEN
					X4S2 = 0.0
				ELSE
				X4S2 = Vol*(Qs2i2 - Qs1i2**2)/DBLE(Natom(2))**2
				END IF
								
				WRITE (84, "(20E17.8)") MDTime, QSS0, QSS1, QSS2, X4S0, X4S1, X4S2
			END IF
		END DO
		CLOSE (84)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(Qs01, Qs02, Qs11, Qs12, Qs21, Qs22)
	END IF

RETURN
END
