!====Get_MSD: compute Mean Squared Displacement and Non-Gaussian Parameter
SUBROUTINE Get_MSD(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, Nt, Nt0, ttel, t, dt, I
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, dx, dy, dz
	DOUBLE PRECISION r2, r4, a2, Sum1, Sum2
	DOUBLE PRECISION r2x, r4x, a2x, Sum1x, Sum2x
	DOUBLE PRECISION r2y, r4y, a2y, Sum1y, Sum2y
	DOUBLE PRECISION r2z, r4z, a2z, Sum1z, Sum2z
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: FR2(:), FR2x(:), FR2y(:), FR2z(:)
	DOUBLE PRECISION, ALLOCATABLE:: FR4(:), FR4x(:), FR4y(:), FR4z(:)
		
	SAVE Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0
	SAVE FR2, FR2x, FR2y, FR2z, FR4, FR4x, FR4y, FR4z
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
				
				Sum1 = 0.0
				Sum1x = 0.0
				Sum1y = 0.0
				Sum1z = 0.0
				Sum2 = 0.0
				Sum2x = 0.0
				Sum2y = 0.0
				Sum2z = 0.0
				DO I = 1, NPart
					dx = MX(I) - RX0(I, t)
					dy = MY(I) - RY0(I, t)
					dz = MZ(I) - RZ0(I, t)
					dx = dx - Rate*(dt-1)*Delt*RY0(I, t)
					r2x = dx**2
					r2y = dy**2
					r2z = dz**2
					r2 = r2x + r2y + r2z
					r4x = r2x**2
					r4y = r2y**2
					r4z = r2z**2
					r4 = r2**2
					Sum1 = Sum1 + r2
					Sum1x = Sum1x + r2x
					Sum1y = Sum1y + r2y
					Sum1z = Sum1z + r2z
					Sum2 = Sum2 + r4
					Sum2x = Sum2x + r4x
					Sum2y = Sum2y + r4y
					Sum2z = Sum2z + r4z
				END DO
				FR2(dt) = FR2(dt) + Sum1
				FR2x(dt) = FR2x(dt) + Sum1x
				FR2y(dt) = FR2y(dt) + Sum1y
				FR2z(dt) = FR2z(dt) + Sum1z
				FR4(dt) = FR4(dt) + Sum2
				FR4x(dt) = FR4x(dt) + Sum2x
				FR4y(dt) = FR4y(dt) + Sum2y 
				FR4z(dt) = FR4z(dt) + Sum2z
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0
		
!		===allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(FR2(TMax), FR2x(TMax), FR2y(TMax), FR2z(TMax))
		ALLOCATE(FR4(TMax), FR4x(TMax), FR4y(TMax), FR4z(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			FR2(I) = 0.0
			FR2x(I) = 0.0
			FR2y(I) = 0.0
			FR2z(I) = 0.0
			FR4(I) = 0.0
			FR4x(I) = 0.0
			FR4y(I) = 0.0
			FR4z(I) = 0.0
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
				r2x = FR2x(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				r2y = FR2y(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				r2z = FR2z(I)/(DBLE(NPart)*DBLE(NDyna(I)))
!				====NGP
				a2 = 0.0
				a2x = 0.0
				a2y = 0.0
				a2z = 0.0
				r4 = FR4(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				r4x = FR4x(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				r4y = FR4y(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				r4z = FR4z(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				IF (r2>0.0) a2 = 0.6*r4/r2**2 - 1
				IF (r2x>0.0) a2x = 1.0/3*r4x/r2x**2 - 1
				IF (r2y>0.0) a2y = 1.0/3*r4y/r2y**2 - 1
				IF (r2z>0.0) a2z = 1.0/3*r4z/r2z**2 - 1
				WRITE (81, "(20E17.8)") MDTime, r2, a2, r2x, a2x, r2y, a2y, r2z, a2z
			END IF
		END DO
		CLOSE (81)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(FR2, FR2x, FR2y, FR2z)
		DEALLOCATE(FR4, FR4x, FR4y, FR4z)
	END IF

RETURN
END
