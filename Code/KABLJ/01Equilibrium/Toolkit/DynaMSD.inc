!====Get_MSD: compute Mean Squared Displacement and Non-Gaussian Parameter
SUBROUTINE Get_MSD(Switch)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Switch, Nt, Nt0, ttel, t, dt, I, Natom(0:2)
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, dx, dy, dz, r2, r4
	DOUBLE PRECISION r02, r04, a02, Sum01, Sum02
	DOUBLE PRECISION r12, r14, a12, Sum11, Sum12
	DOUBLE PRECISION r22, r24, a22, Sum21, Sum22
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: FR02(:), FR04(:)
	DOUBLE PRECISION, ALLOCATABLE:: FR12(:), FR14(:)
	DOUBLE PRECISION, ALLOCATABLE:: FR22(:), FR24(:)
		
	SAVE Natom, Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0
	SAVE FR02, FR04, FR12, FR14, FR22, FR24
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
				
				Sum01 = 0.0
				Sum02 = 0.0
				Sum11 = 0.0
				Sum12 = 0.0
				Sum21 = 0.0
				Sum22 = 0.0

				DO I = 1, NPart
					dx = RX(I) - RX0(I, t)
					dy = RY(I) - RY0(I, t)
					dz = RZ(I) - RZ0(I, t)
					r2 = dx**2 + dy**2 + dz**2
					r4 = r2**2
					Sum01 = Sum01 + r2
					Sum02 = Sum02 + r4
					IF (PType(I)==1) THEN
						Sum11 = Sum11 + r2
						Sum12 = Sum12 + r4
					END IF
					IF (PType(I)==2) THEN
						Sum21 = Sum21 + r2
						Sum22 = Sum22 + r4
					END IF
				END DO
				FR02(dt) = FR02(dt) + Sum01
				FR04(dt) = FR04(dt) + Sum02
				FR12(dt) = FR12(dt) + Sum11
				FR14(dt) = FR14(dt) + Sum12
				FR22(dt) = FR22(dt) + Sum21
				FR24(dt) = FR24(dt) + Sum22
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Natom(0) = NPart
		Natom(1) = NPart1
		Natom(2) = NPart2
		
		Nt = 0
		Nt0 = 0
		
!		===allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(FR02(TMax), FR04(TMax))
		ALLOCATE(FR12(TMax), FR14(TMax))
		ALLOCATE(FR22(TMax), FR24(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			FR02(I) = 0.0
			FR04(I) = 0.0
			FR12(I) = 0.0
			FR14(I) = 0.0
			FR22(I) = 0.0
			FR24(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate MSD_', A1)") FType
	ELSE IF (Switch==2) THEN
!		====output results to file
		OPEN (UNIT=81, FILE='MSD_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyna(I)>0) THEN
!				====MSD
				r02 = FR02(I)/(DBLE(Natom(0))*DBLE(NDyna(I)))
				r12 = FR12(I)/(DBLE(Natom(1))*DBLE(NDyna(I)))
				
				IF(Natom(2) == 0)THEN
					r22 = 0.0
				ELSE
					r22 = FR22(I)/(DBLE(Natom(2))*DBLE(NDyna(I)))				
				END IF
					
!				====NGP
				r04 = FR04(I)/(DBLE(Natom(0))*DBLE(NDyna(I)))
				r14 = FR14(I)/(DBLE(Natom(1))*DBLE(NDyna(I)))

				IF(Natom(2) == 0)THEN
					r24 = 0.0
				ELSE
					r24 = FR24(I)/(DBLE(Natom(2))*DBLE(NDyna(I)))
				END IF
								
				IF (r02>0.0) THEN
					a02 = 0.6*r04/r02**2 - 1
				ELSE
					a02 = 0.0
				END IF
				IF (r12>0.0) THEN
					a12 = 0.6*r14/r12**2 - 1
				ELSE
					a12 = 0.0
				END IF
				IF (r22>0.0) THEN
					a22 = 0.6*r24/r22**2 - 1
				ELSE
					a22 = 0.0
				END IF

				WRITE (81, "(20E17.8)") MDTime, r02, r12, r22, a02, a12, a22
			END IF
		END DO
		CLOSE (81)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(FR02, FR04, FR12, FR14, FR22, FR24)
	END IF

RETURN
END
