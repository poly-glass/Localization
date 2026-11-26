!====SetFSQ: global variables and arrays for FSQ
MODULE SetFSQ
	
	INTEGER, PARAMETER:: QIMax=20
	
	INTEGER Qn
	
	DOUBLE PRECISION Qone
	DOUBLE PRECISION, ALLOCATABLE:: Qx(:), Qy(:), Qz(:)
	
END MODULE

!====Get_FSQ: compute Self-Intermediate Scattering Function Fs(q,t)
SUBROUTINE Get_FSQ(Switch)
	
	USE SetCOM
	USE SetFSQ
	IMPLICIT NONE
	
	INTEGER Switch, Nt, Nt0, ttel, t, dt, I, K
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, dx, dy, dz
	DOUBLE PRECISION SetCOS, SumCOS, FSQ, XSQ, FSQ1, FSQ2
	DOUBLE PRECISION SetCOSx, SumCOSx, FSQx, XSQx, FSQ1x, FSQ2x
	DOUBLE PRECISION SetCOSy, SumCOSy, FSQy, XSQy, FSQ1y, FSQ2y
	DOUBLE PRECISION SetCOSz, SumCOSz, FSQz, XSQz, FSQ1z, FSQ2z
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS1(:), FCOS1x(:), FCOS1y(:), FCOS1z(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS2(:), FCOS2x(:), FCOS2y(:), FCOS2z(:)
	
	SAVE Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0
	SAVE FCOS1, FCOS1x, FCOS1y, FCOS1z, FCOS2, FCOS2x, FCOS2y, FCOS2z
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
				
				SumCOS = 0.0
				SumCOSx = 0.0
				SumCOSy = 0.0
				SumCOSz = 0.0
				DO I = 1, NPart
					dx = MX(I) - RX0(I, t)
					dy = MY(I) - RY0(I, t)
					dz = MZ(I) - RZ0(I, t)
					dx = dx - Rate*(dt-1)*Delt*RY0(I, t)
					DO K = 1, Qn
						SetCOS = COS(Qx(K)*dx + Qy(K)*dy + Qz(K)*dz)
						SumCOS = SumCOS + SetCOS
					END DO
					SetCOSx = COS(Qone*dx)
					SetCOSy = COS(Qone*dy)
					SetCOSz = COS(Qone*dz)
					SumCOSx = SumCOSx + SetCOSx
					SumCOSy = SumCOSy + SetCOSy
					SumCOSz = SumCOSz + SetCOSz
				END DO
				FCOS1(dt) = FCOS1(dt) + (SumCOS/Qn)
				FCOS1x(dt) = FCOS1x(dt) + (SumCOSx)
				FCOS1y(dt) = FCOS1y(dt) + (SumCOSy)
				FCOS1z(dt) = FCOS1z(dt) + (SumCOSz)
				FCOS2(dt) = FCOS2(dt) + (SumCOS/Qn)**2
				FCOS2x(dt) = FCOS2x(dt) + (SumCOSx)**2
				FCOS2y(dt) = FCOS2y(dt) + (SumCOSy)**2
				FCOS2z(dt) = FCOS2z(dt) + (SumCOSz)**2
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0

!		===allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(FCOS1(TMax), FCOS1x(TMax), FCOS1y(TMax), FCOS1z(TMax))
		ALLOCATE(FCOS2(TMax), FCOS2x(TMax), FCOS2y(TMax), FCOS2z(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			FCOS1(I) = 0.0
			FCOS1x(I) = 0.0
			FCOS1y(I) = 0.0
			FCOS1z(I) = 0.0
			FCOS2(I) = 0.0
			FCOS2x(I) = 0.0
			FCOS2y(I) = 0.0
			FCOS2z(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate FSQ_', A1)") FType
				
		ALLOCATE(Qx(QIMax), Qy(QIMax), Qz(QIMax))

!		====set wave vectors matching the given wave number for Fs(q,t)
		CALL Set_FSQ
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=85, FILE='FSQ_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyna(I)>0) THEN
!				====Fs(q, t)
				FSQ = FCOS1(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				FSQx = FCOS1x(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				FSQy = FCOS1y(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				FSQz = FCOS1z(I)/(DBLE(NPart)*DBLE(NDyna(I)))
!				====Xs(q, t)
				FSQ1 = FCOS1(I)/DBLE(NDyna(I))
				FSQ1x = FCOS1x(I)/DBLE(NDyna(I))
				FSQ1y = FCOS1y(I)/DBLE(NDyna(I))
				FSQ1z = FCOS1z(I)/DBLE(NDyna(I))
				FSQ2 = FCOS2(I)/DBLE(NDyna(I))
				FSQ2x = FCOS2x(I)/DBLE(NDyna(I))
				FSQ2y = FCOS2y(I)/DBLE(NDyna(I))
				FSQ2z = FCOS2z(I)/DBLE(NDyna(I))
				XSQ = (FSQ2 - FSQ1**2)/DBLE(NPart)
				XSQx = (FSQ2x - FSQ1x**2)/DBLE(NPart)
				XSQy = (FSQ2y - FSQ1y**2)/DBLE(NPart)
				XSQz = (FSQ2z - FSQ1z**2)/DBLE(NPart)
				WRITE (85, "(20E17.8)") MDTime, FSQ, XSQ, FSQx, XSQx, FSQy, XSQy, FSQz, XSQz
			END IF
		END DO
		CLOSE (85)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(FCOS1, FCOS1x, FCOS1y, FCOS1z)
		DEALLOCATE(FCOS2, FCOS2x, FCOS2y, FCOS2z)
		DEALLOCATE(Qx, Qy, Qz)
	END IF

RETURN
END

!====SUBROUTINE: Set_FSQ--set up wave vectors for Fs(q,t)
SUBROUTINE Set_FSQ
	
	USE SetCOM
	USE SetFSQ
	IMPLICIT NONE 
	
	INTEGER QiSet, Kx, Ky, Kz, K2, Q2

	DOUBLE PRECISION QrSet, DelQ, DelMin
!===================================================	
!	====read the target wavenumber
	OPEN (UNIT=21, FILE='Qset.txt')
	READ (21, *) QrSet
	CLOSE (21)
	
	Qone = 2*PI/Box*INT(ANINT(QrSet*Box/(2*PI)))
	
	QiSet = INT(QrSet*Box/(2*PI)) + 1
	
!	====find wave number matching the real wave number
	DelMin = 1.0D10
	DO Kx = 0, QiSet
		DO Ky = -QiSet, QiSet
			DO Kz = -QiSet, QiSet
				K2 = Kx**2 + Ky**2 + Kz**2
				DelQ = ABS(2*PI/Box*SQRT(DBLE(K2)) - QrSet)
				IF (DelQ<DelMin) THEN
					DelMin = DelQ
					Q2 = K2
				END IF
			END DO
		END DO
	END DO
	
!	====set wave vectors matching the wave number
	Qn = 0
	DO Kx = 0, QiSet
		DO Ky = -QiSet, QiSet
			DO Kz = -QiSet, QiSet
				K2 = Kx**2 + Ky**2 + Kz**2
				IF (K2==Q2 .AND. Qn<QIMax) THEN
					Qn = Qn + 1
					Qx(Qn) = 2*PI/Box*Kx
					Qy(Qn) = 2*PI/Box*Ky
					Qz(Qn) = 2*PI/Box*Kz
				END IF
			END DO
		END DO
	END DO
	
	WRITE (12, "('    Chosen Q         :', F10.4)") Qone
	WRITE (12, "('    Chosen Q         :', F10.4)") 2*PI/Box*SQRT(DBLE(Q2))
	WRITE (12, "('    Number of vectors:', I10)")  Qn

RETURN
END
