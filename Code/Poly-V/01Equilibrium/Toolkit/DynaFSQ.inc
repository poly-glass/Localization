!====SetFSQ: global variables and arrays for FSQ
MODULE SetFSQ
	
	INTEGER, PARAMETER:: QIMax=20
	
	INTEGER Qn
	
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
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS1(:), FCOS2(:)
	
	SAVE Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0, FCOS1, FCOS2
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
				
				SumCOS = 0.0
				DO I = 1, NPart
					dx = RX(I) - RX0(I, t)
					dy = RY(I) - RY0(I, t)
					dz = RZ(I) - RZ0(I, t)
					DO K = 1, Qn
						SetCOS = COS(Qx(K)*dx + Qy(K)*dy + Qz(K)*dz)
						SumCOS = SumCOS + SetCOS
					END DO
				END DO
				FCOS1(dt) = FCOS1(dt) + (SumCOS/Qn)
				FCOS2(dt) = FCOS2(dt) + (SumCOS/Qn)**2
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0

!		===allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(FCOS1(TMax), FCOS2(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			FCOS1(I) = 0.0
			FCOS2(I) = 0.0
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
!				====Xs(q, t)
				FSQ1 = FCOS1(I)/DBLE(NDyna(I))
				FSQ2 = FCOS2(I)/DBLE(NDyna(I))
				XSQ = (FSQ2 - FSQ1**2)/DBLE(NPart)
				WRITE (85, "(10E17.8)") MDTime, FSQ, XSQ
			END IF
		END DO
		CLOSE (85)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(FCOS1, FCOS2)
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

	WRITE (12, "('    Chosen Q         :', F10.4)") 2*PI/Box*SQRT(DBLE(Q2))
	WRITE (12, "('    Number of vectors:', I10)")  Qn

RETURN
END
