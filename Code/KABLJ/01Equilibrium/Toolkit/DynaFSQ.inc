!====SetFSQ: global variables and arrays for FSQ
MODULE SetFSQ
	
	INTEGER, PARAMETER:: QIMax=10
	
	INTEGER Qn
	
	DOUBLE PRECISION, ALLOCATABLE:: Qx(:), Qy(:), Qz(:)
	
END MODULE

!====Get_FSQ: compute Self-Intermediate Scattering Function Fs(q,t)
SUBROUTINE Get_FSQ(Switch)
	
	USE SetCOM
	USE SetFSQ
	IMPLICIT NONE
	
	INTEGER Switch, Nt, Nt0, ttel, t, dt, I, K, Natom(0:2)
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	
	DOUBLE PRECISION MDTime, dx, dy, dz, SetCOS
	DOUBLE PRECISION SumCOS0, FSQ0, XSQ0, FSQ01, FSQ02
	DOUBLE PRECISION SumCOS1, FSQ1, XSQ1, FSQ11, FSQ12
	DOUBLE PRECISION SumCOS2, FSQ2, XSQ2, FSQ21, FSQ22
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS01(:), FCOS02(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS11(:), FCOS12(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS21(:), FCOS22(:)
	
	SAVE Natom, Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0
	SAVE FCOS01, FCOS02, FCOS11, FCOS12, FCOS21, FCOS22
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
				
				SumCOS0 = 0.0
				SumCOS1 = 0.0
				SumCOS2 = 0.0
				DO I = 1, NPart
					dx = RX(I) - RX0(I, t)
					dy = RY(I) - RY0(I, t)
					dz = RZ(I) - RZ0(I, t)
					DO K = 1, Qn
						SetCOS = COS(Qx(K)*dx + Qy(K)*dy + Qz(K)*dz)
						SumCOS0 = SumCOS0 + SetCOS
						IF(PType(I)==1) THEN
							SumCOS1 = SumCOS1 + SetCOS
						END IF
						IF(PType(I)==2) THEN
							SumCOS2 = SumCOS2 + SetCOS
						END IF
					END DO
				END DO
				FCOS01(dt) = FCOS01(dt) + (SumCOS0/Qn)
				FCOS02(dt) = FCOS02(dt) + (SumCOS0/Qn)**2
				FCOS11(dt) = FCOS11(dt) + (SumCOS1/Qn)
				FCOS12(dt) = FCOS12(dt) + (SumCOS1/Qn)**2
				FCOS21(dt) = FCOS21(dt) + (SumCOS2/Qn)
				FCOS22(dt) = FCOS22(dt) + (SumCOS2/Qn)**2
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
		ALLOCATE(FCOS01(TMax), FCOS02(TMax))
		ALLOCATE(FCOS11(TMax), FCOS12(TMax))
		ALLOCATE(FCOS21(TMax), FCOS22(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			FCOS01(I) = 0.0
			FCOS02(I) = 0.0
			FCOS11(I) = 0.0
			FCOS12(I) = 0.0
			FCOS21(I) = 0.0
			FCOS22(I) = 0.0
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

				FSQ0 = FCOS01(I)/(DBLE(Natom(0))*DBLE(NDyna(I)))
				FSQ1 = FCOS11(I)/(DBLE(Natom(1))*DBLE(NDyna(I)))
				
				IF(Natom(2) == 0)THEN
					FSQ2 = 0.0
				ELSE
					FSQ2 = FCOS21(I)/(DBLE(Natom(2))*DBLE(NDyna(I)))
				END IF
										
!				====Xs(q, t)
				FSQ01 = FCOS01(I)/DBLE(NDyna(I))
				FSQ02 = FCOS02(I)/DBLE(NDyna(I))
				FSQ11 = FCOS11(I)/DBLE(NDyna(I))
				FSQ12 = FCOS12(I)/DBLE(NDyna(I))
				FSQ21 = FCOS21(I)/DBLE(NDyna(I))
				FSQ22 = FCOS22(I)/DBLE(NDyna(I))
				XSQ0 = (FSQ02 - FSQ01**2)/DBLE(Natom(0))
				XSQ1 = (FSQ12 - FSQ11**2)/DBLE(Natom(1))
				
				IF(Natom(2) == 0)THEN
					XSQ2 = 0.0
				ELSE
					XSQ2 = (FSQ22 - FSQ21**2)/DBLE(Natom(2))
				END IF
					
				WRITE (85, "(20E17.8)") MDTime, FSQ0, FSQ1, FSQ2, XSQ0, XSQ1, XSQ2
			END IF
		END DO
		CLOSE (85)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(FCOS01, FCOS02, FCOS11, FCOS12, FCOS21, FCOS22)
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
