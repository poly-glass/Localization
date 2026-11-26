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
	DOUBLE PRECISION SetCOS, SumCOS0, FSQ0, XSQ0, FSQ01, FSQ02
	DOUBLE PRECISION SetCOSx, SumCOS0x, FSQ0x, XSQ0x, FSQ01x, FSQ02x
	DOUBLE PRECISION SetCOSy, SumCOS0y, FSQ0y, XSQ0y, FSQ01y, FSQ02y
	DOUBLE PRECISION SetCOSz, SumCOS0z, FSQ0z, XSQ0z, FSQ01z, FSQ02z
	DOUBLE PRECISION SumCOS1, FSQ1, XSQ1, FSQ11, FSQ12
	DOUBLE PRECISION SumCOS1x, FSQ1x, XSQ1x, FSQ11x, FSQ12x
	DOUBLE PRECISION SumCOS1y, FSQ1y, XSQ1y, FSQ11y, FSQ12y
	DOUBLE PRECISION SumCOS1z, FSQ1z, XSQ1z, FSQ11z, FSQ12z
	DOUBLE PRECISION SumCOS2, FSQ2, XSQ2, FSQ21, FSQ22
	DOUBLE PRECISION SumCOS2x, FSQ2x, XSQ2x, FSQ21x, FSQ22x
	DOUBLE PRECISION SumCOS2y, FSQ2y, XSQ2y, FSQ21y, FSQ22y
	DOUBLE PRECISION SumCOS2z, FSQ2z, XSQ2z, FSQ21z, FSQ22z

	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS01(:), FCOS01x(:), FCOS01y(:), FCOS01z(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS02(:), FCOS02x(:), FCOS02y(:), FCOS02z(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS11(:), FCOS11x(:), FCOS11y(:), FCOS11z(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS12(:), FCOS12x(:), FCOS12y(:), FCOS12z(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS21(:), FCOS21x(:), FCOS21y(:), FCOS21z(:)
	DOUBLE PRECISION, ALLOCATABLE:: FCOS22(:), FCOS22x(:), FCOS22y(:), FCOS22z(:)
	
	SAVE Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0
	SAVE FCOS01, FCOS01x, FCOS01y, FCOS01z, FCOS02, FCOS02x, FCOS02y, FCOS02z
	SAVE FCOS11, FCOS11x, FCOS11y, FCOS11z, FCOS12, FCOS12x, FCOS12y, FCOS12z
	SAVE FCOS21, FCOS21x, FCOS21y, FCOS21z, FCOS22, FCOS22x, FCOS22y, FCOS22z
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
				
				SumCOS0 = 0.0
				SumCOS1 = 0.0
				SumCOS2 = 0.0
				SumCOS0x = 0.0
				SumCOS0y = 0.0
				SumCOS0z = 0.0
				SumCOS1x = 0.0
				SumCOS1y = 0.0
				SumCOS1z = 0.0
				SumCOS2x = 0.0
				SumCOS2y = 0.0
				SumCOS2z = 0.0

				DO I = 1, NPart
					dx = MX(I) - RX0(I, t)
					dy = MY(I) - RY0(I, t)
					dz = MZ(I) - RZ0(I, t)
					dx = dx - Rate*(dt-1)*Delt*RY0(I, t)
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
					SetCOSx = COS(Qone*dx)
					SetCOSy = COS(Qone*dy)
					SetCOSz = COS(Qone*dz)
					SumCOS0x = SumCOS0x + SetCOSx
					SumCOS0y = SumCOS0y + SetCOSy
					SumCOS0z = SumCOS0z + SetCOSz
					IF(PType(I)==1) THEN
						SumCOS1x = SumCOS1x + SetCOSx
						SumCOS1y = SumCOS1y + SetCOSy
						SumCOS1z = SumCOS1z + SetCOSz
					END IF
					IF(PType(I)==2) THEN
						SumCOS2x = SumCOS2x + SetCOSx
						SumCOS2y = SumCOS2y + SetCOSy
						SumCOS2z = SumCOS2z + SetCOSz
					END IF
				END DO
				FCOS01(dt) = FCOS01(dt) + (SumCOS0/Qn)
				FCOS01x(dt) = FCOS01x(dt) + (SumCOS0x)
				FCOS01y(dt) = FCOS01y(dt) + (SumCOS0y)
				FCOS01z(dt) = FCOS01z(dt) + (SumCOS0z)
				FCOS02(dt) = FCOS02(dt) + (SumCOS0/Qn)**2
				FCOS02x(dt) = FCOS02x(dt) + (SumCOS0x)**2
				FCOS02y(dt) = FCOS02y(dt) + (SumCOS0y)**2
				FCOS02z(dt) = FCOS02z(dt) + (SumCOS0z)**2
				!====Fs(q, t) for particle A
				FCOS11(dt) = FCOS11(dt) + (SumCOS1/Qn)
				FCOS11x(dt) = FCOS11x(dt) + (SumCOS1x)
				FCOS11y(dt) = FCOS11y(dt) + (SumCOS1y)
				FCOS11z(dt) = FCOS11z(dt) + (SumCOS1z)
				FCOS12(dt) = FCOS12(dt) + (SumCOS1/Qn)**2
				FCOS12x(dt) = FCOS12x(dt) + (SumCOS1x)**2
				FCOS12y(dt) = FCOS12y(dt) + (SumCOS1y)**2
				FCOS12z(dt) = FCOS12z(dt) + (SumCOS1z)**2
				!====Fs(q, t) for particle B
				FCOS21(dt) = FCOS21(dt) + (SumCOS2/Qn)
				FCOS21x(dt) = FCOS21x(dt) + (SumCOS2x)
				FCOS21y(dt) = FCOS21y(dt) + (SumCOS2y)
				FCOS21z(dt) = FCOS21z(dt) + (SumCOS2z)
				FCOS22(dt) = FCOS22(dt) + (SumCOS2/Qn)**2
				FCOS22x(dt) = FCOS22x(dt) + (SumCOS2x)**2
				FCOS22y(dt) = FCOS22y(dt) + (SumCOS2y)**2
				FCOS22z(dt) = FCOS22z(dt) + (SumCOS2z)**2
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		Nt0 = 0

!		===allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(FCOS01(TMax), FCOS01x(TMax), FCOS01y(TMax), FCOS01z(TMax))
		ALLOCATE(FCOS02(TMax), FCOS02x(TMax), FCOS02y(TMax), FCOS02z(TMax))
		ALLOCATE(FCOS11(TMax), FCOS11x(TMax), FCOS11y(TMax), FCOS11z(TMax))
		ALLOCATE(FCOS12(TMax), FCOS12x(TMax), FCOS12y(TMax), FCOS12z(TMax))
		ALLOCATE(FCOS21(TMax), FCOS21x(TMax), FCOS21y(TMax), FCOS21z(TMax))
		ALLOCATE(FCOS22(TMax), FCOS22x(TMax), FCOS22y(TMax), FCOS22z(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			FCOS01(I) = 0.0
			FCOS01x(I) = 0.0
			FCOS01y(I) = 0.0
			FCOS01z(I) = 0.0
			FCOS02(I) = 0.0
			FCOS02x(I) = 0.0
			FCOS02y(I) = 0.0
			FCOS02z(I) = 0.0
			FCOS11(I) = 0.0
			FCOS11x(I) = 0.0
			FCOS11y(I) = 0.0
			FCOS11z(I) = 0.0
			FCOS12(I) = 0.0
			FCOS12x(I) = 0.0
			FCOS12y(I) = 0.0
			FCOS12z(I) = 0.0
			FCOS21(I) = 0.0
			FCOS21x(I) = 0.0
			FCOS21y(I) = 0.0
			FCOS21z(I) = 0.0
			FCOS22(I) = 0.0
			FCOS22x(I) = 0.0
			FCOS22y(I) = 0.0
			FCOS22z(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate FSQ_', A1)") FType
				
		ALLOCATE(Qx(QIMax), Qy(QIMax), Qz(QIMax))

!		====set wave vectors matching the given wave number for Fs(q,t)
		CALL Set_FSQ
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=85, FILE='FSQ0_'//FType)
		OPEN (UNIT=86, FILE='FSQ1_'//FType)
		OPEN (UNIT=87, FILE='FSQ2_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyna(I)>0) THEN
!				====Fs(q, t)
				FSQ0 = FCOS01(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				FSQ0x = FCOS01x(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				FSQ0y = FCOS01y(I)/(DBLE(NPart)*DBLE(NDyna(I)))
				FSQ0z = FCOS01z(I)/(DBLE(NPart)*DBLE(NDyna(I)))
!				====Xs(q, t)
				FSQ01 = FCOS01(I)/DBLE(NDyna(I))
				FSQ01x = FCOS01x(I)/DBLE(NDyna(I))
				FSQ01y = FCOS01y(I)/DBLE(NDyna(I))
				FSQ01z = FCOS01z(I)/DBLE(NDyna(I))
				FSQ02 = FCOS02(I)/DBLE(NDyna(I))
				FSQ02x = FCOS02x(I)/DBLE(NDyna(I))
				FSQ02y = FCOS02y(I)/DBLE(NDyna(I))
				FSQ02z = FCOS02z(I)/DBLE(NDyna(I))
				XSQ0 = (FSQ02 - FSQ01**2)/DBLE(NPart)
				XSQ0x = (FSQ02x - FSQ01x**2)/DBLE(NPart)
				XSQ0y = (FSQ02y - FSQ01y**2)/DBLE(NPart)
				XSQ0z = (FSQ02z - FSQ01z**2)/DBLE(NPart)
				WRITE (85, "(20E17.8)") MDTime, FSQ0, XSQ0, FSQ0x, XSQ0x, FSQ0y, XSQ0y, FSQ0z, XSQ0z
!				====Fs(q, t) for particle A
				FSQ1 = FCOS11(I)/(DBLE(NPart1)*DBLE(NDyna(I)))
				FSQ1x = FCOS11x(I)/(DBLE(NPart1)*DBLE(NDyna(I)))
				FSQ1y = FCOS11y(I)/(DBLE(NPart1)*DBLE(NDyna(I)))
				FSQ1z = FCOS11z(I)/(DBLE(NPart1)*DBLE(NDyna(I)))
!				====Xs(q, t) for particle A
				FSQ11 = FCOS11(I)/DBLE(NDyna(I))
				FSQ11x = FCOS11x(I)/DBLE(NDyna(I))
				FSQ11y = FCOS11y(I)/DBLE(NDyna(I))
				FSQ11z = FCOS11z(I)/DBLE(NDyna(I))
				FSQ12 = FCOS12(I)/DBLE(NDyna(I))
				FSQ12x = FCOS12x(I)/DBLE(NDyna(I))
				FSQ12y = FCOS12y(I)/DBLE(NDyna(I))
				FSQ12z = FCOS12z(I)/DBLE(NDyna(I))
				XSQ1 = (FSQ12 - FSQ11**2)/DBLE(NPart1)
				XSQ1x = (FSQ12x - FSQ11x**2)/DBLE(NPart1)
				XSQ1y = (FSQ12y - FSQ11y**2)/DBLE(NPart1)
				XSQ1z = (FSQ12z - FSQ11z**2)/DBLE(NPart1)
				WRITE (86, "(20E17.8)") MDTime, FSQ1, XSQ1, FSQ1x, XSQ1x, FSQ1y, XSQ1y, FSQ1z, XSQ1z
!				====Fs(q, t) for particle B
				FSQ2 = FCOS21(I)/(DBLE(NPart2)*DBLE(NDyna(I)))
				FSQ2x = FCOS21x(I)/(DBLE(NPart2)*DBLE(NDyna(I)))
				FSQ2y = FCOS21y(I)/(DBLE(NPart2)*DBLE(NDyna(I)))
				FSQ2z = FCOS21z(I)/(DBLE(NPart2)*DBLE(NDyna(I)))
!				====Xs(q, t) for particle B
				FSQ21 = FCOS21(I)/DBLE(NDyna(I))
				FSQ21x = FCOS21x(I)/DBLE(NDyna(I))
				FSQ21y = FCOS21y(I)/DBLE(NDyna(I))
				FSQ21z = FCOS21z(I)/DBLE(NDyna(I))
				FSQ22 = FCOS22(I)/DBLE(NDyna(I))
				FSQ22x = FCOS22x(I)/DBLE(NDyna(I))
				FSQ22y = FCOS22y(I)/DBLE(NDyna(I))
				FSQ22z = FCOS22z(I)/DBLE(NDyna(I))
				XSQ2 = (FSQ22 - FSQ21**2)/DBLE(NPart2)
				XSQ2x = (FSQ22x - FSQ21x**2)/DBLE(NPart2)
				XSQ2y = (FSQ22y - FSQ21y**2)/DBLE(NPart2)
				XSQ2z = (FSQ22z - FSQ21z**2)/DBLE(NPart2)
				WRITE (87, "(20E17.8)") MDTime, FSQ2, XSQ2, FSQ2x, XSQ2x, FSQ2y, XSQ2y, FSQ2z, XSQ2z
			END IF
		END DO
		CLOSE (85)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(FCOS01, FCOS01x, FCOS01y, FCOS01z)
		DEALLOCATE(FCOS02, FCOS02x, FCOS02y, FCOS02z)
		DEALLOCATE(FCOS11, FCOS11x, FCOS11y, FCOS11z)
		DEALLOCATE(FCOS12, FCOS12x, FCOS12y, FCOS12z)
		DEALLOCATE(FCOS21, FCOS21x, FCOS21y, FCOS21z)
		DEALLOCATE(FCOS22, FCOS22x, FCOS22y, FCOS22z)
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
