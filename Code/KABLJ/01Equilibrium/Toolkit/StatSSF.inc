!====SetSSF: global variables and arrays for SSF
MODULE SetSSF
	
	INTEGER, PARAMETER:: QMax=10000, QIMax=20
	
	INTEGER NQ
	INTEGER, ALLOCATABLE:: Qn(:), Qc(:)

	DOUBLE PRECISION, PARAMETER:: Qrh=15.0
	
	DOUBLE PRECISION, ALLOCATABLE:: Qx(:,:), Qy(:,:), Qz(:,:)
	
END MODULE

!=====Get_SSF: compute Static Structure Factor
SUBROUTINE Get_SSF(Switch)
	
	USE SetCOM
	USE SetSSF
	IMPLICIT NONE
	
	INTEGER Switch, NStat, I, J, K, Natom(0:2)
	
	DOUBLE PRECISION Kx, Ky, Kz, SetCOS, SetSIN, SSF00, SSF11, SSF22, SSF12
	DOUBLE PRECISION SumCOS0, SumSIN0, SumCOS1, SumSIN1, SumCOS2, SumSIN2
	DOUBLE PRECISION, ALLOCATABLE:: Q(:), S00(:), S11(:), S22(:), S12(:)
	
	SAVE Natom, NStat, Q, S00, S11, S22, S12
!==========================================================
IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1		
!		====whole system
		DO K = 1, NQ
			DO J = 1, Qn(K)
				Kx = Qx(K, J)
				Ky = Qy(K, J)
				Kz = Qz(K, J)

				SumCOS1 = 0.0
				SumSIN1 = 0.0
				SumCOS0 = 0.0
				SumSIN0 = 0.0
				SumCOS2 = 0.0
				SumSIN2 = 0.0
				
				DO I = 1, NPart
					SetCOS = COS(Kx*RX(I) + Ky*RY(I) + Kz*RZ(I))
					SetSIN = SIN(Kx*RX(I) + Ky*RY(I) + Kz*RZ(I))
					SumCOS0 = SumCOS0 + SetCOS
 					SumSIN0 = SumSIN0 + SetSIN

					IF (PType(I)==1) THEN
						SumCOS1 = SumCOS1 + SetCOS
						SumSIN1 = SumSIN1 + SetSIN
						
					END IF
					IF (PType(I)==2) THEN
						SumCOS2 = SumCOS2 + SetCOS
						SumSIN2 = SumSIN2 + SetSIN

					END IF					
				END DO
				S00(K) = S00(K) + SumCOS0*SumCOS0 + SumSIN0*SumSIN0
				S12(K) = S12(K) + SumCOS1*SumCOS2 + SumSIN1*SumSIN2
				S22(K) = S22(K) + SumCOS2*SumCOS2 + SumSIN2*SumSIN2
				S11(K) = S11(K) + SumCOS1*SumCOS1 + SumSIN1*SumSIN1
			END DO
		END DO


	ELSE IF (Switch==0) THEN 
!		====initialization
		Natom(0) = NPart
		Natom(1) = NPart1
		Natom(2) = NPart2
		
		NStat = 0
		
!		====allocate arrays
		ALLOCATE(Qn(QMax), Qc(QMax))
		ALLOCATE(Qx(QMax, QIMax), Qy(QMax, QIMax), Qz(QMax, QIMax))

!		====set wavenumber for S(q)
		CALL Set_SSF
		
!		====allocate arrays
		ALLOCATE(Q(NQ), S00(NQ), S11(NQ), S22(NQ), S12(NQ))
		
		DO I = 1, NQ
			Q(I) = 2*PI/Box*DBLE(Qc(I))**0.5
			S00(I) = 0.0
			S11(I) = 0.0
			S22(I) = 0.0
			S12(I) = 0.0
		END DO

		WRITE (12, "('  To calculate SSF')")
		WRITE (12, "('    Number of Q        : ', I10)") NQ
		WRITE (12, "('    Minimum & maximum Q: ', 2F10.4)") Q(1), Q(NQ)
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=65, FILE='SSF')
		DO I = 1, NQ
			SSF00 = S00(I)/(DBLE(NStat)*DBLE(Qn(I))*(DBLE(Natom(0))*DBLE(Natom(0)))**0.5)
			SSF11 = S11(I)/(DBLE(NStat)*DBLE(Qn(I))*(DBLE(Natom(1))*DBLE(Natom(1)))**0.5)
			
			IF(Natom(2) == 0)THEN
				SSF22 = 0.0
				SSF12 = 0.0		
			ELSE
				SSF22 = S22(I)/(DBLE(NStat)*DBLE(Qn(I))*(DBLE(Natom(2))*DBLE(Natom(2)))**0.5)
				SSF12 = S12(I)/(DBLE(NStat)*DBLE(Qn(I))*(DBLE(Natom(1))*DBLE(Natom(2)))**0.5)
			END IF
			
			WRITE(65, "(10E17.8)") Q(I), SSF00, SSF11, SSF22, SSF12
		END DO
		CLOSE (65)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(Qn, Qc)
		DEALLOCATE(Qx, Qy, Qz)
		DEALLOCATE(Q, S00, S11, S22, S12)
	END IF

RETURN
END 

!====SUBROUTINE: Set_SSF--set up wave vectors for static structure factor
SUBROUTINE Set_SSF
	
	USE SetCOM
	USE SetSSF
	IMPLICIT NONE
	
	INTEGER I, J, K, Qih, Nk, Kx, Ky, Kz, K2, Ks
	INTEGER, ALLOCATABLE:: Kn(:), Kc(:)
!===================================================================
	Qih = INT(Qrh*Box/(2*PI)) + 1
	
!	====allocate arrays
	ALLOCATE(Kn(QMax), Kc(QMax))
	
!	====set the square of wave numbers up to Qih^2
	Ks = 8
	K = 0
	DO I = Ks, Ks*Qih
		K = K + 1
		Kc(K) = INT(ANINT((DBLE(I)/DBLE(Ks))**2))
	END DO
	
!	====remove those the same as previous ones
	Nk = 1 
	DO I = 2, K
		IF (Kc(I)-Kc(I-1)>0) THEN
			Nk = Nk + 1
			Kc(Nk) = Kc(I)
		END IF
	END DO

!	====find wave vectors matching each wave number
	DO I = 1, Nk
		Kn(I) = 0
		DO Kx = 0, Qih
			DO Ky = -Qih, Qih
				DO Kz = -Qih, Qih
					K2 = Kx**2 + Ky**2 + Kz**2
					IF (K2==Kc(I) .AND. Kn(I)<QIMax) THEN
						Kn(I) = Kn(I) + 1
						Qx(I, Kn(I)) = Kx
						Qy(I, Kn(I)) = Ky
						Qz(I, Kn(I)) = Kz
					END IF
				END DO
			END DO
		END DO
	END DO 

!	====set number of wave numbers and wave vectors for S(q)
	NQ = 0
	DO I = 1, Nk
		IF (Kn(I)>0) THEN
			NQ = NQ + 1
			Qn(NQ) = Kn(I)
			Qc(NQ) = Kc(I)
			DO J = 1, Kn(I)
				Qx(NQ, J) = 2*PI/Box*Qx(I, J)
				Qy(NQ, J) = 2*PI/Box*Qy(I, J)
				Qz(NQ, J) = 2*PI/Box*Qz(I, J)
			END DO
		END IF
	END DO
	
!	====free memory
	DEALLOCATE(Kn, Kc)

RETURN
END
