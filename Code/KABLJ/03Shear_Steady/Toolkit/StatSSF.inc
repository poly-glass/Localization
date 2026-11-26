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
	
	INTEGER Switch, NStat, NDist, I, J, K, U, V
	
	DOUBLE PRECISION Kx, Ky, Kz, SumCOS, SumSIN, SSF
	DOUBLE PRECISION, ALLOCATABLE:: dx(:), dy(:), dz(:)
	DOUBLE PRECISION, ALLOCATABLE:: Q(:), Sbb(:)
	
	SAVE NStat, NDist, Q, Sbb, dx, dy, dz
!==========================================================
	IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
		I = 0
		DO U = 1, NPart
			DO V = U + 1, NPart
				I = I + 1
				dx(I) = MX(U) - MX(V)
				dy(I) = MY(U) - MY(V)
				dz(I) = MZ(U) - MZ(V)
				dx(I) = dx(I) - DelX*ANINT(dy(I)/Box)
				dx(I) = dx(I) - Box*ANINT(dx(I)/Box)
				dy(I) = dy(I) - Box*ANINT(dy(I)/Box)
				dz(I) = dz(I) - Box*ANINT(dz(I)/Box)
			END DO
		END DO
		
		DO K = 1, NQ
			DO J = 1, Qn(K)
				Kx = Qx(K, J)
				Ky = Qy(K, J)
				Kz = Qz(K, J)
				
				SumCOS = 0.0
				DO I = 1, NDist
					SumCOS = SumCOS + COS(Kx*dx(I) + Ky*dy(I) + Kz*dz(I))
				END DO
				Sbb(K) = Sbb(K) + 2*SumCOS
			END DO
		END DO
	ELSE IF (Switch==0) THEN 
!		====initialization
		NStat = 0
		NDist = NPart*(NPart - 1)/2
		
!		====allocate arrays
		ALLOCATE(dx(NDist), dy(NDist), dz(NDist))
		ALLOCATE(Qn(QMax), Qc(QMax))
		ALLOCATE(Qx(QMax, QIMax), Qy(QMax, QIMax), Qz(QMax, QIMax))

!		====set wavenumber for S(q)
		CALL Set_SSF
		
!		====allocate arrays
		ALLOCATE(Q(NQ), Sbb(NQ))
		
		DO I = 1, NQ
			Q(I) = 2*PI/Box*DBLE(Qc(I))**0.5
			Sbb(I) = 0.0
		END DO

		WRITE (12, "('  To calculate SSF')")
		WRITE (12, "('    Number of Q        : ', I10)") NQ
		WRITE (12, "('    Minimum & maximum Q: ', 2F10.4)") Q(1), Q(NQ)
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=65, FILE='SSF')
		DO I = 1, NQ
			SSF = Sbb(I)/(DBLE(NStat)*DBLE(Qn(I))*DBLE(NPart)) + 1
			WRITE(65, "(10E17.8)") Q(I), SSF
		END DO
		CLOSE (65)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(dx, dy, dz)
		DEALLOCATE(Qn, Qc)
		DEALLOCATE(Qx, Qy, Qz)
		DEALLOCATE(Q, Sbb)
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
