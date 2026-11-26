!====SetPCS: global variables and arrays for PCS
MODULE SetPCS
	
	INTEGER Nmb
	INTEGER, ALLOCATABLE:: Pmb(:)
	
	DOUBLE PRECISION Rcut, Rcut2
	DOUBLE PRECISION, ALLOCATABLE:: r2(:)
	
END MODULE

!====Get_PCS: compute average string length
SUBROUTINE Get_PCS(Switch)
	
	USE SetCOM
	USE SetPCS
	IMPLICIT NONE
	
	INTEGER Switch, Nt, Nt0, ttel, t, dt, I, J, K, IT, LK, LIT, NCluster, N1, N2
	INTEGER, ALLOCATABLE:: NDyna(:), ttv0(:)
	INTEGER, ALLOCATABLE:: NIT(:), L(:), NumN(:)

	DOUBLE PRECISION MDTime, dx, dy, dz, dr2, PCS, FN1, FN2
	DOUBLE PRECISION, ALLOCATABLE:: RX0(:,:), RY0(:,:), RZ0(:,:)
	DOUBLE PRECISION, ALLOCATABLE:: ST(:),F1(:), F2(:)
	
	SAVE Nt, Nt0, NDyna, ttv0, RX0, RY0, RZ0, ST, NIT, L,NumN, F1, F2
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

				DO I = 1, NPart
					dx = RX(I) - RX0(I, t)
					dy = RY(I) - RY0(I, t)
					dz = RZ(I) - RZ0(I, t)
					r2(I) = dx**2 + dy**2 + dz**2
				END DO

!				====string-like collective motion
!				====get mobile particles
				CALL Get_Pmb1
				
!				====sort the clusters
				DO I = 1, Nmb
					L(I) = I
				END DO
				
				DO I = 1, Nmb - 1
					IF (I==L(I)) THEN
						J = I
						DO K = I + 1, Nmb
							LK = L(K)
							IF (LK==K) THEN
								dx = RX(Pmb(J)) - RX(Pmb(K))
								dy = RY(Pmb(J)) - RY(Pmb(K))
								dz = RZ(Pmb(J)) - RZ(Pmb(K))
								dx = dx - BOX*ANINT(dx/BOX)
								dy = dy - BOX*ANINT(dy/BOX)
								dz = dz - BOX*ANINT(dz/BOX)
								dr2 = dx**2 + dy**2 + dz**2
								IF (dr2<Rcut2) THEN
									L(K) = L(J)
									L(J) = LK
								END IF
							END IF
						END DO
						J = L(J)
10						IF (J .NE. I) THEN
							DO K = I + 1, Nmb
								LK = L(K)
								IF (LK==K) THEN
									dx = RX(Pmb(J)) - RX(Pmb(K))
									dy = RY(Pmb(J)) - RY(Pmb(K))
									dz = RZ(Pmb(J)) - RZ(Pmb(K))
									dx = dx - BOX*ANINT(dx/BOX)
									dy = dy - BOX*ANINT(dy/BOX)
									dz = dz - BOX*ANINT(dz/BOX)
									dr2 = dx**2 + dy**2 + dz**2
									IF (dr2<Rcut2) THEN
										L(K) = L(J)
										L(J) = LK
									END IF
								END IF
							END DO
							J = L(J)
							GOTO 10
						END IF
					END IF
				END DO
								
!				====count number in a cluster containing particle IT
				DO IT = 1, Nmb
					NIT(IT) = 1
					LIT = L(IT)
30					IF (LIT .NE. IT) THEN
						NIT(IT) = NIT(IT) + 1
						LIT = L(LIT)
						GOTO 30
					END IF
				END DO
		
!				====count number of total particles for a given cluster size
				DO I = 1, Nmb
					NumN(I) = 0
				END DO
				
				DO IT = 1, Nmb
					DO I = 1, Nmb
						IF (NIT(IT)==I) THEN
							NumN(I) = NumN(I) + 1
						END IF
					END DO
				END DO
				
				NCluster = 0
				DO I = 1, Nmb
					NCluster = NCluster + NumN(I)/I
				END DO
				ST(dt) = ST(dt) + DBLE(Nmb)/DBLE(NCluster)			
!				====statistics of mobile particles
				N1 = 0
				N2 = 0
				DO I = 1, Nmb
					IF (PType(Pmb(I))==1) N1 = N1 + 1
					IF (PType(Pmb(I))==2) N2 = N2 + 1
				END DO
				F1(dt) = F1(dt) + DBLE(N1)/DBLE(Nmb)
				F2(dt) = F2(dt) + DBLE(N2)/DBLE(Nmb)
			END IF
		END DO
	ELSE IF (Switch==0) THEN
!		====initialization
		Nt = 0
		
!		===allocate arrays
		ALLOCATE(NDyna(TMax), ttv0(T0Max))
		ALLOCATE(RX0(NPart, T0Max), RY0(NPart, T0Max), RZ0(NPart, T0Max))
		ALLOCATE(ST(TMax), F1(TMax), F2(TMax))
		
		DO I = 1, TMax
			NDyna(I) = 0
			ST(I) = 0.0
			F1(I) = 0.0
			F2(I) = 0.0
		END DO
		
		!		====read cutcoff
		Rcut = 1.4
		Rcut2 = Rcut**2

!		====number of mobile particles and threshold for mobile cluster
		Nmb = INT(0.05*NPart)

!		===allocate arrays
		ALLOCATE(Pmb(NPart), r2(NPart), NIT(Nmb), L(Nmb), NumN(Nmb))

		WRITE (12, "('  To calculate PCS_', A1)") FType
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=87, FILE='PCS_'//FType)
		DO I = 1, TMax
			MDTime = Delt*(I - 1)
			IF (NDyna(I)>0) THEN
				PCS = ST(I)/DBLE(NDyna(I))/ST(1)
				FN1 = F1(I)/DBLE(NDyna(I))
				FN2 = F2(I)/DBLE(NDyna(I))
				WRITE (87, "(10E17.8)") MDTime, PCS, FN1, FN2
			END IF
		END DO
		CLOSE (87)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(NDyna, ttv0)
		DEALLOCATE(RX0, RY0, RZ0)
		DEALLOCATE(ST, F1, F2)
		DEALLOCATE(Pmb, r2, NIT, L, NumN)
	END IF

RETURN
END

!====SUBROUTINE: Get_Pmb--sort particles into descending order according to distance
SUBROUTINE Get_Pmb1
	
	USE SetCOM
	USE SetPCS
	IMPLICIT NONE
	
	INTEGER N, LOGNB2, NN, H, I, J, K, L, Pmbtmp
	DOUBLE PRECISION r2tmp 
!=====================================================
	N = NPart
	
!	====initialize particle index
	DO I = 1, N
		Pmb(I) = I
	END DO

!	====sorting starts
	LOGNB2 = INT(LOG(DBLE(N))/0.69314718 + 1.0D-5)
	NN = N
	DO H = 1, LOGNB2
		NN = NN/2 
		K = N - NN
		DO J = 1, K
			I = J
20			CONTINUE
			L = I + NN
			IF (r2(L)>r2(I)) THEN
				r2tmp = r2(I)
				r2(I) = r2(L)
				r2(L) = r2tmp
				Pmbtmp = Pmb(I)
				Pmb(I) = Pmb(L)
				Pmb(L) = Pmbtmp
				I = I - NN
				IF(I>=1) GOTO 20
			END IF
		END DO
	END DO
	
RETURN
END
