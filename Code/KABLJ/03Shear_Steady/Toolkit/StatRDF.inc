!====Get_RDF: compute Radial Distribution Function
SUBROUTINE Get_RDF(Switch)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER Switch, NStat, I, J, NBin, Bin, Natom(0:2)
	
	DOUBLE PRECISION Rmax, DelR, dx, dy, dz, r2
	DOUBLE PRECISION RDF00, RDF11, RDF22, RDF12
	DOUBLE PRECISION, ALLOCATABLE:: RBin(:), VBin(:), G00(:)
	DOUBLE PRECISION, ALLOCATABLE:: G11(:), G22(:), G12(:)

	SAVE Natom, NStat, NBin, DelR, RBin, VBin, G00, G11, G22, G12
!=========================================================
IF (Switch==1) THEN
!		====sample
		NStat = NStat + 1
		
!		====whole system
		DO I = 1, NPart
			DO J = I + 1, NPart
				dx = RX(I) - RX(J)
				dy = RY(I) - RY(J)
				dz = RZ(I) - RZ(J)
				dx = dx - Box*ANINT(dx/Box)
				dy = dy - Box*ANINT(dy/Box)
				dz = dz - Box*ANINT(dz/Box)
				r2 = dx**2 + dy**2 + dz**2
				Bin = INT(r2**0.5/DelR)
				IF(Bin<=NBin) THEN
					G00(Bin) = G00(Bin) + 1
				END IF

				IF((Bin<=NBin) .and. (PType(I)==1 .and. PType(J)==1)) THEN
					G11(Bin) = G11(Bin) + 1
					!WRITE(*,*) I, J, PType(I), PType(J), 'Type1'
				END IF

				IF((Bin<=NBin) .and. (PType(I)==1 .and. PType(J)==2)) THEN
					G12(Bin) = G12(Bin) + 1
					!WRITE(*,*) I, J, PType(I), PType(J), 'Type2'
				END IF
					
				IF((Bin<=NBin) .and. (PType(I)==2 .and. PType(J)==2)) THEN
					G22(Bin) = G22(Bin) + 1
					!WRITE(*,*) I, J, PType(I), PType(J), 'Type3'
				END IF

			END DO
		END DO

		
		
	ELSE IF (Switch==0) THEN
!		====initialization
		Natom(0) = NPart
		Natom(1) = NPart1
		Natom(2) = NPart2
		
		NStat = 0
		
!		====set basic info
		Rmax = 0.50D0*Box
		DelR = 0.01D0
		
!		====number of bins
		NBin = INT(Rmax/DelR)
		
!		====allocate arrays
		ALLOCATE(RBin(0:NBin), VBin(0:NBin), G00(0:NBin))
		ALLOCATE(G11(0:NBin), G22(0:NBin), G12(0:NBin))

!		====bin position and shell volume
		DO I = 0, NBin
			RBin(I) = (I + 0.5D0)*DelR
			VBin(I) = 4*PI*DelR**3*(DBLE(I+1)**3 - DBLE(I)**3)/3
			G00(I) = 0.0
			G11(I) = 0.0
			G22(I) = 0.0
			G12(I) = 0.0
		END DO
		
		WRITE (12, "('  To calculate RDF')")
	ELSE IF (Switch==2) THEN
!		====output
		OPEN (UNIT=25, FILE='RDF')
		DO I = 0, NBin - 1
			RDF00 = 2*Vol*G00(I)/(VBin(I)*DBLE(NStat)*Natom(0)*Natom(0))
			RDF11 = 2*Vol*G11(I)/(VBin(I)*DBLE(NStat)*Natom(1)*Natom(1))
			
			IF(Natom(2) == 0)THEN
				RDF22 = 0.0
				RDF12 = 0.0
			ELSE
				RDF22 = 2*Vol*G22(I)/(VBin(I)*DBLE(NStat)*Natom(2)*Natom(2))
				RDF12 = 1*Vol*G12(I)/(VBin(I)*DBLE(NStat)*Natom(1)*Natom(2))
			END IF
				
			WRITE (25, "(10E17.8)") RBin(I), RDF00, RDF11, RDF22, RDF12
		END DO
		CLOSE (25)
	ELSE IF (Switch==3) THEN
!		====free memory
		DEALLOCATE(RBin, VBin, G00)
		DEALLOCATE(G11, G22, G12)
	END IF

RETURN
END
