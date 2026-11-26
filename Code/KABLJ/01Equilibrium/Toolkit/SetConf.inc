!====ReadHeader: read header in the DCD file
SUBROUTINE ReadHeader(NFrame, NAtom)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER NFrame, NAtom, I, dummyi
	
	REAL dummyr
	
	CHARACTER*4 dummyc
!================================================================
	READ (30) dummyc, NFrame, (dummyi, I=1, 8), dummyc, (dummyi, I=1, 9)
	READ (30) dummyi, dummyr
	READ (30) NAtom
	
	WRITE(*,*) NAtom
RETURN
END

!====Get_Conf: get configuration
SUBROUTINE Get_Conf
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER I, P, B
	INTEGER NumB, NumC
	
	DOUBLE PRECISION D
!================================================================
!	====read particle coordinates
	READ (30) (D, I = 1, 6)
	READ (30) (WX(I), I = 1, NPart)
	READ (30) (WY(I), I = 1, NPart)
	READ (30) (WZ(I), I = 1, NPart)
		
	DO I = 1, NPart
		RX(I) = DBLE(WX(I))
		RY(I) = DBLE(WY(I))
		RZ(I) = DBLE(WZ(I))
	END DO
	
!	====classify particles
	NumB = 0
	NumC = 0
	
	DO I = 1, NPart
		IF (PType(I)==1) THEN
			NumB = NumB + 1
			BX(NumB) = RX(I)
			BY(NumB) = RY(I)
			BZ(NumB) = RZ(I)
		END IF
		IF (PType(I)==2) THEN
			NumC = NumC + 1
			CX(NumC) = RX(I)
			CY(NumC) = RY(I)
			CZ(NumC) = RZ(I)
		END IF
	END DO
	
RETURN
END

!====Set_Part: set particle properties
SUBROUTINE Set_Part
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER I, NumI
!================================================================
!	====read particle types
	DO I = 1, NPart
		IF (I <= NPart1) THEN
			PType(I) = 1
		ELSE
			PType(I) = 2
		END IF
	END DO
		
RETURN
END

!====SetMemory: allocate memory for all variables
SUBROUTINE SetMemory
	
	USE SetCOM
	IMPLICIT NONE
	
	ALLOCATE(WX(NPart), WY(NPart), WZ(NPart))
	ALLOCATE(RX(NPart), RY(NPart), RZ(NPart))
	ALLOCATE(BX(NPart1), BY(NPart1), BZ(NPart1))
	ALLOCATE(CX(NPart2), CY(NPart2), CZ(NPart2))
	ALLOCATE(PType(NPart))

RETURN
END

!====FreeMemory: free memory for all variables
SUBROUTINE FreeMemory
	
	USE SetCOM
	IMPLICIT NONE
	
	DEALLOCATE(WX, WY, WZ)
	DEALLOCATE(RX, RY, RZ)
	DEALLOCATE(BX, BY, BZ)
	DEALLOCATE(CX, CY, CZ)
	DEALLOCATE(PType)

RETURN
END
