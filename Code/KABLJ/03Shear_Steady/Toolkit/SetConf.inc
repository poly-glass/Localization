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
	
RETURN
END

!====Get_Conf: get configuration
SUBROUTINE Get_Conf(Frame)
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER Frame, I, P, B
	
	DOUBLE PRECISION Dim(6), Lx, Ly, Lz, XY, XZ, YZ
	DOUBLE PRECISION dx, dy, dz
!================================================================
!	====Dim(1:6) = a, cos(gamma), b, cos(beta), cos(alpha), c
!	====read particle coordinates
	READ (30) (Dim(I), I = 1, 6)
	READ (30) (AX(I), I = 1, NPart)
	READ (30) (AY(I), I = 1, NPart)
	READ (30) (AZ(I), I = 1, NPart)
	
	Lx = Dim(1)
	XY = Dim(3)*Dim(2)
	XZ = Dim(6)*Dim(4)
	Ly = (Dim(3)**2 - XY**2)**0.5
	YZ = (Dim(3)*Dim(6)*Dim(5) - XY*XZ)/Ly
	Lz = (Dim(6)**2 - XZ**2 - YZ**2)**0.5
	
	IF (XY>0.0) THEN
		DelX = XY
	ELSE
		DelX = XY + Ly
	END IF
		
	DO I = 1, NPart
		RX(I) = DBLE(AX(I))
		RY(I) = DBLE(AY(I))
		RZ(I) = DBLE(AZ(I))
	END DO
	
!	====wrap particles in a cubic box
	DO I = 1, NPart
		RX(I) = RX(I) - DelX*ANINT(RY(I)/Ly)
		RX(I) = RX(I) - Lx*ANINT(RX(I)/Lx)
		RY(I) = RY(I) - Ly*ANINT(RY(I)/Ly)
		RZ(I) = RZ(I) - Lz*ANINT(RZ(I)/Lz)
		MX(I) = RX(I)
		MY(I) = RY(I)
		MZ(I) = RZ(I)
	END DO
	
!	====unwrap particles based on the first frame
	IF (Frame==1) THEN
		DO I = 1, NPart
			OX(I) = MX(I)
			OY(I) = MY(I)
			OZ(I) = MZ(I)
		END DO
	ELSE
		DO I = 1, NPart
			dx = MX(I) - OX(I)
			dy = MY(I) - OY(I)
			dz = MZ(I) - OZ(I)
			dx = dx - DelX*ANINT(dy/Ly)
			dx = dx - Lx*ANINT(dx/Lx)
			dy = dy - Ly*ANINT(dy/Ly)
			dz = dz - Lz*ANINT(dz/Lz)
			MX(I) = OX(I) + dx
			MY(I) = OY(I) + dy
			MZ(I) = OZ(I) + dz
 			OX(I) = MX(I)
			OY(I) = MY(I)
			OZ(I) = MZ(I)
		END DO
	END IF
	
	DO I = 1, NPart
		MX(I) = MX(I) + 0.5*Lx
		MY(I) = MY(I) + 0.5*Ly
		MZ(I) = MZ(I) + 0.5*Lz
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
	
	ALLOCATE(AX(NPart), AY(NPart), AZ(NPart))
	ALLOCATE(RX(NPart), RY(NPart), RZ(NPart))
	ALLOCATE(MX(NPart), MY(NPart), MZ(NPart))
	ALLOCATE(OX(NPart), OY(NPart), OZ(NPart))
	ALLOCATE(PType(NPart))

RETURN
END

!====FreeMemory: free memory for all variables
SUBROUTINE FreeMemory
	
	USE SetCOM
	IMPLICIT NONE
	
	DEALLOCATE(AX, AY, AZ)
	DEALLOCATE(RX, RY, RZ)
	DEALLOCATE(MX, MY, MZ)
	DEALLOCATE(OX, OY, OZ)
	DEALLOCATE(PType)

RETURN
END
