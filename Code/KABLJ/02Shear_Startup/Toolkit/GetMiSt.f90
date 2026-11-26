!********************************************************
!====  read pressure tensors generated from LAMMPS   ====
!====       produce stress at specified times        ====
!====         Wensheng Xu, CIAC, June 2020           ====
!********************************************************

!====define global variables and arrays
MODULE SetCOM
	
	INTEGER (KIND=8) Gnum, NComp, NInter, NStep, TStep, Length
	INTEGER (KIND=8), ALLOCATABLE:: StepList(:), Step(:)
	
	DOUBLE PRECISION, PARAMETER:: DelG=0.1, DStrain=0.025

	DOUBLE PRECISION Rate, Delt, Temp
	DOUBLE PRECISION, ALLOCATABLE:: Tensor(:,:), Stress(:)
	
END MODULE

!====main program
PROGRAM MAIN
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER (KIND=8) I
	
	DOUBLE PRECISION Start, Final
!=================================================================
	OPEN (UNIT=12, FILE='LogMiSt')
	
	CALL CPU_TIME(Start)
	CALL PrintDate
	
	OPEN (UNIT=24, FILE='Stress')
	
	CALL INIT

	DO I = 0, Gnum
		TStep = StepList(I)
		CALL Get_Stress
	END DO
	
	CALL FreeMemory
	
	CLOSE (24)
	
	CALL CPU_TIME(Final)
	WRITE (12, "('Total time is', F9.2, ' minutes')") (Final - Start)/60
	CLOSE (12)
	
END

!====INIT: read input parameters and set initial information
SUBROUTINE INIT
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER (KIND=8) I
	
	DOUBLE PRECISION Gmax, TProdt
!================================================================
!	====read basic info
	OPEN (UNIT=15, FILE='../Conf/Basic.txt')
	READ (15, *) Rate, Delt, TProdt
	CLOSE (15)
	
!	====set number of points to calculate
	Gmax = ANINT(Rate*TProdt)
	Gnum = INT(Gmax/DelG)
	
!	====allocate arrays
	ALLOCATE(StepList(0:Gnum))
	
	DO I = 0, Gnum
		StepList(I) = int(I*DelG/(Rate*Delt))
	END DO

!	====get thermo info
	CALL Get_Thermo
	
!	====get all tensors
	CALL Get_Tensor

RETURN
END

!====Get_Thermo: get thermo info
SUBROUTINE Get_Thermo
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER (KIND=8) I, K, NItem, NLine, NProp, MDStep, Tnum
	
	DOUBLE PRECISION MDTime, Strain
	DOUBLE PRECISION, ALLOCATABLE:: Thermo(:)
!================================================================	
!	====read thermo quantities
	OPEN (UNIT=14, FILE='../Conf/InfoThermo')
	READ (14, *) NItem, NLine
	CLOSE (14)
	
	NProp = NItem - 1

!	====allocate arrays
	ALLOCATE(Thermo(NProp))
		
	DO K = 1, NProp
		Thermo(K) = 0
	END DO
	
	Tnum = 0
	Temp = 0.0
	
	OPEN (UNIT=35, FILE='../Conf/Thermo.txt')
	OPEN (UNIT=36, FILE='Thermo')
	READ (35, *)
	DO I = 1, NLine
		READ (35, *) MDStep, Thermo(:)
		MDTime = MDStep*Delt
		Strain = MDTime*Rate
		WRITE (36, "(F16.3, E21.12, 20E17.8)") MDTime, Strain, Thermo(:)
		IF (I>NLine/2) THEN
			Tnum = Tnum + 1
			Temp = Temp + Thermo(1)
		END IF
	END DO
	CLOSE (35)
	CLOSE (36)
	
	Temp = Temp/Tnum
			
!	====free memory
	DEALLOCATE(Thermo)
	
RETURN
END

!====Get_Tensor: get all tensors
SUBROUTINE Get_Tensor
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER (KIND=8) NItem, I, K, Tnum
	
	DOUBLE PRECISION Viscosity
	DOUBLE PRECISION, ALLOCATABLE:: TensorMean(:)
!================================================================	
!	====read thermo quantities
	OPEN (UNIT=14, FILE='../Conf/InfoTensor')
	READ (14, *) NItem, NStep
	CLOSE (14)
	
	NComp = NItem - 1
	
	ALLOCATE(Step(NStep), Tensor(NStep, NComp), Stress(NComp))
	
	WRITE (12, "('  Tensor number              : ', I10)") NStep
	
	OPEN (UNIT=35, FILE='../Conf/Tensor.txt')
	READ (35, *)
	DO I = 1, NStep
		READ (35, *) Step(I), Tensor(I, :)
	END DO
	CLOSE (35)
	
	NInter = Step(2) - Step(1)
	
	WRITE (12, "('  Time steps between tensors : ', I10)") NInter

!	====get Length for averages
	Length = INT(DStrain/(Rate*Delt*NInter))
	
!	====get average values
	ALLOCATE(TensorMean(NComp))
	
	DO K = 1, NComp
		Tnum = 0
		TensorMean(K) = 0.0
		DO I = NStep/2, NStep
			Tnum = Tnum + 1
			TensorMean(K) = TensorMean(K) + Tensor(I, K)
		END DO
		TensorMean(K) = TensorMean(K)/Tnum
	END DO
	
	Viscosity = -TensorMean(1)/Rate

!	====save results
	OPEN (UNIT=21, FILE='Tensor')
	WRITE (21, "(20E17.8)") Temp, Viscosity, TensorMean(:)
	CLOSE (21)

!	====free memory
	DEALLOCATE(TensorMean)
	
RETURN
END

!====Get_Stress: get averaged stress
SUBROUTINE Get_Stress
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER (KIND=8) I, J, K, N, Diff, DMax, StepLo, StepHi
	
	DOUBLE PRECISION MDTime, Strain
!================================================================
	DMax = 1000000
	DO I = 1, NStep
		Diff = ABS(TStep-Step(I))
		IF (Diff<DMax) THEN
			DMax = Diff
			K = I
		END IF
	END DO
	
	DO J = 1, NComp
		Stress(J) = 0.0
		N = 0
		DO I = K - Length, K + Length
			IF (I>=1 .AND. I<=NStep) THEN
				N = N + 1
				Stress(J) = Stress(J) + Tensor(I, J)
				IF (N==1) THEN
					StepLo = Step(I)
					StepHi = Step(I)
				ELSE
					StepHi = Step(I)
				END IF
			END IF
		END DO
		Stress(J) = -Stress(J)/N
	END DO
	
	WRITE (12, *) TStep, K, N, StepLo, StepHi
!	WRITE (*, *) TStep, K, N, StepLo, StepHi
	
	MDTime = TStep*Delt
	Strain = MDTime*Rate
	WRITE (24, "(F16.3, E21.12, 50E17.8)") MDTime, Strain, Stress(:)
	
RETURN
END

!====FreeMemory: free memory for all variables
SUBROUTINE FreeMemory
	
	USE SetCOM
	IMPLICIT NONE
	
	DEALLOCATE(StepList, Step, Tensor, Stress)

RETURN
END

!====Subroutines
INCLUDE "SetDate.inc"
