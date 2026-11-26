!********************************************************
!====      read a DCD file and perform analysis      ====
!====        Wensheng Xu, CIAC, September 2019       ====
!********************************************************

!====SetCOM: define global variables and arrays
MODULE SetCOM
	
	INTEGER NPart, NBead, NPoly
	
	REAL, ALLOCATABLE:: AX(:), AY(:), AZ(:)
	
	DOUBLE PRECISION, PARAMETER:: PI=4.0*ATAN(1.0)

	DOUBLE PRECISION Box, Vol, Rho
	DOUBLE PRECISION, ALLOCATABLE:: RX(:), RY(:), RZ(:)
	DOUBLE PRECISION, ALLOCATABLE:: PX(:), PY(:), PZ(:)
	
END MODULE

!====main program
PROGRAM MAIN
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER NFrame, I, FStatic, Foutput
	
	DOUBLE PRECISION TimeStart, TimeFinal
!=================================================================
	CALL CPU_TIME(TimeStart)
	
	OPEN (UNIT=12, FILE='LogStat')
	CALL PrintDate
	
!	====open DCD file
	OPEN (UNIT=30, FILE='../Conf/CoordL.dcd', STATUS='old', FORM='unformatted')
	
!	====set basic information
	CALL INIT(NFrame)
	
	FStatic = INT(ANINT(DBLE(NFrame)/1000))
	Foutput = INT(ANINT(DBLE(NFrame)/10))
	
	IF (FStatic==0) FStatic = 1

!	====start of calculations
	CALL Get_Static(0)
	
	WRITE (12, "(/, '---->>calculation begins')")
	DO I = 1, NFrame
!		====get confiruration
		CALL Get_Conf
		IF (MOD(I, 10)==0) WRITE (*, *) I, NFrame
!		====static properties
		IF (MOD(I, FStatic)==0) CALL Get_Static(1)
		
!		====write intermediate results
		IF (MOD(I, Foutput)==0) THEN
			CALL Get_Static(2)
			WRITE (12, "('  #frames completed: ', I8)") I
		END IF
	END DO
	CALL Get_Static(2)
	CALL Get_Static(3)
	WRITE (12, "('<<----calculation completed', /)")
	
!	====free memory
	CALL FreeMemory
	
	CALL CPU_TIME(TimeFinal)
	WRITE (12, "('Total time is ', F8.2, ' minutes')") &
		(TimeFinal - TimeStart)/60.0
	
	CLOSE (30)
	CLOSE (12)

END

!====INIT: read input parameters and set initial information
SUBROUTINE INIT(NFrame)
	
	USE SetCOM
	IMPLICIT NONE
	
	INTEGER NFrame, NAtom
!================================================================
!	====read system information
	OPEN (UNIT=15, FILE='../Conf/Basic.txt')
	READ (15, *)
	READ (15, *) NPart, NBead, NPoly
	READ (15, *)
	READ (15, *) Box
	CLOSE (15)
	
	Vol = Box**3
	Rho = NPart/Vol
	
!	====allocate arrays
	CALL SetMemory
	
!	====print basic information
	WRITE (12, "('---->>Basic information')")
	WRITE (12, "('  Total particles        : ', I10)") NPart
	WRITE (12, "('  Chain length           : ', I10)") NBead
	WRITE (12, "('  Chain number           : ', I10)") NPoly
	WRITE (12, "('  Box dimension          : ', F10.4)") Box
	WRITE (12, "('  Number density         : ', F10.4, /)") Rho

!	====read header in the DCD file
	CALL ReadHeader(NFrame, NAtom)
	
	IF (NAtom .NE. NPart) STOP 'NAtom from DCD file is inconsistent with NPart!'
	
	WRITE (12, "('  Number of frames in DCD: ', I10)") NFrame
	WRITE (12, "('  Number of atoms  in DCD: ', I10, /)") NAtom
	
!	====get thermo info
	CALL Get_Thermo
	
RETURN
END

!====Get_Thermo: get thermo info
SUBROUTINE Get_Thermo
	
	USE SetCOM
	IMPLICIT NONE

	INTEGER I, K, NItem, NStep, NProp, Step

	DOUBLE PRECISION Temp, Pres
	DOUBLE PRECISION, ALLOCATABLE:: ThermoI(:), Thermo(:)
!================================================================	
!	====read thermo quantities
	OPEN (UNIT=14, FILE='../Conf/InfoThermo')
	READ (14, *) NItem, NStep
	CLOSE (14)
	
	NProp = NItem - 1

!	====allocate arrays
	ALLOCATE(ThermoI(NProp), Thermo(NProp))
		
	DO K = 1, NProp
		Thermo(K) = 0
	END DO
	
	OPEN (UNIT=35, FILE='../Conf/Thermo.txt')
	READ (35, *)
	DO I = 1, NStep
		READ (35, *) Step, ThermoI(:)
		DO K = 1, NProp
			Thermo(K) = Thermo(K) + ThermoI(K)
		END DO
	END DO
	CLOSE (35)
	
	DO K = 1, NProp
		Thermo(K) = Thermo(K)/NStep
	END DO
	
	Temp = Thermo(1)
	Pres = Thermo(2)
	
!	====print thermo information
	WRITE (12, "('  Data info in Thermo    : ', I10, I12)") NItem, NStep
	WRITE (12, "('  Temperature            : ', F10.4)") Temp
	WRITE (12, "('  Pressure               : ', F10.4, /)") Pres
	
!	====save results
	OPEN (UNIT=21, FILE='Thermo')
	WRITE (21, "(20E17.8)") Rho, Thermo(:)
	CLOSE (21)
	
!	====free memory
	DEALLOCATE(ThermoI, Thermo)
	
RETURN
END

!====Get_Static: call subroutines for computing static properties
SUBROUTINE Get_Static(Switch)

	IMPLICIT NONE
	
	INTEGER Switch
!=================================================================
	CALL Get_RDF(Switch)
	CALL Get_SSF(Switch)
	CALL Get_PSP(Switch)
	CALL Get_MID(Switch)
	CALL Get_LB2(Switch)
	CALL Get_ANG(Switch)

RETURN
END

!====Subroutines
INCLUDE "SetDate.inc"
INCLUDE "SetConf.inc"
INCLUDE "StatRDF.inc"
INCLUDE "StatSSF.inc"
INCLUDE "StatPSP.inc"
INCLUDE "StatMID.inc"
INCLUDE "StatLB2.inc"
INCLUDE "StatANG.inc"
