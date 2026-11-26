!====PrintDate: print current time
SUBROUTINE PrintDate

	IMPLICIT NONE

	INTEGER DateInfo(8)
	
	CHARACTER*30 Date, Time, Zone 
!==================================================
	CALL DATE_AND_TIME(Date, Time, Zone, DateInfo)
	WRITE (12, "('Date & Time: ', I0.4, '-', I0.2, &
		'-', I0.2, ' ', I0.2, ':', I0.2, /)") &
		DateInfo(1), DateInfo(2), DateInfo(3), DateInfo(5), DateInfo(6)

RETURN
END
