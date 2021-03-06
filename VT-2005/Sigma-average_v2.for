PROGRAM SIGMAPROFILEV2
!******************************************************************************************
!CREATED USING DIGITAL VISUAL FORTRAN 6.0 (2006)
!
!	This program reads the modified COSMO output file (in text format) from Accelrys' 
! Materials Studio DMol3 and averages the surface segment charge densities per 
! Klamt (1995), Klamt et al (1998), Klamt et al (2000), Lin and Sandler (2002) to 
! establish the segment charges for the "sigma-profile".  This program creates 
! a text file that MS Excel can read and plot.
!
!	THIS PROGRAM WRITTEN BY:
!		RICHARD OLDLAND (roldland@vt.edu)               MIKE ZWOLAK (zwolak@caltech.edu)
!		DEPARTMENT OF CHEMICAL ENGINEERING              PHYSICS DEPARTMENT
!		VIRGINIA TECH                                   CALIFORNIA INSTITUTE OF TECHNOLOGY
!		BLACKSBURG, VA 24060                            PASADENA, CA 91125
!
!	EDITED MARCH 2006 BY:
!		ERIC MULLINS (pmullins@vt.edu)
!		DEPARTMENT OF CHEMICAL ENGINEERING
!		VIRGINIA TECH
!		BLACKSBURG, VA 24060  
!
!
!	VALUES READ FROM THE DATA FILE:
!	ATOM = ATOM NUMBER IN MOLECULE
!	POSXAU = X-CORDINATE OF THE SEGMENT POSITION IN ATOMIC UNITS
!	POSYAU = Y-CORDINATE OF THE SEGMENT POSITION IN ATOMIC UNITS
!	POSZAU = Z-CORDINATE OF THE SEGMENT POSITION IN ATOMIC UNITS
!	POSXA = X-CORDINATE OF THE SEGMENT POSITION IN ANGSTROMS
!	POSYA = Y-CORDINATE OF THE SEGMENT POSITION IN ANGSTROMS
!	POSZA = Z-CORDINATE OF THE SEGMENT POSITION IN ANGSTROMS
!	A = SURFACE SEGMENT AREA; APPROXIMATED AS CIRCULAR (ANGSTROMS SQUARED)
!	CHG = CHARGE OF THE SURFACE SEGMENT (SIGMA, e)
!	SIGMA = RATIO OF SURFACE CHARGE TO AREA (SIGMA/A, e/A**2)
!	POTENT = SURFACE POTENTIAL
!
!	FROM SEGMENT CHARGE AVERAGING:
!	REFF = RADIUS OF THE AREA THAT AFFECTS THE SURFACE SEGMENT CHARGE (ANGSTROMS)
!       DMN = DISTANCE BETWEEN CALCULATED SEGEMENT AND SEGEMENTS AFFECTING IT (ANGSTROMS)
!	RAD = RADIUS OF THE PARTICULAR SURFACE SEGMENT (ASSUMED CIRCULAR, ANGSTROMS SQUARED)
!	SIGMANEW = NEW AVERAGED SIGMA VALUE
!	SIGMASUM = SUMMATION OF EFFECTS FROM OTHER SURFACE SEGMENTS
!	NORMDIST = NORMALIZATION FACTOR
!	NORMSUM = SUMMATION OF ALL NORMALIZATION FACTORS
!
!	REQUIRED INPUT (ON PROMPT):
!		NAME OF FILE (INCLUDING LOCATION AND EXTENSION)
!		NAME OF CHEMICAL (THIS WILL APPEAR IN THE OUTPUT FILE)
!		NUMSEGMENT = THE NUMBER OF SURFACE SEGMENTS
!
!  THE OUTPUT FILE "'CHEMICAL'SIGMA-PROFILE.TXT" IS THE SORTED SIGMA PROFILE
!
!  LITERATURE CITED:
!  Klamt, A. Conductor-like Screening Model for Real Solvents: A New Approach to the
!       Quantitative Calculation of Solvation Phenomena. J. Phys. Chem 1995, 99, 2224.
!  Klamt, A.; Jonas, V.; Burger, T.; Lohrenz, J. Refinement and Parameterization of 
!	COSMO-RS. J. Phys. Chem A 1998, 102, 5074.
!  Klamt, A.; Eckert, F.; COSMO-RS: A Novel and Efficient Method for the a Priori 
!	Prediction of Thermophysical Data of Liquids.  Fluid Phase Equilibria 2000, 
!	172, 43.
!  Lin, S.T.; Sandler, S. A Priori Phase Equilibrium Prediction from a Segment 
!       Contribution Solvation Model. Ind. Eng. Chem. Res, 2002, 41, 899 
!  Lin, S.T.;  Quantum Mechanical Approaches to the Prediction of Phase Equilibria: 
!	Solvation Thermodynamics and Group Contribution Methods, PhD. Dissertation, 
!	University of Delaware, Newark, DE, 2000
!******************************************************************************************
IMPLICIT NONE 

CHARACTER(128):: FILEINDEX, FILEOUTPUT
CHARACTER(16):: CHEMICAL
CHARACTER(1) :: RESPONSE
CHARACTER (256) :: FILENAME
INTEGER :: I, J, K, F, M, N, O, DUMBI, TMP, NUMSEGMENT
INTEGER, DIMENSION (:), ALLOCATABLE :: ATOM
INTEGER, DIMENSION (1500) :: FILEIN, FILEOUT
REAL*8 :: REFF, PI, DUMMY
REAL*8, DIMENSION (:), ALLOCATABLE :: POSXAU, POSYAU, POSZAU, POSXA, POSYA, POSZA, A
REAL*8, DIMENSION (:), ALLOCATABLE :: CHG, SIGMA, POTENT, SIGMANEW, SIGMASUM, RAD, NORMDIST
REAL*8, DIMENSION (:), ALLOCATABLE :: NORMSUM, DMN
REAL*8, DIMENSION(1:51) :: CHGDEN,SP

! ESTABLISH INPUT FILE UNIT NUMBERS
DO N = 1, 1500
   FILEIN(N) = N+20
END DO

! ESTABLISH OUTPUT FILE UNIT NUMBERS
DO O = 1, 1500
   FILEOUT(O) = O+1521
END DO

!ESTABLISH CONSTANTS
PI = 3.14159265358979
REFF = 0.81764200000000

!REPETITION LOOP TO CALCULATE MULTIPLE P(S)
DO F = 1, 1500
    
   ! INPUT VARIABLES DEFINED FROM ARRAYS FILENAME AND NUMBSEGMENT
   !FILELOOP = 'C:\COSMO\VT-2006_AMBER\'//FILENAME(F)
   !FILEOUTPUT = 'C:\PROFILES\VT-2006_AMBER\'//FILENAME(F)
   !NUMSEGMENT = NUMBSEGMENT(F)

   !ESTABLISHING THE COSMO FILE TO READ
   WRITE(*,*) "TYPE THE NAME OF THE FILE YOU WISH TO READ IN,"
   WRITE(*,*) "INCLUDING LOCATION (MAX 256 CHARACTERS), AND HIT ENTER"
   READ (*,*) FILENAME

   !ESTABLISH THE CHEMICAL NAME
   WRITE(*,*) "TYPE IN THE NAME OF THE CHEMICAL (MAX 16 CHARACTERS), AND HIT ENTER"
   WRITE(*,*) "SIGMA PROFILE OUTPUT FILE WILL BE CREATED IN C:\PROFILES\"
   WRITE(*,*) "IF THIS DIRECTORY DOES NOT EXIST, PLEASE CREATE IT AT THIS TIME."
   READ (*,*) CHEMICAL

   FILEOUTPUT = 'PROFILE.TXT'

   !OPEN THE COSMO FILE WITH ALL SEGMENTS AND CORRESPONDING CHARGE DENSITIES
   OPEN(UNIT=FILEIN(F), FILE = FILENAME, STATUS = "OLD", ACTION = "READ", POSITION = "REWIND")

   !ESTABLISH THE NUMBER OF SURFACE SEGMENTS AND ALLOCATE THE ARRAYS
   WRITE(*,*) "TYPE THE NUMBER OF SURFACE SEGMENTS, FROM THE COSMO OUTPUT, AND HIT ENTER"
   READ (*,*) NUMSEGMENT

   ALLOCATE(ATOM(NUMSEGMENT), POSXAU(NUMSEGMENT), POSYAU(NUMSEGMENT), &
        POSZAU(NUMSEGMENT), POSXA(NUMSEGMENT), POSYA(NUMSEGMENT), & 
		POSZA(NUMSEGMENT), A(NUMSEGMENT), CHG(NUMSEGMENT), SIGMA(NUMSEGMENT), &
		POTENT(NUMSEGMENT), SIGMANEW(NUMSEGMENT), SIGMASUM(NUMSEGMENT), &
		RAD(NUMSEGMENT), NORMDIST(NUMSEGMENT), NORMSUM(NUMSEGMENT), DMN(NUMSEGMENT))

   DO I = 1, NUMSEGMENT
      READ(FILEIN(F),*) DUMBI,ATOM(I),POSXAU(I),POSYAU(I),POSZAU(I),CHG(I),A(I),SIGMA(I),POTENT(I),DUMMY
	  !CONVERT THE POSITIONS FROM ATOMIC UNITS TO ANGSTROMS AND ASSIGN NEW ARRAYS
	  POSXA(I) = POSXAU(I) * 0.529177249 
	  POSYA(I) = POSYAU(I) * 0.529177249
	  POSZA(I) = POSZAU(I) * 0.529177249

! CONVERTING THE UNITS FOR RECENT DMOL3 VERSIONS (4.3 AND ABOVE)
	  A(I) = A(I) * 0.529177249**2
	  SIGMA(I) = SIGMA(I)/(0.529177249**2)

	  RAD(I) = SQRT(A(I)/PI)
   
   
   END DO
   
   !CLOSE COSMO FILE
   CLOSE(FILEIN(F))


   !BEGIN AVERAGING SURFACE CHARGES
   DO J=1, NUMSEGMENT
	  SIGMANEW(J) = 0.D0
	  NORMSUM(J)=0.D0
	  
	  DO K=1, NUMSEGMENT
		 DMN(K) = SQRT((POSXA(K)-POSXA(J))**2+(POSYA(K)-POSYA(J))**2+ &
			     (POSZA(K)-POSZA(J))**2)
		 SIGMASUM(K)= SIGMA(K)*(RAD(K)**2*REFF**2)/(RAD(K)**2+REFF**2)* &
			          DEXP(-(DMN(K)**2)/(RAD(K)**2+REFF**2))
		 NORMDIST(K) =(RAD(K)**2*REFF**2)/(RAD(K)**2+REFF**2)* &
			          DEXP(-(DMN(K)**2)/(RAD(K)**2+REFF**2))
		 NORMSUM(J) = NORMSUM(J) + NORMDIST(K)
		 SIGMANEW(J) = SIGMANEW(J) + SIGMASUM(K)
      END DO
	
	SIGMANEW(J) = SIGMANEW(J)/NORMSUM(J)
   END DO

   !CONTAINS AVERAGED SIGMA-PROFILE 
   OPEN (FILEOUT(F), FILE = FILEOUTPUT)

   !SETTING CHGDEN MATRIX
   DO J=1,51
	  SP(J)=0.D0
	  CHGDEN(J) = -0.025D0+0.001D0*DBLE(J-1)
   END DO

   !SIGMA PROFILE SORTING TAKEN FROM LIN DISSERTATION**
   DO J=1,NUMSEGMENT
	  TMP=INT((SIGMANEW(J)-CHGDEN(1))/0.001D0)
	  SP(TMP+1)=SP(TMP+1)+A(J)*(CHGDEN(TMP+2)-SIGMANEW(J))/0.001D0
	  SP(TMP+2)=SP(TMP+2)+A(J)*(SIGMANEW(J)-CHGDEN(TMP+1))/0.001D0
   END DO

   DO J=1,51
	  WRITE(FILEOUT(F),*) CHGDEN(J),SP(J)
   END DO

   CLOSE(FILEOUT(F))

     DEALLOCATE(ATOM, POSXAU, POSYAU, POSZAU, &
    	POSXA, POSYA, POSZA, A, &
	    CHG, SIGMA, POTENT, SIGMANEW, &
	    SIGMASUM, RAD, NORMDIST, NORMSUM, &
	    DMN)

   !REPEAT SIGMA PROFILE CALCULATION FOR ANOTHER COMPOUND
   WRITE (*,*) "DO YOU WISH TO CALCULATE ANOTHER SIGMA PROFILE (Y or N)"
   READ (*,*) RESPONSE

   IF (RESPONSE=="N") EXIT


END DO
 

END PROGRAM SIGMAPROFILEV2
