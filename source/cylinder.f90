MODULE CYLINDERUTILS
  ! utilities for dealing with cylindrical sterics
  USE GENUTILS
  USE TRICUBSPLINEUTILS
  IMPLICIT NONE  
  ! cutoff phi value to take advantage of symmetry but avoid cusp at 0 near gimbal lock
  DOUBLE PRECISION, PARAMETER :: PHICUT = -0.1D0
  ! the spline data object for the nucleosome cylinders
  ! and for the segment cylinders
  TYPE(SPDATA3D), POINTER :: NUCCYL,SEGCYL, NUCSEGCYL
  LOGICAL :: NUCCYLSET = .FALSE., SEGCYLSET = .FALSE., NUCSEGCYLSET = .FALSE.

CONTAINS
  SUBROUTINE SETUPCYLARRAYS(DONUCCYL,DOSEGCYL,DONUCSEGCYL,NUCRAD,NUCHEIGHT,SEGRAD,SEGHEIGHT)
    ! set up the various arrays for cylindrical sterics
    ! the parameters are three logical flags for whether each type of
    ! array should be set up
    USE KEYS, ONLY : CYLINDERFILE,CYLINDERFROMFILE,SEGCYLFILE,SEGCYLFROMFILE,NUCSEGCYLFILE,NUCSEGCYLFROMFILE, CYLINDERPTS
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: DONUCCYL, DOSEGCYL, DONUCSEGCYL
    DOUBLE PRECISION, INTENT(IN) :: NUCRAD,NUCHEIGHT,SEGRAD,SEGHEIGHT
    DOUBLE PRECISION :: RADS(2), HEIGHTS(2)

    IF (DONUCCYL) THEN
       ! set up data arrays for cylindrical sterics of nucleosomes
       ALLOCATE(NUCCYL); NUCCYLSET = .TRUE.
       RADS = (/NUCRAD,NUCRAD/); HEIGHTS = (/NUCHEIGHT,NUCHEIGHT/)
       CALL SETUPCYLINDERDATA(NUCCYL,CYLINDERFILE,CYLINDERFROMFILE,RADS,HEIGHTS,CYLINDERPTS)
    ENDIF

    IF (DOSEGCYL) THEN
       ! set up data arrays for cylindrical sterics of DNA segments
       ALLOCATE(SEGCYL); SEGCYLSET = .TRUE.
       RADS = (/SEGRAD,SEGRAD/); HEIGHTS = (/SEGHEIGHT,SEGHEIGHT/)
       CALL SETUPCYLINDERDATA(SEGCYL,SEGCYLFILE,SEGCYLFROMFILE,RADS,HEIGHTS,CYLINDERPTS)
    ENDIF

    IF (DONUCSEGCYL) THEN
       ! set up data arrays for cylindrical sterics of between DNA and nucs
       ALLOCATE(NUCSEGCYL); NUCSEGCYLSET = .TRUE.
       RADS = (/SEGRAD,NUCRAD/); HEIGHTS = (/SEGHEIGHT,NUCHEIGHT/)
       CALL SETUPCYLINDERDATA(NUCSEGCYL,NUCSEGCYLFILE,NUCSEGCYLFROMFILE,RADS,HEIGHTS,CYLINDERPTS)
    ENDIF
  END SUBROUTINE SETUPCYLARRAYS

  SUBROUTINE CLEANUPCYLARRAYS
    ! clean up any previously set up cylinder data arrays
    IMPLICIT NONE

    IF (NUCCYLSET) THEN
       CALL CLEANUPSPLINE(NUCCYL)
       DEALLOCATE(NUCCYL); NULLIFY(NUCCYL); NUCCYLSET = .FALSE.
    ENDIF
    IF (SEGCYLSET) THEN       
       CALL CLEANUPSPLINE(SEGCYL)
       DEALLOCATE(SEGCYL); NULLIFY(SEGCYL); SEGCYLSET = .FALSE.
    ENDIF
    IF (NUCSEGCYLSET) THEN
       CALL CLEANUPSPLINE(NUCSEGCYL)
       DEALLOCATE(NUCSEGCYL); NULLIFY(NUCSEGCYL); NUCSEGCYLSET = .FALSE.
    ENDIF

  END SUBROUTINE CLEANUPCYLARRAYS

  SUBROUTINE SETUPCYLINDERDATA(SPCYL,FILENAME,FROMFILE,RADS,HEIGHTS,NPT)
    ! set up the tabulation for cylindrical sterics
    ! if FROMFILE is true, just read in the spline stuff from that file
    ! otherwise, set it up from scratch and output the data into the file    
    IMPLICIT NONE
    TYPE(SPDATA3D), POINTER :: SPCYL
    LOGICAL, INTENT(IN):: FROMFILE
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: RADS(2), HEIGHTS(2)
    INTEGER, INTENT(IN) :: NPT(3)
    DOUBLE PRECISION :: RADFILE(2), HEIGHTFILE(2)
    TYPE(SPDATA3D), POINTER :: SP
    INTEGER :: N1,N2,N3
    LOGICAL :: READFILE, fileexists
    CHARACTER*10 :: RADSTR1, RADSTR2, HSTR1, HSTR2
    CHARACTER*100 :: PARAMSTR, CYLFILE, CYLFILEOUT

    ! replace the last instance of "?" in the file name with parameter data
    CYLFILE = FILENAME
    WRITE(RADSTR1,'(F10.2)') RADS(1)
    WRITE(RADSTR2,'(F10.2)') RADS(2)
    WRITE(HSTR1,'(F10.2)') HEIGHTS(1)
    WRITE(HSTR2,'(F10.2)') HEIGHTS(2)
    PARAMSTR = 'R'//TRIM(ADJUSTL(RADSTR1))//'H'//TRIM(ADJUSTL(HSTR1))//&
         & 'R'//TRIM(ADJUSTL(RADSTR2))//'H'//TRIM(ADJUSTL(HSTR2))

    CALL REPLACESUBSTR(CYLFILE,'?',TRIM(PARAMSTR))

    READFILE = FROMFILE
    IF (FROMFILE) THEN
       INQUIRE(FILE=CYLFILE,EXIST=FILEEXISTS) 
       READFILE = FILEEXISTS
       IF (.NOT.FILEEXISTS) THEN
          PRINT*, 'Failed to find cylinder file ', trim(adjustl(CYLFILE)), '. Will generate from scratch instead.'
       ENDIF
    ENDIF

    IF (READFILE) THEN
       PRINT*, 'READING CYLINDER DATA FROM FILE:', TRIM(ADJUSTL(CYLFILE))
       OPEN(unit=77,file=CYLFILE,FORM='UNFORMATTED',STATUS='OLD')
       READ(77) N1, N2, N3   
       READ(77) RADFILE,HEIGHTFILE
       PRINT*, 'Cylinder data pts used:', N1, N2, N3
       PRINT '(A,4G10.5)', ' Cylinder radius and height:', RADFILE, HEIGHTFILE     
       IF(ANY(ABS(RADFILE-RADS).GT.10*EPSILON(1D0))) THEN
          PRINT*, 'WARNING: mismatched cylinder radii. Will recalculate cylinder data.', RADFILE,RADS
       ELSEIF (ANY(ABS(HEIGHTFILE-HEIGHTS).GT.10*EPSILON(1D0))) THEN
          PRINT*, 'WARNING: mismatched cylinder heights. Will recalculate cylinder data.', HEIGHTFILE, HEIGHTS
       ELSEIF (NPT(1).NE.N1.OR.NPT(2).NE.N2.OR.NPT(3).NE.N3) THEN
          PRINT*, 'WARNING: wrong number of of points for cylinder tables. Will recalculate cylinder data.'
          PRINT*, NPT
          PRINT*, N1, N2, N3
       ELSE
          CALL SETUPSPLINE(SPCYL,N1,N2,N3)
          READ(77) SPCYL%X1, SPCYL%X2, SPCYL%X3
          READ(77) SPCYL%DATAMAT
          READ(77) SPCYL%SPCOEFF
          READ(77) BINVMAT
          CLOSE(77)
          RETURN
       ENDIF
       CLOSE(77)
    ENDIF

    ! ------------ calculate cylinder data -----------
    CALL SETBINVMAT   

    ! prepare spline arrays
    N1 = NPT(1); N2 = NPT(2); N3 = NPT(3)
    CALL SETUPSPLINE(SPCYL,N1,N2,N3)

    ! tabulate the cylinder distances
    PRINT*, 'TABULATING MINIMAL APPROACH DISTANCES...'
    CALL TABULATECYLDIST(RADS(1),RADS(2),HEIGHTS(1),HEIGHTS(2),1D-8,N1,N2,N3,&
         & SPCYL%DATAMAT,SPCYL%X1,SPCYL%X2,SPCYL%X3)

    ! get all the spline coefficients
    PRINT*, 'SETTING UP CUBIC INTERPOLATION...'
    CALL CREATESPLINE(SPCYL)

    ! ------- output cylinder data to file -----------
    IF (READFILE) THEN
       CYLFILEOUT = 'new-'//CYLFILE
    ELSE
       CYLFILEOUT = CYLFILE
    ENDIF

    PRINT*, 'OUTPUTTING DATA TO FILE:', CYLFILEOUT
    OPEN(unit=77,file=CYLFILEOUT,FORM='UNFORMATTED',STATUS='UNKNOWN')
    write(77) N1, N2, N3  
    WRITE(77) RADS, HEIGHTS
    WRITE(77) SPCYL%X1, SPCYL%X2, SPCYL%X3
    WRITE(77) SPCYL%DATAMAT
    WRITE(77) SPCYL%SPCOEFF
    WRITE(77) BINVMAT
    CLOSE(77)

END SUBROUTINE SETUPCYLINDERDATA

  SUBROUTINE TABULATECYLDIST(RA,RB,LA,LB,TOL,N1,N2,N3,TABLE,X1,X2,X3)
    ! tabulate the minimal center-center distances that avoid cylinder overlap
    ! as a function of cylinder orientations
    ! N1,N2,N3 are the number of point in each dimension (cos(theta), alpha, cos(beta))
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RA, RB, LA, LB, TOL
    INTEGER, INTENT(IN) :: N1, N2, N3
    DOUBLE PRECISION, INTENT(OUT) :: TABLE(N1,N2,N3)
    DOUBLE PRECISION, INTENT(OUT) :: X1(N1),X2(N2),X3(N3)
    INTEGER :: I, J, K

    ! note: alpha only goes from 0 to pi, since results should be symmetric around pi
    print*, 'Tabulating cylindrical sterics for R, L, npt:', RA, LA, RB, LB, N1, N2, N3
    DO I = 1,N1
       X1(I) = DBLE(I-1)/DBLE(N1-1)*2-1
       PRINT*, 'I:', I
       DO J = 1,N2
          X2(J) = DBLE(J-1)/DBLE(N2-1)*2-1
          !PRINT*, 'I,J:', I, J
          DO K = 1,N3
             !PRINT*, 'I,J,K:', I, J,K
             X3(K) = -PI+DBLE(K-1)/DBLE(N3-1)*PI*2
             CALL GETMINCENTDIST(RA,RB,LA,LB,X1(I),X2(J),X3(K),TOL,TABLE(I,J,K))
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE TABULATECYLDIST

  SUBROUTINE GETMINCENTDIST(RA,RB,LA,LB,CTHETA1,CTHETA2,PHI,TOL,DIST)
    ! for a given orientation of two cylinders, get the minimal distance
    ! between their centers such that they do not overlap
    ! the orientation is given by: 
    ! CTHETA1,2 = cos of angle between each cylinder axis and the vector joining their centers
    ! PHI is the dihedral angle between the cylinder axis around the joining vector
    ! TOL is the relative tolerance on how well the distance is determined
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RA, RB, LA, LB, CTHETA1,CTHETA2,PHI, TOL
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    DOUBLE PRECISION :: CA, SA, SB, NA(3), NB(3), CENTA(3), CENTB(3)
    DOUBLE PRECISION :: ST, UPB, LOWB
    LOGICAL :: INTER
    INTEGER, PARAMETER :: MAXITER = 100
    INTEGER :: COUNT

    ! get the cylinder axes
    CA = -SIN(PHI); SA = COS(PHI); SB = SQRT(1-CTHETA2*CTHETA2)    

    CENTA = (/0D0,0D0,0D0/)
    NA = (/SQRT(1-CTHETA1**2),0D0,CTHETA1/)
    NB = (/SB*SA,-SB*CA,-CTHETA2/)

    ! set the upper and lower bound on the distance
    UPB = SQRT((LA/2)**2 + RA**2)  + SQRT((LB/2)**2 + RB**2) + 10*EPSILON(1D0)
    LOWB = MIN(RA,LA/2) + MIN(RB,LB/2)- 10*EPSILON(1D0)   

    CENTB = (/0D0,0D0,LOWB/)
    INTER = CYLINDERINTERSECT(RA,RB,LA,LB,CENTA,NA,CENTB,NB)

    IF (.NOT.INTER) THEN
       PRINT*, 'ERROR IN GETMINCENTDIST: bad lower bound', LOWB, INTER
       print*, 'RA, RB, LA, LB:', RA,RB,LA,LB
       PRINT*, 'THETA, ALPHA, BETA:', cTHETA1,CTHETA2,PHI
       PRINT*, 'CENTA:', CENTA
       PRINT*, 'NA:', NA
       PRINT*, 'CENTB:', CENTB
       PRINT*, 'NB:', NB
       STOP 1
    ENDIF
    CENTB = (/0D0,0D0,UPB/)
    INTER = CYLINDERINTERSECT(RA,RB,LA,LB,CENTA,NA,CENTB,NB)
    IF (INTER) THEN
       PRINT*, 'ERROR IN GETMINCENTDIST: bad upper bound', UPB, INTER
       print*, 'RA, RB, LA, LB:', RA,RB,LA,LB
       PRINT*, 'THETA, ALPHA, BETA:', cTHETA1,CTHETA2,PHI
       PRINT*, 'CENTA:', CENTA
       PRINT*, 'NA:', NA
       PRINT*, 'CENTB:', CENTB
       PRINT*, 'NB:', NB
       STOP 1
    ENDIF

    ! use bisection to find cutoff distance
    DIST = (UPB+LOWB)/2
    COUNT = 0
    DO WHILE ((UPB-LOWB)/DIST.GT.TOL)
       CENTB = (/0D0,0D0,DIST/)       
       INTER = CYLINDERINTERSECT(RA,RB,LA,LB,CENTA,NA,CENTB,NB)
       IF (DIST.LT.MIN(RA,LA/2)+MIN(RB,LB/2).AND..NOT.INTER) THEN
          PRINT*, 'ERROR IN GETMINCENTDIST: cylinders are too close together but not intersecting'
          print*, RA,LA,RB,LB,DIST,INTER
          STOP 1
       ENDIF
       IF (INTER) THEN
          LOWB = DIST
       ELSE
          UPB = DIST
       ENDIF
       DIST = (UPB+LOWB)/2
       COUNT = COUNT + 1
       IF (COUNT.GT.MAXITER) THEN
          PRINT*, 'ERROR IN GETMINCENTDIST: too many iterations', LOWB,UPB
       ENDIF
    END DO

  END SUBROUTINE GETMINCENTDIST

  LOGICAL FUNCTION CYLINDERINTERSECT(RA,RB,LA,LB,CA,NA,CB,NB)
    ! check if two cylinders intersect
    ! RA, RB are the radii
    ! LA,LB are the cylinder lengths
    ! CA,NA are the center and normalized axis of the 1st cylinder
    ! CB,NB are for 2nd cylinder
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RA, RB, LA, LB
    DOUBLE PRECISION, INTENT(IN) :: CA(3), NA(3), CB(3), NB(3)
    DOUBLE PRECISION :: DIST, CENTA(3),CENTB(3), CENTB2(3), TMIN
    INTEGER :: I, J
    LOGICAL :: INSEG1, INSEG2


    CYLINDERINTERSECT = .FALSE.

    ! get distance between line segments
    CALL LINESEGDIST(NA,CA,NB,CB,LA,LB,DIST,INSEG1,INSEG2)

    IF (DIST.GT.(RA+RB)**2) then
       RETURN ! lower bound is above 0       
    ELSE IF (INSEG1.AND.INSEG2) THEN
       CYLINDERINTERSECT = .TRUE. ! shells intersect
       RETURN
    ENDIF


    DO I = -1,1,2
       ! check if axis B crosses any disc on A
       CENTA = CA + I*NA*LA/2
       CYLINDERINTERSECT = LINEDISCINTERSECT(RA,CENTA,NA,NB,CB,LB)       
       IF (CYLINDERINTERSECT) RETURN

       ! check if axis A crosses any disc on B
       CENTB = CB + I*NB*LB/2
       CYLINDERINTERSECT = LINEDISCINTERSECT(RB,CENTB,NB,NA,CA,LA)       

       IF (CYLINDERINTERSECT) RETURN       

       ! check if any pair of discs intersects
       DO J = -1,1,2
          CENTB2 = CB + J*NB*LB/2         

          CYLINDERINTERSECT =  DISCDISCINTERSECT(RA,CENTA,NA,RB,CENTB2,NB)         
          IF (CYLINDERINTERSECT) RETURN
       ENDDO

       ! check for circle-shell intersections
       ! get distance between circle A and axis B      
       CALL CIRCLELINEDIST(RA,CENTA,NA,NB,CB,DIST,TMIN)
       IF (DIST.LT.RB*RB.AND.ABS(TMIN).LT.LB/2) THEN
          CYLINDERINTERSECT = .TRUE.; RETURN
       ENDIF

       ! get distance between circle B and axis A
       CALL CIRCLELINEDIST(RB,CENTB,NB,NA,CA,DIST,TMIN)
       IF (DIST.LT.RA*RA.AND.ABS(TMIN).LT.LA/2) THEN
          CYLINDERINTERSECT = .TRUE.; RETURN
       ENDIF
    ENDDO

  END FUNCTION CYLINDERINTERSECT

  LOGICAL FUNCTION LINEDISCINTERSECT(RD,CD,ND,M,B,L)
    ! check if a line intersects a disc
    ! RD, CD, ND are radius, center and normal of disc
    ! line is M*T+B; length of line is L
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: CD(3), ND(3), M(3), B(3),L,RD
    DOUBLE PRECISION :: T, BC(3), PT(3)

    ! get point where line passes through disc plane
    BC = B-CD
    T = -DOT_PRODUCT(BC,ND)/DOT_PRODUCT(M,ND)

    LINEDISCINTERSECT = .FALSE.
    IF (ABS(T).LE.L/2) THEN
       ! check distance of intersection from disc center
       PT = M*T + BC
       IF (DOT_PRODUCT(PT,PT).LE.RD*RD) THEN
          LINEDISCINTERSECT = .TRUE.
       ENDIF
    ENDIF

  END FUNCTION LINEDISCINTERSECT

  LOGICAL FUNCTION DISCDISCINTERSECT(RA,CA,NA,RB,CB,NB)
    ! check if two discs intersect
    ! assume na, nb are normalized    
    IMPLICIT NONE
    DOUBLE PRECISION :: RA, RB, CA(3), NA(3), CB(3), NB(3)
    DOUBLE PRECISION, PARAMETER :: TINY = 1D-10
    DOUBLE PRECISION :: CBA(3), V(3), U(3), ST, CT, DIFF(3)
    double precision :: DIST, T,PT1(3)

    DISCDISCINTERSECT = .FALSE.
    CBA = CB-CA   

    IF (DOT_PRODUCT(CBA,CBA).GT.(RA+RB)**2) RETURN

    IF (ABS(ABS(DOT_PRODUCT(NA,NB))-1).LT.TINY) THEN
       ! discs are parallel
       DISCDISCINTERSECT = ABS(DOT_PRODUCT(CBA,NA)).LT.TINY           
    ELSE       
       ! consider case where circle A passes through disc B
       CALL CROSS_PRODUCT(NA,NB,U); U = U/NORM(U)
       CALL CROSS_PRODUCT(NA,U,V); 

       ST = DOT_PRODUCT(CBA,NB)/RA/DOT_PRODUCT(V,NB)        
       IF (ABS(ST).GT.1) RETURN ! circle never passes through plane of disc
       CT = SQRT(1-ST**2)

       IF (CT.EQ.0D0) THEN
          ! circle A only hits disc B at one point
          PT1 = RA*CT*U + RA*ST*V -CBA
          DISCDISCINTERSECT = DOT_PRODUCT(PT1,PT1).LE.RB*RB
          RETURN
       ENDIF

       PT1 = -RA*CT*U + RA*ST*V +CA
       DIFF = 2*RA*CT*U      

       ! check if each relevant endpoint on circle A is within the disc
       ! and otherwise, whether some point on the line between them is within the disc
       IF (DOT_PRODUCT(PT1-CB,PT1-CB).LE.RB*RB) THEN
          DISCDISCINTERSECT = .TRUE.; RETURN          
       ELSEIF (DOT_PRODUCT(PT1+DIFF-CB,PT1+DIFF-CB).LE.RB*RB) THEN          
          DISCDISCINTERSECT = .TRUE.; RETURN
       ELSE
          calL ptlinedist(cb,diff,pt1,DIST,T)          
          IF (DIST.LE.RB*RB.AND.T.LT.1.AND.T.GT.0) THEN
             DISCDISCINTERSECT = .TRUE.; RETURN
          ENDIF
       ENDIF

    ENDIF
  END FUNCTION DISCDISCINTERSECT

  SUBROUTINE PTLINEDIST(PT,M,B,DIST,T)
    ! give the squared distance between a point and a line
    ! and the position on the line L = M*t + B where the distance vector hits
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: PT(3), M(3),B(3)
    DOUBLE PRECISION, INTENT(OUT) :: DIST, T
    DOUBLE PRECISION :: DIFF(3)

    T = DOT_PRODUCT(M,PT-B)/DOT_PRODUCT(M,M)
    DIFF = M*T+B-PT
    DIST = DOT_PRODUCT(DIFF,DIFF)

  END SUBROUTINE PTLINEDIST

    SUBROUTINE LINESEGDIST(M,B,U,V,LENA,LENB,DIST,INSEG1,INSEG2)
    ! calculate the minimal squared distance between 2 line segments
    ! line segment 1 has center B, slope M, length LENA
    ! line segment 2 has center V, slope U, length LENB
    ! INSEG is true if the points at minimal separation fall within (not at ends) of each segment
    ! see cylinder-cylinder notes from 10/19/2009

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: M(3), B(3),U(3), V(3),LENA,LENB
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    LOGICAL, INTENT(OUT) :: INSEG1,INSEG2
    DOUBLE PRECISION :: BV(3), UU, UM, MM, MBV, TMP, T, S, DIFF(3)
    DOUBLE PRECISION :: Z, CMP, SMM, SUU

    INSEG1 = .FALSE.; INSEG2 = .FALSE.
    BV = B - V
    UU = DOT_PRODUCT(U,U)
    UM = DOT_PRODUCT(U,M)
    MM = DOT_PRODUCT(M,M)
    MBV = DOT_PRODUCT(M,BV)

    IF (UM.EQ.0D0) THEN
       ! lines are perpendicular
       T = MBV/MM
       S = DOT_PRODUCT(U,BV)/UU
    ELSE
       TMP = (UM-UU*MM/UM)       
       IF (TMP.EQ.0D0) THEN
          ! lines are parallel and do not intersect       
          ! T = -MBV/MM
          ! INSEG2 = .TRUE.
          ! INSEG1 = ABS(T).LT.LENA/2
          ! DIST = DOT_PRODUCT(M*T+BV,M*T+BV)   
          
          SMM = SQRT(MM); SUU = SQRT(UU)
          Z = -MBV/SMM
          CMP = SMM*LENA/2 + SUU*LENB/2

          IF (ABS(Z).LT.CMP) THEN
             ! distance btwn segments is distance btwn lines
             DIST = DOT_PRODUCT(BV,BV)-MM*LENA*LENA/4
             INSEG1 = .TRUE.; INSEG2 = .TRUE.
          ELSE ! segments are offset
             T = SIGN(LENA/2,Z)
             S = SIGN(LENB/2,Z)
             DIFF = M*T+U*S+BV
             DIST = DOT_PRODUCT(DIFF,DIFF)
          ENDIF
          RETURN
       ENDIF

       T= (UU/UM*MBV-DOT_PRODUCT(BV,U))/TMP
       S = (MM*T + MBV)/UM
    ENDIF

    INSEG1 = ABS(T).LT.LENA/2
    IF (.NOT.INSEG1) THEN   
       T = SIGN(LENA/2,T-LENA/2)     
       CALL PTLINEDIST(M*T + B,U,V,DIST,S)
    ENDIF
    
    INSEG2 = ABS(S).LT.LENB/2
    IF (.NOT.INSEG2) THEN
       S = SIGN(LENB/2,S-LENB/2)
       CALL PTLINEDIST(U*S + V,M,B,DIST,T)
    ENDIF
    
    INSEG1 = ABS(T).LT.LENA/2
    IF (.NOT.INSEG1) THEN
       T = SIGN(LENA/2,T-LENA/2) 
    ENDIF   

    DIFF = M*T + B - U*S - V
    DIST = DOT_PRODUCT(DIFF,DIFF)    
  END SUBROUTINE LINESEGDIST

  SUBROUTINE CIRCLELINEDIST(RA,CA,NA,MV,BV,DIST,TMIN)
    ! minimial squared distance between a circle of radius RA, centered at CA
    ! with normal (normalized) given by NA
    ! and the line L(t) = MV*t + BV
    ! also return the point at which minimal distance happens
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: CA(3), NA(3), MV(3), BV(3), RA
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    DOUBLE PRECISION :: A0,A1,A2,A3,A4,A5,A6,ACOEFF(7)
    DOUBLE PRECISION :: DV(3), EV(3), FV(3)
    DOUBLE PRECISION :: AX, BX, CX, TMIN, FA, FD, DFC, DFA,T
    INTEGER :: TC
    DOUBLE PRECISION :: DBRENT
    EXTERNAL DBRENT
    DOUBLE PRECISION :: C1, B0, B1, B2, B3,B4,D0,D1,D2
    DOUBLE PRECISION :: DSCR, SDSCR    

    DV = BV-CA
    EV = MV - DOT_PRODUCT(NA,MV)*NA
    FV = DV - DOT_PRODUCT(NA,DV)*NA

    A6 = DOT_PRODUCT(MV,MV)
    A5 = 2*DOT_PRODUCT(DV,MV)
    A4 = DOT_PRODUCT(DV,DV) + RA**2
    A3 = -2*RA
    A2 = DOT_PRODUCT(EV,EV)
    A1 = 2*DOT_PRODUCT(EV,FV)
    A0 = DOT_PRODUCT(FV,FV)

    ACOEFF = (/A0,A1,A2,A3,A4,A5,A6/)

    ! bracket the interval
    BX = DOT_PRODUCT(CA-BV,MV)/DOT_PRODUCT(MV,MV)
    AX = BX - RA; CX = BX + RA

    ! -------------
    ! get first minimum of function
    ! -------------
    DIST = DBRENT(AX,BX,CX,CLFUNC,10,ACOEFF,1d-7,TMIN)
    call cLfunc(TMIN,acoeff,FD,dfA)    

    ! get coefficients for checking for 2nd minimum
    C1 = FD-A4
    B4 = A6**2
    B3 = 2*A6*A5
    B2=A5**2 - 2*A6*C1 - A3**2*A2
    B1 = -2*A5*C1-A3**2*A1
    B0 = C1**2 - A3**2*A0    
    
    ! deflate to 2nd degree polynomial by dividing by the double root
    D2 = B4
    D1 = B3 + 2*B4*TMIN
    D0 = B2+3*B4*TMIN**2+2*B3*TMIN

    ! find roots of 2nd degree polynomial (if it exists); use this to bracket
    DSCR = D1**2 - 4*D0*D2   
    IF (DSCR.GT.0D0) THEN
       SDSCR = SQRT(DSCR)
       IF (D2.GT.0) THEN
          AX = (-D1 - SDSCR)/(2*D2)
          CX = (-D1 + SDSCR)/(2*D2)
       ELSE
          AX = (-D1 + SDSCR)/(2*D2)
          CX = (-D1 - SDSCR)/(2*D2)
       ENDIF
       BX = (AX+CX)/2       

       ! get 2nd minimum
       DIST = DBRENT(AX,BX,CX,CLFUNC,10,ACOEFF,1d-7,TMIN)
       
    endif

  END SUBROUTINE CIRCLELINEDIST

  SUBROUTINE CLFUNC(T,ACOEFF,F,DF)
    ! get the squared distance btwn a circle and a line using precalculated coefficients
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: T,ACOEFF(:)
    DOUBLE PRECISION, INTENT(OUT) :: F, DF
    DOUBLE PRECISION :: TMP

    TMP = ACOEFF(3)*T**2 + ACOEFF(2)*T+ACOEFF(1)

    IF (TMP.LT.-100*EPSILON(1D0)) THEN
       PRINT*, 'ERROR IN CLFUNC: bad value inside square root', TMP
       STOP
    ELSE
       TMP = SQRT(MAX(0D0,TMP))
    ENDIF

    F = ACOEFF(7)*T**2 + ACOEFF(6)*T + ACOEFF(5) + ACOEFF(4)*TMP
    DF = 2*ACOEFF(7)*T + ACOEFF(6) + ACOEFF(4)/2*(2*ACOEFF(3)*T + ACOEFF(2))/TMP

  END SUBROUTINE CLFUNC
END MODULE CYLINDERUTILS
