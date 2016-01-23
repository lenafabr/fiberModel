MODULE BASINHOPUTILS
  ! utilities for implementing a basin-hopping global optimization search
  USE MT19937, only : GRND
  IMPLICIT NONE

CONTAINS
  SUBROUTINE BASINHOP(CHAINP,DBPT,NHOP,BHTOL,STARTHOPIN)
    ! run a basinhopping calculation for NHOP hops, starting with the given configuration
    ! save results into the given database (dbpt)
    ! bhtol is the minimization tolerance at each basinhopping step
    ! start from hop # starthop and go all the way to NHOP (default starthop=1)
    USE KEYS, ONLY : SAMEECUT, SAMEDISTCUT,NEWDATAFILE, &
         & BHCHECKRANGE,SEGHOPEVERY, BHFACC, BHFSAME, MAXHOPTRY, BHTEMP,&
         & STARTRRANGE, STARTTRANGE, STARTARANGE, OUTFILE
    USE CHAINUTILS, ONLY : CHAIN, CREATECHAIN, COPYCHAIN, CLEANUPCHAIN, OUTPUTCHAIN, SAMECONFIG
    USE DATABASEUTILS, ONLY : DATACHAIN, READDATABASE, DUMPDATABASE, ORDEREDINSERT
    USE QUATUTILS, ONLY : PI
    USE OPTIMIZEUTILS, ONLY : BFGS
    USE HELIXUTILS, ONLY : REGULARIZEHELCOORDS
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(DATACHAIN), POINTER :: DBPT
    INTEGER, INTENT(IN) :: NHOP
    DOUBLE PRECISION, INTENT(IN) :: BHTOL
    INTEGER, INTENT(IN),OPTIONAL  :: STARTHOPIN
    INTEGER :: I, HOP, STPTYPE, TRY, RANK, STARTHOP
    TYPE(CHAIN), TARGET :: PREVCHAIN
    TYPE(CHAIN), POINTER :: PREVCHAINP
    LOGICAL :: SUCCESS, ACCEPT, ISSAME
    DOUBLE PRECISION :: TRANGE(2), RRANGE(2), ARANGE
    INTEGER :: NSAME(2), NACC(2),NTRIAL, NTRIALTYPE(2)
    DOUBLE PRECISION :: FACC(2), FSAME(2), R

    PREVCHAINP=>PREVCHAIN
    TRANGE = STARTTRANGE; RRANGE = STARTRRANGE; ARANGE = STARTARANGE
    NTRIAL = 0; NACC = 0; NSAME = 0; FACC = 0D0; FSAME = 0D0
    NTRIALTYPE = 0;

    IF (PRESENT(STARTHOPIN)) THEN
       STARTHOP = STARTHOPIN
    ELSE
       STARTHOP = 1
    ENDIF

    CALL CREATECHAIN(PREVCHAINP)

    PRINT '(A,6G20.10)', 'Initial chain coordinates prior to basinhopping: ', CHAINP%VEC(1:6)

    ! minimize the chain
    CALL BFGS(CHAINP,BHTOL,SUCCESS)       

    IF (SUCCESS) THEN
       ! standardize the coordinate ranges (no negative height, etc)
       CALL REGULARIZEHELCOORDS(CHAINP%VEC(1:6))
       CALL ORDEREDINSERT(DBPT,CHAINP,RANK,.TRUE.,SAMEECUT, SAMEDISTCUT)
       IF (RANK.GT.0) THEN
          PRINT*, 'Inserting chain config into database with rank ', RANK
          ! output updated database into file
          CALL DUMPDATABASE(DBPT,NEWDATAFILE,0)          
       ENDIF
    ELSE       
       PRINT*, 'BASINHOP: failed to minimize chain config all the way'
    ENDIF

    CALL COPYCHAIN(CHAINP,PREVCHAINP,.TRUE.)

    DO HOP = STARTHOP,NHOP

       STPTYPE = 1 ! steps that move nucleosomes
       IF (SEGHOPEVERY.GT.0) THEN
          IF (MOD(HOP-STARTHOP+1,SEGHOPEVERY).EQ.0) THEN
             STPTYPE = 2 ! steps that move linker dna only
          ENDIF
       ENDIF

       DO TRY = 1,MAXHOPTRY
          ! take a MC step in regular helix coordinates
          CALL HELIXMCSTEP(CHAINP,STPTYPE,RRANGE(STPTYPE),TRANGE(STPTYPE),ARANGE)

          print '(A,I5,A,6G12.5)', 'Helix parameters after MC hop ', HOP, ': ', CHAINP%VEC(1:6)
          ! Minimize the chain
          CALL BFGS(CHAINP,BHTOL,SUCCESS)
          IF (SUCCESS) THEN
             EXIT
          ELSE
             PRINT*, 'FAILED to minimize chain all the way. Try again.', TRY
             CHAINP%VEC = PREVCHAINP%VEC
          ENDIF
       ENDDO

       IF (TRY.GT.MAXHOPTRY) THEN
          PRINT*, 'PROBLEM IN BASINHOP: failed to minimize chain for too many tries'
          STOP 1
       ENDIF

       ! decide whether or not to accept the step
       ACCEPT = .FALSE.
       IF (CHAINP%ENERGY.LT.PREVCHAINP%ENERGY) THEN
          ACCEPT = .TRUE.
       ELSE
          R = GRND()
          IF (R < EXP((PREVCHAINP%ENERGY-CHAINP%ENERGY)/BHTEMP)) THEN
             ACCEPT = .TRUE.
          ENDIF
       ENDIF

       ! check if this is the same as the previous minimum
       CALL SAMECONFIG(CHAINP,PREVCHAINP,SAMEECUT, SAMEDISTCUT,ISSAME)

       IF (ISSAME) THEN  
          print*, 'New minimum is same as previous'
          ! Number of steps that failed to change minimum
          NSAME(STPTYPE) = NSAME(STPTYPE) + 1
       ENDIF
       IF (ACCEPT) THEN
          NACC(STPTYPE) = NACC(STPTYPE) + 1             
       END IF
       NTRIAL = NTRIAL + 1     
       NTRIALTYPE(STPTYPE) = NTRIALTYPE(STPTYPE) + 1

       PRINT '(A,L1,A,2G15.5)', 'Accepting chain? ', ACCEPT, ' Previous energy, new energy: ', &
            & PREVCHAINP%ENERGY, CHAINP%ENERGY

       ! standardize the coordinate ranges (no negative height, etc)
       CALL REGULARIZEHELCOORDS(CHAINP%VEC(1:6))
       ! insert configuration into database
       CALL ORDEREDINSERT(DBPT,CHAINP,RANK,.TRUE.,SAMEECUT, SAMEDISTCUT)       
       IF (RANK.GT.0) THEN
          PRINT*, 'Inserting chain config into database with rank ', RANK
          ! output updated database into file
          CALL DUMPDATABASE(DBPT,NEWDATAFILE,HOP)

          IF (RANK.EQ.1) THEN
             PRINT*, 'New lowest energy structure: ', CHAINP%ENERGY, ' Outputting to file: ', TRIM(OUTFILE)
             CALL OUTPUTCHAIN(CHAINP,OUTFILE)
          ENDIF
       ENDIF

       IF (ACCEPT) THEN
          CALL COPYCHAIN(CHAINP,PREVCHAINP,.FALSE.)
       ELSE
          CALL COPYCHAIN(PREVCHAINP,CHAINP,.FALSE.)
       ENDIF

       FSAME = DBLE(NSAME)/NTRIALTYPE
       FACC = DBLE(NACC)/NTRIALTYPE

       PRINT '(A,I6,6G12.5)', 'HOP, MINE, FSAME, FACC, RRANGE,TRANGE,ARANGE:', HOP, &
            & DBPT%MINE, FSAME(1), FACC(1), RRANGE(1), TRANGE(1), ARANGE


       ! update ranges 
       IF (BHCHECKRANGE > 0 .and. MOD(HOP-STARTHOP+1,BHCHECKRANGE).EQ.0) THEN 
          ! update ranges aiming for a specific frequency of stepping to new minima
          DO I = 1,2          
             IF (NTRIALTYPE(I).GT.0) THEN                          
                IF (FSAME(I).LT.BHFSAME-0.1)  THEN
                   TRANGE(I) = TRANGE(I)/2
                   IF (I.EQ.1) ARANGE = ARANGE /2
                   RRANGE(I) = RRANGE(I)/2
                ELSEIF (FSAME(I).GT.BHFSAME+0.1) THEN
                   TRANGE(I) = MIN(TRANGE(I)*2,2*PI)
                   IF (I.EQ.1) ARANGE = MIN(ARANGE*2,2D0)
                   RRANGE(I) = RRANGE(I)*2                
                ENDIF
             ENDIF
          ENDDO

          ! update ranges aiming for a specific frequency of accepted structures
          DO I = 1,2
             IF (NTRIALTYPE(I).GT.0) THEN                          
                IF (FACC(I).LT.BHFACC-0.1)  THEN
                   TRANGE(I) = TRANGE(I)/2
                   IF (I.EQ.1) ARANGE = ARANGE /2
                   RRANGE(I) = RRANGE(I)/2
                ELSEIF (FACC(I).GT.BHFACC+0.1) THEN
                   TRANGE(I) = MIN(TRANGE(I)*2,2*PI)
                   IF (I.EQ.1) ARANGE = MIN(ARANGE *2,2D0)
                   RRANGE(I) = RRANGE(I)*2                
                ENDIF
             ENDIF
          ENDDO

          ! reset counters 
          NTRIALTYPE = 0; NACC = 0; NSAME = 0; NTRIAL = 0
       ENDIF

    ENDDO

    CALL CLEANUPCHAIN(PREVCHAINP)
  END SUBROUTINE BASINHOP

  SUBROUTINE HELIXMCSTEP(CHAINP,STPTYPE,RRANGE,TRANGE,ARANGE)
    ! take a monte carlo step in the space of helix coordinates
    ! STPTYPE=1 means move the regular helix coordinates; 
    ! redistribute linker beads if the end-to-end vector changes by more than a segment length 
    ! STPTYPE=2 means move the linker beads and angles only
    ! RRANGE is a range for distance hops 
    ! TRANGE is a range for angle hops (max 2pi)
    ! ARANGE is a range for cos(angle) hops (max 2)
    USE QUATUTILS
    USE KEYS, ONLY : FIXHELIXPARAM, FLIPNUCMC, FLIPNUCPROB
    USE CHAINUTILS, ONLY : CHAIN
    USE HELIXUTILS, ONLY : FROMHELIXREP
    USE INPUTSTRUCTS, ONLY : INTERPOLATEBEADS
    USE GENUTILS, ONLY : RANDOMAXIS
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: STPTYPE
    DOUBLE PRECISION, INTENT(IN) :: RRANGE,TRANGE,ARANGE
    DOUBLE PRECISION :: RND, MINT, MAXT, MINH, MAXH, MINR, MAXR, SAVEEUL(3)
    DOUBLE PRECISION :: LINKP(3), LINKM(3), LVECOLD(3), LVECNEW(3), LVECDIFF(3)
    DOUBLE PRECISION :: RANDAX(3), THETA
    INTEGER :: I,J,N
    TYPE(QUATERNION) :: QROT, QSTART
    LOGICAL :: REDISTRIBUTEBEADS

    IF (STPTYPE.EQ.1) THEN
       ! change regular helix coordinates
       ! find the vector from the start to end of the linker  
       CALL FROMHELIXREP(CHAINP)
       LINKP = QUAT2PT(CHAINP%BENDS(2)%ORIENT*PTQUAT(CHAINP%BENDS(2)%POSM)/CHAINP%BENDS(2)%ORIENT)
       LINKM  = QUAT2PT(CHAINP%BENDS(1)%ORIENT*PTQUAT(CHAINP%BENDS(1)%POSP)/CHAINP%BENDS(1)%ORIENT)
       LVECOLD = LINKP-LINKM

       ! change the height per nucleosome
       IF (.NOT.FIXHELIXPARAM(1)) THEN
          MINH = MAX(CHAINP%VEC(1)-RRANGE/2,EPSILON(1D0))
          MAXH = CHAINP%VEC(1)+RRANGE/2
          CHAINP%VEC(1) = GRND()*(MAXH-MINH)+MINH
       ENDIF

       ! change the angle per nucleosome
       IF (.NOT.FIXHELIXPARAM(2)) THEN
          MINT = CHAINP%VEC(2)-TRANGE/2
          MAXT = CHAINP%VEC(2)+TRANGE/2
          RND = GRND()*(MAXT-MINT)+MINT
          CHAINP%VEC(2) = ANGLE2PI(RND)
       ENDIF

       ! change the radius
       IF (.NOT.FIXHELIXPARAM(3)) THEN
          MINR = MAX(0D0,CHAINP%VEC(3) - RRANGE/2)
          MAXR = CHAINP%VEC(3) + RRANGE/2        
          CHAINP%VEC(3) = GRND()*(MAXR-MINR) + MINR                    
       END IF

       SAVEEUL = CHAINP%VEC(4:6)
       ! reorient nucleosomes around random axis
       CALL RANDOMAXIS((/0D0,0D0,1D0/),ARANGE,RANDAX)
       THETA = GRND()*TRANGE-TRANGE/2
       QROT = ROTQUAT(THETA,RANDAX)
       CALL EULER2QUAT(CHAINP%VEC(4:6),QSTART)
       CALL QUAT2EULER(QSTART*QROT,CHAINP%VEC(4:6))
       DO I = 1,3
          IF (FIXHELIXPARAM(I+3)) CHAINP%VEC(I+3) = SAVEEUL(I)
       ENDDO

       IF (FLIPNUCMC) THEN
          ! randomly decide whether to flip nucleosomes upside down
          IF (GRND() < FLIPNUCPROB.AND..NOT.ANY(FIXHELIXPARAM(4:6))) THEN
             print*, 'Flipping nucleosomes', CHAINP%VEC(1:6)
             CHAINP%VEC(4) = PI+CHAINP%VEC(4)
             CHAINP%VEC(5) = PI - CHAINP%VEC(5)
             CHAINP%VEC(6) = PI - CHAINP%VEC(6)
          ENDIF
       ENDIF

       ! find the vector from the start to end of the linker
       CALL FROMHELIXREP(CHAINP)
       LINKP = QUAT2PT(CHAINP%BENDS(2)%ORIENT*PTQUAT(CHAINP%BENDS(2)%POSM)/CHAINP%BENDS(2)%ORIENT)
       LINKM  = QUAT2PT(CHAINP%BENDS(1)%ORIENT*PTQUAT(CHAINP%BENDS(1)%POSP)/CHAINP%BENDS(1)%ORIENT)
       LVECNEW = LINKP-LINKM

       ! redistribute beads  only if vector from start to end of linker
       ! has moved more than a segment length
       LVECDIFF = LVECNEW-LVECOLD
       REDISTRIBUTEBEADS = DOT_PRODUCT(LVECDIFF,LVECDIFF).GT.CHAINP%LS**2

       IF (REDISTRIBUTEBEADS) THEN         
          !print*, 'redistributing beads...', DOT_PRODUCT(LVECDIFF,LVECDIFF)
          CALL INTERPOLATEBEADS(CHAINP)
          CALL FROMHELIXREP(CHAINP)
       ENDIF
    ELSEIF (STPTYPE.EQ.2) THEN
       ! randomly move linker segments
       DO I = 1,CHAINP%NSEG
          ! change segment twist
          CHAINP%VEC(4*I+3) = CHAINP%VEC(4*I+3) + GRND()*TRANGE-TRANGE/2
          ! change bead position
          IF (I.LT.CHAINP%NSEG) THEN
             DO J = 1,3
                CHAINP%VEC(4*I+3+J) = CHAINP%VEC(4*I+3+J) + GRND()*RRANGE-RRANGE/2
             ENDDO
          ENDIF
       ENDDO
    ELSE
       PRINT*, 'ERROR IN HELIXMCSTEP: STPTYPE must be 1 or 2. Actual value: ', STPTYPE
       STOP 1
    ENDIF
  END SUBROUTINE HELIXMCSTEP

END MODULE BASINHOPUTILS
