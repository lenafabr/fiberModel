MODULE INPUTSTRUCTS
  ! utilities to 
  ! create a chain structure from a previouly saved configuration
  ! and for interpolating a linker between 2 bends
  IMPLICIT NONE

CONTAINS
  SUBROUTINE READOUTFILE(CHAINP,INFILE,USEBEADS)
    ! Set the structure of CHAINP based on coordinates from an output file    
    ! Will only read files generated using OUTPUTCHAIN in module CHAINUTILS
    ! Only reads positions of the first 2 bends and the linker between them
    ! Does NOT read in tail positions or filler beads    
    ! Only sets the BEADS / BEADQ / BEND%CENT / BEND%ORIENT coordinates (for first linker and first two bends only)
    ! To convert to the single vector regular helix coordinates used in optimization, run GETHELIXREP after this

    ! if USEBEADS=0, then the beads are interpolated between bends without reference to the input file
    ! if USEBEADS = 1, then the new beads are interpolated between the old beads (can have different number of beads; but even if the same number is used, the spacing between beads will be reregularized, so the exact energy may change a little)
    ! if USEBEADS = 2, then the exact positions of the beads in the file are used (WARNING: must have same current number of beads as in the file)
    USE QUATUTILS
    USE GENUTILS, ONLY : NORM, NORMALIZE, CROSS_PRODUCT
    USE CHAINUTILS, ONLY : CHAIN, CREATECHAIN, BEND, PLACETAILBEADS
    USE INPUTPARAMS, ONLY : READLINE, READU, READF

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    CHARACTER(LEN=*),INTENT(IN) :: INFILE
    INTEGER, INTENT(IN) :: USEBEADS
    CHARACTER*100 :: STR
    DOUBLE PRECISION :: TMP(4),DUMMY1(3),DUMMY2(3)
    INTEGER :: B, I, NUMBENDS, NPT, TAILCT, NITEMS, MAXPT, IND
    DOUBLE PRECISION, ALLOCATABLE :: SPT(:),XPT(:,:), BRANCHPT(:,:)
    TYPE(QUATERNION), ALLOCATABLE :: QPT(:)
    LOGICAL :: NEWBEND, TAILFOUND, REACHEDBEND, FILEEND
    TYPE(BEND), POINTER :: BP
    DOUBLE PRECISION :: ZAX(3), NZ, ROTMAT(3,3)

    IF (.NOT.CHAINP%ARRAYSET) CALL CREATECHAIN(CHAINP)    

    print*, 'Reading structure information from file: ', trim(INFILE), USEBEADS

    TAILFOUND = .FALSE.

    ! read in the bend information from the file
    ! and find the number of beads as well
    OPEN(UNIT=99,FILE=INFILE,STATUS='OLD')        
    STR=''
    B = 1; MAXPT = 0
    NEWBEND = .FALSE.
    NUMBENDS = 0
    TAILCT = 0
    DO        
       CALL READLINE(99,FILEEND,NITEMS)
       IF (NITEMS.EQ.0) EXIT
       CALL READU(STR)
       DO I = 1,3
          CALL READF(DUMMY1(I))
       ENDDO
       DO I = 1,NITEMS-4
          CALL READF(TMP(I))
       ENDDO
      
       IF (FILEEND) EXIT ! end of file reached
       IF (STR.EQ.'B') THEN          
          IF (B.LE.2) THEN
             CHAINP%BENDS(B)%CENT = DUMMY1
             CHAINP%BENDS(B)%ORIENT%W  = TMP(1)
             CHAINP%BENDS(B)%ORIENT%X  = TMP(2)
             CHAINP%BENDS(B)%ORIENT%Y  = TMP(3)
             CHAINP%BENDS(B)%ORIENT%Z  = TMP(4)                  
             B = B + 1   
             TAILCT =0
          ELSE
             EXIT
          ENDIF
       ELSEIF (STR.EQ.'T'.and.CHAINP%BENDS(1)%NTAIL.GT.0) THEN
          BP=>CHAINP%BENDS(B-1)
          ! flexible tail
          TAILCT = TAILCT + 1
          BP%FLEXTAILS(3*(TAILCT-1)+1:3*TAILCT) = QUAT2PT(INVQUAT(BP%ORIENT)*PTQUAT(DUMMY1-BP%CENT)*BP%ORIENT) 
          TAILFOUND = .TRUE.
       ELSEIF (STR.EQ.'A') THEN
          IF (NEWBEND) NUMBENDS = NUMBENDS + 1
          IF (NUMBENDS.EQ.1) MAXPT = MAXPT + 1
          NEWBEND = .FALSE.
       ELSE
          NEWBEND = .TRUE.
       ENDIF
    END DO
    CLOSE(99)
    
    IF (CHAINP%NBEND.LT.2) THEN
       CALL PLACETAILBEADS(CHAINP)
       RETURN
    ENDIF
    
    ! XPT contains bead positions from the file
    ! SPT contains position along chain for each bead position in the file 
    ! (cumulative distance in the piecewise linear path)
    ! BRANCHPT has the branch positions for each linker bead
    ! QPT has the quaternions for each input segment

    ALLOCATE(XPT(MAXPT,3),SPT(MAXPT),QPT(0:MAXPT),BRANCHPT(MAXPT+1,3))
    
    ! read through file again to get bead coordinates
    OPEN(UNIT=99,FILE=INFILE,STATUS='OLD')       
    REACHEDBEND = .FALSE.
    IF (USEBEADS.EQ.1.OR.USEBEADS.EQ.2) THEN
       B = 0 ! current bend
       CALL READLINE(99,FILEEND,NITEMS)
       CALL READU(STR)
       DO I = 1,3
          CALL READF(DUMMY1(I))
       ENDDO
       DO I = 1,3
          CALL READF(DUMMY2(I))
       ENDDO

       DO WHILE (.NOT.REACHEDBEND.OR.STR.NE.'A')
          CALL READLINE(99,FILEEND,NITEMS)
          CALL READU(STR)
          DO I = 1,3
             CALL READF(DUMMY1(I))
          ENDDO
          DO I = 1,3
             CALL READF(DUMMY2(I))
          ENDDO
          
          REACHEDBEND = REACHEDBEND.OR.STR.NE.'A'
       ENDDO      
       
       BRANCHPT(1,:) = DUMMY2

       NPT = 0
       DO WHILE (STR.EQ.'A')
          NPT = NPT + 1          
          XPT(NPT,:) = DUMMY1; 
          IF (NPT.EQ.1) THEN
             SPT(1) = 0D0
          ELSE
             ZAX = XPT(NPT,:)-XPT(NPT-1,:); NZ = NORM(ZAX)
             SPT(NPT) = SPT(NPT-1) + NZ
          ENDIF

          ! get the quaternion
          IF (NPT.GT.1) THEN
             ROTMAT(:,3) = ZAX/NZ
             ROTMAT(:,1) = BRANCHPT(NPT-1,:)-XPT(NPT-1,:)
             CALL NORMALIZE(ROTMAT(:,1))             
             CALL CROSS_PRODUCT(ROTMAT(:,3),ROTMAT(:,1),ROTMAT(:,2))
             QPT(NPT-1) = ROTMAT2QUAT(ROTMAT)             
!             PRINT*, 'TESTXA:', NPT-1, QPT(NPT-1)
          ENDIF
          CALL READLINE(99,FILEEND,NITEMS)
          CALL READU(STR)
          DO I = 1,3
             CALL READF(DUMMY1(I))
          ENDDO
          DO I = 1,3
             CALL READF(BRANCHPT(NPT+1,I))
          ENDDO
       ENDDO
    ENDIF
    
    QPT(0) = CHAINP%BENDS(1)%ORIENT*CHAINP%BENDS(1)%TP
    QPT(NPT) = CHAINP%BENDS(2)%ORIENT*CHAINP%BENDS(2)%TM

    IF (USEBEADS.EQ.0) THEN
       CALL INTERPOLATEBEADS(CHAINP)
    ELSEIF (USEBEADS.EQ.1) THEN
       ! interpolate beads on the chain linker along piecewise path set in the configuration file
       CALL PLACEBEADSONPATH(CHAINP, MAXPT, XPT, QPT)
       CALL PLACETAILBEADS(CHAINP)
    ELSEIF(USEBEADS.EQ.2) THEN
       ! place chain linker beads exactly as set in file
       ! tail beads are still placed in straight line
       IF (MAXPT.NE.CHAINP%NSEG+1) THEN
          PRINT*, 'ERROR IN READOUTFILE: trying to use linker beads directly from file, &
               & but number of beads does not match number of linker segments in chain.', &
               & MAXPT, CHAINP%NSEG
          STOP 1
       ENDIF

       DO I = 1,CHAINP%NSEG
          IND = CHAINP%BENDIND(1)+I
          CHAINP%BEADS(IND,:) = XPT(I,:)          
          CHAINP%BEADQ(IND) = QPT(I)          
       ENDDO

       IND = CHAINP%BENDIND(2)
       CHAINP%BEADS(IND,:) = XPT(MAXPT,:)      
       CHAINP%BEADQ(IND) = CHAINP%BENDS(2)%ORIENT*CHAINP%BENDS(2)%TM
       CALL PLACETAILBEADS(CHAINP)
    ENDIF

    IF (ALLOCATED(XPT)) THEN
       DEALLOCATE(XPT,SPT,QPT)
    ENDIF

    CLOSE(99)
  END SUBROUTINE READOUTFILE

  SUBROUTINE PLACEBEADSONPATH(CHAINP, NPT, XPT, QPT)
    ! given a piecewise linear path of NPT beads at positions XPT
    ! and with orientations QPT
    ! place the beads or the chain along this path
    ! note that the number of beads in the chain does not have to
    ! equal NPT
    USE KEYS, ONLY : BPLEN
    USE CHAINUTILS, ONLY : CHAIN
    USE QUATUTILS

    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: NPT
    DOUBLE PRECISION, INTENT(IN) :: XPT(NPT,3)
    TYPE(QUATERNION), INTENT(IN) :: QPT(0:NPT)
    INTEGER :: I, B, IND, PIECE(CHAINP%NSEG)
    DOUBLE PRECISION :: CUMLEN(NPT), LENPERSEG, PATHCRD, FRAC, DIFF(3)
    DOUBLE PRECISION :: PREVBRANCH, CURBRANCH, TWISTPERLEN, ANG, EULREL(3)
    DOUBLE PRECISION :: CURTWIST, NEXTBRANCH
    TYPE(QUATERNION) :: QREL, QROT, QCUR

    ! cumulative length along the path
    CUMLEN(1) = 0D0
    DO I = 1,NPT-1
       DIFF = XPT(I+1,:)-XPT(I,:)
       CUMLEN(I+1) = CUMLEN(I) + SQRT(DOT_PRODUCT(DIFF, DIFF))
    ENDDO
    ! path length per bead
    LENPERSEG = CUMLEN(NPT)/CHAINP%NSEG

    ! interpolate the bead placement
    IND = CHAINP%BENDIND(1)
    PIECE(1) = 1
    DO B = 1,CHAINP%NSEG
       IF (B.GT.1) PIECE(B) = PIECE(B-1)

       ! coordinate along path
       PATHCRD = LENPERSEG*(B-1)

       ! which piece of the path are we on?
       DO WHILE (CUMLEN(PIECE(B)+1).LT.PATHCRD)
          IF (PIECE(B).GT.NPT-1) THEN
             PRINT*, 'ERROR IN PLACEBEASONPATH: desired chain length is outside of path length'
             STOP 1
          ENDIF
          PIECE(B) = PIECE(B) + 1
       ENDDO

       ! fraction of the way along the piece
       FRAC = (PATHCRD-CUMLEN(PIECE(B)))/(CUMLEN(PIECE(B)+1)-CUMLEN(PIECE(B)))

       ! place the bead at the appropriate position
       CHAINP%BEADS(IND+B,:) = XPT(PIECE(B),:) + FRAC*(XPT(PIECE(B)+1,:)-XPT(PIECE(B),:))       
    ENDDO

    ! final bead position and orientation
    IND = CHAINP%BENDIND(2)
    CHAINP%BEADS(IND,:) = XPT(NPT,:)
    CHAINP%BEADQ(IND) = CHAINP%BENDS(2)%ORIENT*CHAINP%BENDS(2)%TM

    ! interpolate the twist orientations
    IND = CHAINP%BENDIND(1)
    DO B = 1,CHAINP%NSEG
       ! coordinate along path
       PATHCRD = LENPERSEG*(B-1)+LENPERSEG/2

       ! position of previous branch to current piece branch
       ! (each branch is in center of segment)       
       IF (PIECE(B).EQ.1) THEN
          PREVBRANCH = -BPLEN/2
       ELSE
          PREVBRANCH = (CUMLEN(PIECE(B))+CUMLEN(PIECE(B)-1))/2
       ENDIF
       CURBRANCH = (CUMLEN(PIECE(B)+1)+CUMLEN(PIECE(B)))/2

       IF (PIECE(B).EQ.NPT-1) THEN
          NEXTBRANCH = CUMLEN(NPT)+BPLEN/2
       ELSE
          NEXTBRANCH = (CUMLEN(PIECE(B)+2)+CUMLEN(PIECE(B)+1))/2
       ENDIF

       ! find the twist per length in this region
       IF (PATHCRD.LT.CURBRANCH) THEN
          QREL = INVQUAT(QPT(PIECE(B)-1))*QPT(PIECE(B))  
          CALL QUAT2EULER(QREL, EULREL)
          TWISTPERLEN = ANGLE2PI(EULREL(1)+EULREL(3))/(CURBRANCH-PREVBRANCH)
       ELSE
          QREL = INVQUAT(QPT(PIECE(B)))*QPT(PIECE(B)+1)
          CALL QUAT2EULER(QREL, EULREL)
          TWISTPERLEN = ANGLE2PI(EULREL(1)+EULREL(3))/(NEXTBRANCH-CURBRANCH)
       ENDIF

       !print*, 'TESTX2:', B, PIECE(B), PATHCRD, TWISTPERLEN, NEXTBRANCH-CURBRANCH, CURBRANCH
       ! set up orientation
       CALL COORDS2QUAT(0D0,CHAINP%BEADS(IND+B+1,:)-CHAINP%BEADS(IND+B,:),QCUR)

       ! current twist relative to branch of current piece
       QREL = INVQUAT(QPT(PIECE(B)))*QCUR


       CALL QUAT2EULER(QREL,EULREL)
       CURTWIST = ANGLE2PI(EULREL(1)+EULREL(3))             

       ! how much to rotate
       ANG = (PATHCRD - CURBRANCH)*TWISTPERLEN - CURTWIST

       QROT = ROTQUAT(ANG,(/0D0,0D0,1D0/))
       CHAINP%BEADQ(IND+B) = QCUR*QROT
    ENDDO    
    
    CHAINP%BEADQ(CHAINP%BENDIND(2)) = CHAINP%BENDS(2)%ORIENT*CHAINP%BENDS(2)%TM

  END SUBROUTINE PLACEBEADSONPATH

  SUBROUTINE INTERPOLATEBEADS(CHAINP)
    ! given a chain where all bends have already been placed
    ! (eg: from a previously output configuration or according to helix coordinates)
    ! fill in all linker and tail beads
    ! this fills in the BEADS and BEADQ arrays 
    USE KEYS, ONLY : BPLEN, RANDSTART, RSRANGE
    USE CHAINUTILS, ONLY : CHAIN, BEND, PLACETAILBEADS
    USE MT19937, ONLY : GRND ! random number generator
    USE QUATUTILS ! quaternion utilities
    USE GENUTILS, ONLY : NORM, NORMALIZE, CROSS_PRODUCT, RANDOMAXIS
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: B, I, IND, PLACEBEADS
    TYPE(BEND), POINTER :: BP, BP1
    TYPE(QUATERNION) :: QP, QREL, QPREV, QZ, QWANT
    DOUBLE PRECISION :: PP(3), PM(3), VEC(3), ROTMAT(3,3), EULREL(3)
    DOUBLE PRECISION :: DIST, PHI, R1, R2, R3, TWISTERR, AX(3), RANDAX(3)

    DO B = 1,CHAINP%NBEND-1
       BP => CHAINP%BENDS(B)
       IND = CHAINP%BENDIND(B)+1
       PLACEBEADS = CHAINP%BENDIND(B+1)-IND+1 ! number of beads to place

       CALL QUAT2EULER(BP%ORIENT,EULREL)

       ! place the edge beads
       QP = PTQUAT(BP%POSP)
       QP = BP%ORIENT*QP/BP%ORIENT
       PP = BP%CENT + (/QP%X,QP%Y,QP%Z/)
       ! shift by half a basepair
       PP = PP + QUAT2PT(BP%ORIENT*BP%TP*PTQUAT((/0D0,0D0,BPLEN/2/))/(BP%ORIENT*BP%TP))

       BP1 => CHAINP%BENDS(B+1)
       QP = PTQUAT(BP1%POSM)
       QP = BP1%ORIENT*QP/BP1%ORIENT
       PM = BP1%CENT + (/QP%X,QP%Y,QP%Z/)
       PM = PM - QUAT2PT(BP1%ORIENT*BP1%TM*PTQUAT((/0D0,0D0,BPLEN/2/))/(BP1%ORIENT*BP1%TM))      

       CHAINP%BEADS(IND,:) = PP

       ! vector from edge bead to edge bead
       VEC = PM-PP
       
       qprev = bp%orient*bp%tp       

       ! Place all the beads in a straight line from one edge to the other
       DO I = IND+1,CHAINP%BENDIND(B+1)
          IF (I.EQ.CHAINP%BENDIND(B+1)) THEN
             CHAINP%BEADS(I,:) = PM
          ELSE
             CHAINP%BEADS(I,:) = CHAINP%BEADS(I-1,:) + VEC*1D0/(CHAINP%BENDIND(B+1)-IND)
          ENDIF
          IF (RANDSTART.AND.I.NE.CHAINP%BENDIND(B+1)) THEN
             ! random perturbation of each bead position
             R1 = GRND()*RSRANGE/3; R2 = GRND()*RSRANGE/3; R3 = GRND()*RSRANGE/3
             CHAINP%BEADS(I,:) = CHAINP%BEADS(I,:) + (/R1,R2,R3/)
          ENDIF
       ENDDO

       DIST = NORM(PM-CHAINP%BEADS(CHAINP%BENDIND(B+1),:))
       IF (DIST.GT.1D-10) THEN
          PRINT*, 'ERROR IN INTERPOLATEBEADS: end bead does not match up with bottom edge of bend', dist
          STOP 1
       ENDIF

       ! set the bead orientations
       QPREV = BP%ORIENT*BP%TP       
       DO I = IND,CHAINP%BENDIND(B+1)-1                      
          ! start with 0 alpha and gamma angle
          CALL COORDS2QUAT(0D0, CHAINP%BEADS(I+1,:)-CHAINP%BEADS(I,:),QP)
          ! relative quaternion (relative to previous orientation
          IF (I.EQ.IND) THEN
             QREL = INVQUAT(QPREV)*QP
          ELSE
             QREL = INVQUAT(QPREV)*QP
          ENDIF
          ! relative euler angles
          CALL QUAT2EULER(QREL,EULREL)

          ! relative quaternion with gamma angle adjusted to give desired twist
          IF (I.EQ.IND) THEN
             PHI = ANGLE0(EULREL(1)+EULREL(3) - (BPLEN+CHAINP%LS)/2*CHAINP%TWIST)
          ELSE
             PHI = ANGLE0(EULREL(1)+EULREL(3) - CHAINP%LS*CHAINP%TWIST)
          ENDIF
          CALL EULER2QUAT((/EULREL(1),EULREL(2),EULREL(3)-PHI/),QREL)

          ! absolute orientation
          CHAINP%BEADQ(I) = QPREV*QREL                    
          QPREV = CHAINP%BEADQ(I)
       ENDDO
       CHAINP%BEADQ(CHAINP%BENDIND(B+1)) = BP1%ORIENT*BP1%TM

       ! current quaternion of bend edge relative to last segment
       QREL = INVQUAT(CHAINP%BEADQ(CHAINP%BENDIND(B+1)-1))*CHAINP%BEADQ(CHAINP%BENDIND(B+1))
        CALL QUAT2EULER(QREL,EULREL)
       ! desired quaternion of bend edge relative to last segment
       CALL EULER2QUAT((/(BPLEN/2+CHAINP%LS/2)*CHAINP%TWIST-EULREL(3),EULREL(2),EULREL(3)/),QWANT)

       ! how much do we have to rotate the last segment to get the desired orientation?
       QREL = QREL/QWANT
       CALL QUAT2EULER(QREL,EULREL)      
       IF (EULREL(2).NE.0.OR.EULREL(1).NE.0) THEN
          PRINT*, 'WARNING: something is wrong in interpolating beads! Structure will not be ideal'          
       ENDIF

       ! go through and adjust twist to spread the defect out over the chain
        TWISTERR = ANGLE0(EULREL(3))
        !print*, 'testxb:', twisterr, TWISTERR/CHAINP%NSEG
        TWISTERR = TWISTERR/(CHAINP%NSEG+1)

        DO I = IND,CHAINP%BENDIND(B+1)-1 
           CHAINP%BEADQ(I) = CHAINP%BEADQ(I)*ROTQUAT(TWISTERR*(I-IND+1),(/0D0,0D0,1D0/))           
        ENDDO
    ENDDO

    CALL PLACETAILBEADS(CHAINP,RANDSTART,RSRANGE)

  END SUBROUTINE INTERPOLATEBEADS
END MODULE INPUTSTRUCTS
