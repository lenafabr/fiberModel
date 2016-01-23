MODULE HELIXUTILS
  ! utilities for setting up, defining, and dealing with regular helix coordinates
  USE QUATUTILS
  IMPLICIT NONE

CONTAINS
  SUBROUTINE ROTATECHAIN(CHAINP,QROT)
    ! rotate the chain configuration by the given quaternion
    ! This only makes sense for a single-nucleosome calculation (equivalently: nucleosomes are right on top of each other)
    ! NOTE: this acts only on the single vector representation of the chain (CHAINP%VEC)    
    USE CHAINUTILS, ONLY : CHAIN
    USE ENERGYUTILS, ONLY : GETTOPEDGE, GETBOTTOMEDGE
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(QUATERNION), INTENT(IN) :: QROT
    INTEGER :: I
    DOUBLE PRECISION :: P1(3), P2(3), ZVEC(3), EUL(3), VEC(3)
    DOUBLE PRECISION :: PBOT(3), PTOP(3), MATTOP(3,3), MATBOT(3,3)
    TYPE(QUATERNION) :: STARTQ

    IF (CHAINP%VEC(1).NE.0.OR.CHAINP%VEC(2).NE.0) THEN
       PRINT*, 'ERROR IN ROTATECHAIN: cannot rotate the chain for an ordinary &
            & regular fiber configuration. This only works for nucleosomes &
            & right on top of each other.'
       PRINT*, CHAINP%VEC(1:6)
       STOP 1
    ENDIF    

    ! reorient the  segment orientations
    CALL GETTOPEDGE(CHAINP,PTOP,MATTOP)
    CALL GETBOTTOMEDGE(CHAINP,PBOT,MATBOT)
    DO I = 1,CHAINP%NSEG

       IF (I.EQ.1) THEN
          P1 = PTOP
       ELSE
          P1 = CHAINP%VEC(4*I:4*I+2)
       ENDIF
       IF (I.EQ.CHAINP%NSEG) THEN
          P2 = PBOT
       ELSE
          P2 = CHAINP%VEC(4*I+4:4*I+6)
       ENDIF

       ZVEC = P2-P1
       CALL COORDS2QUAT(CHAINP%VEC(4*I+3),ZVEC,STARTQ)
      
       ! convert segment quaternion into euler angles
       CALL QUAT2EULER(QROT*STARTQ,EUL)
       ! reorient segment
       CHAINP%VEC(4*I+3) = EUL(1)+EUL(3)       
    ENDDO

    ! change the nucleosome orientation
    CALL EULER2QUAT(CHAINP%VEC(4:6),STARTQ)
    CALL QUAT2EULER(QROT*STARTQ,CHAINP%VEC(4:6))

    ! move beads
    DO I = 1,CHAINP%NSEG-1
       ! change bead position (rotate relative to nuc center)
       VEC = CHAINP%VEC(4*I+4:4*I+6) - (/CHAINP%VEC(3),0D0,0D0/)
       CHAINP%VEC(4*I+4:4*I+6) = QUAT2PT(QROT*PTQUAT(VEC)/QROT) 
       CHAINP%VEC(4*I+4) = CHAINP%VEC(4*I+4)+CHAINP%VEC(3)
    ENDDO
  END SUBROUTINE ROTATECHAIN


  SUBROUTINE CHECKGIMBAL(CHAINP,GIMBAL,BEADHIT)
    ! check for gimbal lock, which only occurs if some vector points along the negative z axis
    ! and should thus only be a problem for single nucleosome calculations
    ! also check if any 2 beads are bumping into each other
    USE CHAINUTILS, ONLY : CHAIN
    USE ENERGYUTILS, ONLY : GETTOPEDGE, GETBOTTOMEDGE
    USE GENUTILS, ONLY : NORM
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    LOGICAL, INTENT(OUT) :: GIMBAL,BEADHIT
    INTEGER :: I
    DOUBLE PRECISION :: VEC(3), NVEC, PTOP(3), PBOT(3)
    DOUBLE PRECISION :: NV, MATTOP(3,3), MATBOT(3,3)   
    CALL GETTOPEDGE(CHAINP,PTOP,MATTOP)
    CALL GETBOTTOMEDGE(CHAINP,PBOT,MATBOT)
    GIMBAL = .FALSE.; BEADHIT = .FALSE.
    DO I = 1,CHAINP%NSEG
       IF (I.EQ.1) THEN
          VEC = CHAINP%VEC(4*I+4:4*I+6) - PTOP
       ELSEIF (I.EQ.CHAINP%NSEG) THEN
          VEC = PBOT - CHAINP%VEC(4*I:4*I+2)
       ELSE
          VEC = CHAINP%VEC(4*I+4:4*I+6)-CHAINP%VEC(4*I:4*I+2)
       ENDIF
       NV = NORM(VEC)
       IF (NV/CHAINP%LS.LT.1D-3) THEN
          BEADHIT = .TRUE.
          EXIT
       ENDIF
       VEC = VEC/NV
       IF (VEC(1)**2+VEC(2)**2.LT.NZTINY*1000) THEN
          GIMBAL = .TRUE.
          EXIT
       ENDIF
    ENDDO

  END SUBROUTINE CHECKGIMBAL

  SUBROUTINE MAKESTRAIGHTLINKERHELIX(CHAINP,HELCRDRESET,HELCRD)
    ! for a chain where all arrays and parameters have been set up
    ! set the configuration to the straight linker one
    ! optionally, HELCRDRESET is a logical array for which helical coordinates
    ! should be reset to the values in HELCRD

    USE CHAINUTILS, ONLY : CHAIN
    USE INPUTSTRUCTS, ONLY : INTERPOLATEBEADS
    USE KEYS, ONLY : BPLEN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    LOGICAL, INTENT(IN), OPTIONAL :: HELCRDRESET(6)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: HELCRD(6)
    TYPE(QUATERNION) :: QLINK, QROT
    DOUBLE PRECISION :: ZAX(3), PM(3), TRANS(3)
    integer :: i

    IF (.NOT.CHAINP%ARRAYSET) THEN
       PRINT*, 'ERROR IN GETSTRAIGHTLINKERHELIX: chain arrays have not yet been set up.'
       STOP 1
    ENDIF

    ! quaternion associated with twist along the linker
    QLINK = ROTQUAT(CHAINP%TWIST*(CHAINP%LS*CHAINP%NSEG+BPLEN),(/0D0,0D0,1D0/))

    ! get the desired orientation of the 2nd nucleosome relative to the 1st    
    QROT = CHAINP%BENDS(1)%TP*QLINK/CHAINP%BENDS(2)%TM

    ! get the desired translation of the 2nd nucleosome center relative to the 1st
    ZAX = QUAT2PT(CHAINP%BENDS(1)%TP*PTQUAT((/0D0,0D0,1D0/))/CHAINP%BENDS(1)%TP)
    PM = QUAT2PT(QROT*PTQUAT(CHAINP%BENDS(2)%POSM)/QROT)
    TRANS = CHAINP%BENDS(1)%POSP + ZAX*(CHAINP%LS*CHAINP%NSEG +BPLEN)-PM

    ! convert rotation + translation to helical coordinates
    CALL QUAT2SCREW(QROT,TRANS,CHAINP%VEC(1:6))

    IF (PRESENT(HELCRDRESET)) THEN
       IF (ANY(HELCRDRESET).AND..NOT.PRESENT(HELCRD)) THEN
          PRINT*, 'ERROR in MAKESTRAIGHTLINKERHELIX: cannot reset helical coordinates without HELCRD array supplied'
          STOP 1
       ENDIF

       DO I = 1,6
          IF (HELCRDRESET(I)) CHAINP%VEC(I) = HELCRD(I)          
       ENDDO
    ENDIF

    ! place and orient all the bends
    CALL PLACEBENDSFROMHELIX(CHAINP)    

    ! fill in the linker beads
    CALL INTERPOLATEBEADS(CHAINP)
    ! single vector coordinates for linker beads
    CALL GETHELIXREP(CHAINP,.TRUE.)
   
  END SUBROUTINE MAKESTRAIGHTLINKERHELIX

  SUBROUTINE HELIXDIST(CHAINP1,CHAINP2,DIST)
    ! distance metric for comparing two helical configuration
    ! based on the single-vector regular helix coordinates
    USE CHAINUTILS, ONLY : CHAIN
    TYPE(CHAIN), POINTER :: CHAINP1
    TYPE(CHAIN), POINTER :: CHAINP2
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    TYPE(QUATERNION) :: QREL, Q1, Q2
    DOUBLE PRECISION :: THETA
    DOUBLE PRECISION :: POS1(3), POS2(3), DIFF(3)
    INTEGER :: B

    IF (CHAINP1%NSEG.NE.CHAINP2%NSEG) THEN
       PRINT*, 'ERROR IN HELIXDIST: chains have different number of segments', CHAINP1%NSEG, CHAINP2%NSEG
       STOP 1
    ENDIF

    ! get the relative orientation of the nucleomes
    CALL EULER2QUAT(CHAINP1%VEC(4:6),Q1)
    CALL EULER2QUAT(CHAINP2%VEC(4:6),Q2)
    QREL = INVQUAT(Q1)*Q2
    ! angle of rotation
    THETA = ACOS(QREL%W)*2

    ! position of 2nd nucleosome
    POS1 = (/CHAINP1%VEC(3)*COS(CHAINP1%VEC(2)), &
          & CHAINP1%VEC(3)*SIN(CHAINP1%VEC(2)),CHAINP1%VEC(1)/)
    POS2 = (/CHAINP2%VEC(3)*COS(CHAINP2%VEC(2)), &
         & CHAINP2%VEC(3)*SIN(CHAINP2%VEC(2)),CHAINP2%VEC(1)/)
    DIFF = POS2 - POS1

    DIST = DOT_PRODUCT(DIFF,DIFF)/CHAINP1%LS**2+ THETA**2

    ! positions of linker beads
    DO B = 1,CHAINP1%NSEG
       DIFF = CHAINP1%VEC(4*B+4:4*B+6) - CHAINP2%VEC(4*B+4:4*B+6)
       DIST = DIST + DOT_PRODUCT(DIFF,DIFF)/CHAINP1%LS**2
    ENDDO

    DIST = SQRT(DIST)
  END SUBROUTINE HELIXDIST

  SUBROUTINE OUTPUTREPLIC(CHAINP,NBEND,OUTFILE,REPLICP)
    ! Replicate the regular chain structure to have a given number of bends
    ! Dump the result in OUTFILE
    ! if REPLICP is supplied, return a pointer to the replicated chain object
    ! otherwise, clean up the chain object at end of subroutine
    USE CHAINUTILS, ONLY : CHAIN, CREATECHAIN, OUTPUTCHAIN, CLEANUPCHAIN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: NBEND
    CHARACTER (LEN=*) :: OUTFILE
    TYPE(CHAIN), POINTER, OPTIONAL :: REPLICP
    TYPE(CHAIN), POINTER :: CP

    PRINT*, 'Replicating chain to ', NBEND, ' bends and outputting into ', TRIM(OUTFILE)
    IF (PRESENT(REPLICP)) THEN
       ALLOCATE(REPLICP)
       CALL CREATECHAIN(REPLICP,NBEND)
       REPLICP%VEC = CHAINP%VEC
       CALL FROMHELIXREP(REPLICP)
       CALL OUTPUTCHAIN(REPLICP,OUTFILE)
    ELSE
       ALLOCATE(CP)
       CALL CREATECHAIN(CP,NBEND)
       CP%VEC = CHAINP%VEC
       CALL FROMHELIXREP(CP)
       CALL OUTPUTCHAIN(CP,OUTFILE)

       CALL CLEANUPCHAIN(CP)
       DEALLOCATE(CP)
    ENDIF


  END SUBROUTINE OUTPUTREPLIC

  SUBROUTINE GETNUCDISTS(CHAINP,DISTS)
    ! find the distances between the 0th nucleosome and each subsequent one
    ! approximate distances  between cylindrical shells 
    ! (specifically: the distance between the nucleosome centers minus
    ! the minimum distance they can approach without touching)
    USE CHAINUTILS, ONLY : CHAIN
    USE TRICUBSPLINEUTILS, ONLY : EVALUATESPLINE
    USE CYLINDERUTILS, ONLY : NUCCYL, NUCCYLSET
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: DISTS(CHAINP%NINTERACT)
    INTEGER :: NNUC
    DOUBLE PRECISION :: SNT, CNT, SB, CB, SA, CA, NTA, SNTA, CNTA, H, R
    DOUBLE PRECISION :: IJ(3), KJ(3), LK(3)
    DOUBLE PRECISION :: MINDIST, ACTDIST, CT1, CT2, PHI
    
    IF (.NOT.NUCCYLSET) THEN
       PRINT*, 'ERROR IN GETNUCDISTS: cannot get nucleosome distances as cylinder array is not set.'
       DISTS = 0D0
       RETURN
    ENDIF
       
    
    DO NNUC = 1,CHAINP%NINTERACT       
       SNT = SIN(NNUC*CHAINP%VEC(2)); CNT = COS(NNUC*CHAINP%VEC(2))
       SB = SIN(CHAINP%VEC(5)); CB = COS(CHAINP%VEC(5))
       SA = SIN(CHAINP%VEC(4)); CA = COS(CHAINP%VEC(4))
       NTA = NNUC*CHAINP%VEC(2) + CHAINP%VEC(4); SNTA = SIN(NTA); CNTA = COS(NTA)
       H = CHAINP%VEC(1); R = CHAINP%VEC(3); 
       
       IJ = (/SB*SA,-SB*CA,CB/)
       KJ = (/R*(CNT-1),R*SNT,H*NNUC/)
       LK = (/SB*SNTA,-SB*CNTA,CB/)

       ! get angles between cylinders and connecting vector
       CALL GETANGLE(IJ,KJ,CT1)
       CALL GETANGLE(-KJ,LK,CT2)
       ! get dihedral angle between cylinders
       CALL GETDIHEDRAL(IJ,-KJ,LK,PHI)   

       ! interpolate the minimal separation distance 
       CALL EVALUATESPLINE(NUCCYL,CT1,CT2,PHI,MINDIST)
       ! actual cylinder separation distance
       ACTDIST = SQRT(DOT_PRODUCT(KJ,KJ))

       ! approximate distance between shells
       DISTS(NNUC) = ACTDIST-MINDIST
    ENDDO
    
  END SUBROUTINE GETNUCDISTS

  SUBROUTINE GETFIBERDIAM(CHAINP)
    ! get the overall diameter of the fiber, using cylindrical-shaped nucleosomes
    ! use pre-tabulated cylindrical sterics with a reference cylinder
    USE CYLINDERUTILS, ONLY : NUCCYL, NUCCYLSET
    USE TRICUBSPLINEUTILS, ONLY : EVALUATESPLINE, SPLINETEST
    USE CYLINDERUTILS, ONLY : GETMINCENTDIST
    USE CHAINUTILS, ONLY : CHAIN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION :: A, B, CT1, CT2, PHI, CA, SA, CB, SB
    DOUBLE PRECISION :: DMDC1,DMDC2, DMDP,MINDIST, tmp

    IF (.NOT.NUCCYLSET) THEN
       PRINT*, 'Cannot find fiber diameter because cylinder arrays are not set.'
       CHAINP%DIAM = 0D0
       RETURN
    ENDIF

    A = CHAINP%VEC(4); B = CHAINP%VEC(5)
    CA = COS(A); SA = SIN(A); CB = COS(B); SB = SIN(B)
    CT2 = 0D0
    CT1 = SB*SA
    PHI = ATAN2(-SB*CA,CB)

    CALL EVALUATESPLINE(NUCCYL,CT1,CT2,PHI,MINDIST,DMDC1,DMDC2,DMDP)

    chainp%DIAM = 2*(CHAINP%VEC(3) + MINDIST-CHAINP%BENDS(1)%RADIUS)
    
  END SUBROUTINE GETFIBERDIAM

  SUBROUTINE FROMHELIXREP(CHAINP)
    ! convert from the single vector representation of the regular helix
    ! (6 helix coordinates followed by absolute alpha+gamma and bead
    ! position for each linker segment)
    ! to an ordinary chain representation   
    USE KEYS, ONLY : BPLEN
    USE CHAINUTILS, ONLY : CHAIN, BEND, PLACETAILBEADS
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: I,B,IND, IND0, IND2, I0
    TYPE(BEND), POINTER :: BP, BP2
    TYPE(QUATERNION) :: QROT,QP
    DOUBLE PRECISION :: SHIFT, ZAX(3)

    IF (.NOT.CHAINP%arraySET) THEN
       PRINT*, 'ERROR IN FROMHELIXREP: arrays have not been allocated'
       STOP 1   
    ENDIF

    ! place and orient all the bends
    CALL PLACEBENDSFROMHELIX(CHAINP)

    ! place the fiber tails (exiting straight from the bend)
    CALL PLACETAILBEADS(CHAINP,.FALSE.)
    
    IF (CHAINP%NBEND.LT.2) RETURN

    ! place the first linker
    IND = CHAINP%BENDIND(1)+1
    BP=>CHAINP%BENDS(1)

    ! first edge bead of linker
    CHAINP%BEADS(IND,:) = BP%CENT + QUAT2PT(BP%ORIENT*PTQUAT(BP%POSP)/BP%ORIENT)

    ! shift by half a bp along tangent
    QP = BP%ORIENT*BP%TP
    CHAINP%BEADS(IND,:) = CHAINP%BEADS(IND,:) + QUAT2PT(QP*PTQUAT((/0D0,0D0,BPLEN/2/))/QP)

    DO I = 1,CHAINP%NSEG-1
       CHAINP%BEADS(IND+I,:) = CHAINP%VEC(4*I+4:4*I+6)

       ! orientation according to following segment
       ZAX = CHAINP%BEADS(IND+I,:) - CHAINP%BEADS(IND+I-1,:)
       CALL COORDS2QUAT(CHAINP%VEC(4*I+3),ZAX,CHAINP%BEADQ(IND+I-1))
    ENDDO

    ! place final bead
    IND2 = CHAINP%BENDIND(2)
    BP2 => CHAINP%BENDS(2)
    CHAINP%BEADS(IND2,:) = BP2%CENT + QUAT2PT(BP2%ORIENT*PTQUAT(BP2%POSM)/BP2%ORIENT)
    ! shift by half a bp along tangent
    QP = BP2%ORIENT*BP2%TM
    CHAINP%BEADS(IND2,:) =  CHAINP%BEADS(IND2,:) - QUAT2PT(QP*PTQUAT((/0D0,0D0,BPLEN/2/))/QP)

    ! orientation of last segment
    ZAX = CHAINP%BEADS(IND2,:) - CHAINP%BEADS(IND2-1,:)
    CALL COORDS2QUAT(CHAINP%VEC(4*CHAINP%NSEG+3),ZAX,CHAINP%BEADQ(IND2-1))
    ! orientation at last bead matches bottom end of bend
    CHAINP%BEADQ(IND2) = QP

    IND0 = CHAINP%BENDIND(1)+1
    DO B = 2,CHAINP%NBEND-1
       ! rotate and translate everything along the helix
       IND = CHAINP%BENDIND(B)+1
       IND2 = CHAINP%BENDIND(B+1)      

       ! rotation around fiber axis
       QROT = ROTQUAT(CHAINP%VEC(2)*(B-1),(/0D0,0D0,1D0/))
       ! shift along fiber axis
       SHIFT = CHAINP%VEC(1)*(B-1)

       DO I = IND,IND2
          I0 = I-IND+IND0
          CHAINP%BEADQ(I) = QROT*CHAINP%BEADQ(I0)

          CHAINP%BEADS(I,:) = QUAT2PT(QROT*PTQUAT(CHAINP%BEADS(I0,:))/QROT)
          CHAINP%BEADS(I,3) = CHAINP%BEADS(I,3)+SHIFT
       ENDDO
    ENDDO

  END SUBROUTINE FROMHELIXREP

  SUBROUTINE MAKEHELIX(CHAINP,HPERNUC,THETA,RAD,EUL)
    ! create a helical structure with the given parameters
    ! parameters are:
    ! height per nucleosome along fiber axis, 
    ! angle per nucleosome around the axis
    ! radius of the helix
    ! euler angles defining 1st nucleosome orienation relative to fiber
    ! WARNING: CHAINP must already have all its arrays allocated and 
    ! parameters set (eg: with CREATECHAIN subroutine); 
    ! this subroutine only sets actual coordinates
    USE CHAINUTILS, ONLY : CHAIN
    USE INPUTSTRUCTS, ONLY : INTERPOLATEBEADS

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(IN) :: HPERNUC, THETA, RAD, EUL(3)

    IF (.NOT.CHAINP%ARRAYSET) THEN
       PRINT*, 'ERROR IN MAKEHELIX: arrays for chain must already be set before putting in helix coordinates. Run CREATECHAIN first'
       STOP 1
    ENDIF

    ! set up the regular helix coordinates used for minimization / searching 
    ! configuration space

    CHAINP%VEC(1:6) = (/HPERNUC,THETA,RAD,EUL(1),EUL(2),EUL(3)/)       

    ! place the bends
    CALL PLACEBENDSFROMHELIX(CHAINP)

    ! place the beads
    CALL INTERPOLATEBEADS(CHAINP)
    
  END SUBROUTINE MAKEHELIX

  SUBROUTINE PLACEBENDSFROMHELIX(CHAINP)
    ! place and orient all the bends in the chain according to the regular
    ! helix coordinates stored in CHAINP%VEC(1:6)
    USE CHAINUTILS, ONLY : CHAIN, BEND

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    TYPE(BEND), POINTER :: BP
    DOUBLE PRECISION :: H, T, R
    INTEGER :: B

    H = CHAINP%VEC(1); T = CHAINP%VEC(2); R = CHAINP%VEC(3)    

    DO B = 1,CHAINP%NBEND
       BP=>CHAINP%BENDS(B)

       BP%CENT = (/R*COS((B-1)*T),R*SIN((B-1)*T),H*(B-1)/)
       CALL EULER2QUAT((/(B-1)*CHAINP%VEC(2)+CHAINP%VEC(4),CHAINP%VEC(5),CHAINP%VEC(6)/),BP%ORIENT)
    ENDDO
  END SUBROUTINE PLACEBENDSFROMHELIX

  SUBROUTINE GETHELIXREP(CHAINP, LINKERONLY)
    ! convert the BEADS, BEADQ, BEND%CENT, BEND%ORIENT general 
    ! representation of the chain that is used for output and input 
    ! into the regular helix coordinates used for minimization 
    ! (6 coordinates for helix geometry, 4 coordinates per linker segment 
    ! giving alpha+gamma euler angle and bead position)
    ! if LINKERONLY is supplied and true, then assume the first 6 coordinates
    ! of the regular helix representation are already set and only 
    ! set the linker DNA coordinates
    ! Only the first 2 bends of the chain are used for getting 
    ! the helix coordinates

    USE CHAINUTILS, ONLY : CHAIN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    LOGICAL, INTENT(IN), OPTIONAL :: LINKERONLY
    LOGICAL :: SETHELCRD
    DOUBLE PRECISION :: TRANS(3), HELCRD(6), EUL(3)
    TYPE(QUATERNION) :: QROT
    INTEGER :: I, IND

    IF (CHAINP%NBEND.LT.2.OR..NOT.CHAINP%ARRAYSET) THEN
       PRINT*, 'ERROR IN GETHELIXREP: chainp either has less than 2 bends or has not had its arrays allocated'
       PRINT*, CHAINP%NBEND, CHAINP%ARRAYSET
       STOP 1
    ENDIF

    IF (PRESENT(LINKERONLY)) THEN
       SETHELCRD = .NOT.LINKERONLY
    ELSE
       SETHELCRD = .TRUE.
    ENDIF

    IF (SETHELCRD) THEN
       ! get the regular helix coordinates
       ! based on the relative orientation and position of the first 2 bends

       ! Relative translation and orientation
       TRANS = CHAINP%BENDS(2)%CENT-CHAINP%BENDS(1)%CENT
       ! convert translation to nucleosome coordinate system
       TRANS = QUAT2PT(INVQUAT(CHAINP%BENDS(1)%ORIENT)*PTQUAT(TRANS)*CHAINP%BENDS(1)%ORIENT)
       QROT = INVQUAT(CHAINP%BENDS(1)%ORIENT)*CHAINP%BENDS(2)%ORIENT

       ! convert to helical coordinates
       CALL QUAT2SCREW(QROT, TRANS, HELCRD)
       
       ! regularize the helix coordinates (no negative heights, etc)
       CALL REGULARIZEHELCOORDS(HELCRD)

       CHAINP%VEC(1:6) = HELCRD
    ENDIF

    ! Set the linker DNA coordinates
    IND = CHAINP%BENDIND(1)+1
    DO I = 1,CHAINP%NSEG
       ! convert quaternion for the segment into euler angles
       CALL QUAT2EULER(CHAINP%BEADQ(IND+I-1),EUL)
       CHAINP%VEC(4*I+3) = EUL(1)+EUL(3)
       ! bead position
       IF (I.LT. CHAINP%NSEG) THEN
          CHAINP%VEC(4*I+4:4*I+6) = CHAINP%BEADS(IND+I,:)
       ENDIF
    ENDDO
  END SUBROUTINE GETHELIXREP

  SUBROUTINE REGULARIZEHELCOORDS(VEC)
    ! regularize the various euler angles for the regular helix coords
    ! to make them fall within the ordinary range of euler angles
    ! also get rid of negative radii and negative heights
    DOUBLE PRECISION, intent(inout) :: VEC(6)
    DOUBLE PRECISION :: ALPHA, BETA, GAMMA

    ! make sure height is positive
    IF (VEC(1).LT.0) THEN       
       VEC(1) = ABS(VEC(1))
       VEC(2) = -VEC(2)
       VEC(5) = PI-VEC(5)
       VEC(4) = PI-VEC(4)
       VEC(6) = PI + VEC(6)
    ENDIF

    ! make sure radius is positive
    IF (VEC(3).LT.0) THEN       
       VEC(3) = -VEC(3)
       VEC(4) = VEC(4) + PI
    ENDIF

    ! make sure theta angle is between 0 and 2pi
    VEC(2) = ANGLE2PI(VEC(2))

    ! fix up the euler angles
    ALPHA = VEC(4); BETA = VEC(5); GAMMA = VEC(6)
    IF (BETA.GT.PI) THEN
       BETA = 2*PI-BETA
       ALPHA = ALPHA + PI
       GAMMA = GAMMA + PI
    ELSEIF (BETA.LT.0) THEN
       BETA = -BETA
       ALPHA = ALPHA + PI
       GAMMA = GAMMA + PI
    ENDIF
    ALPHA = ANGLE2PI(ALPHA)
    GAMMA = ANGLE2PI(GAMMA)
    
    VEC(4:6) = (/ALPHA,BETA,GAMMA/)   
  END SUBROUTINE REGULARIZEHELCOORDS

  ! ---------- for debugging only ----------
  ! input and output of helix coords directly
  SUBROUTINE DUMPHELIXCRD(CHAINP,FILENAME)    
    ! dump regular helix coordinats into a file
    USE CHAINUTILS, ONLY : CHAIN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    CHARACTER(LEN=*) :: FILENAME
    INTEGER :: I

    OPEN(UNIT=77,FILE=FILENAME)
    DO I = 1,CHAINP%NCRD
       WRITE(77,'(G25.15)') CHAINP%VEC(I)
    ENDDO
    CLOSE(77)
  END SUBROUTINE DUMPHELIXCRD

  SUBROUTINE HELCRDFROMFILE(CHAINP,FILENAME)
    ! for a chain that is already set up, replace the coordinates with those from a file    
    USE CHAINUTILS, ONLY : CHAIN
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    CHARACTER(LEN=*) :: FILENAME
    INTEGER :: I

    OPEN(UNIT=77,FILE=FILENAME,STATUS='OLD')
    DO I = 1,CHAINP%NCRD
       READ(77,'(G25.15)') CHAINP%VEC(I)
    ENDDO
    CLOSE(77)

    CALL FROMHELIXREP(CHAINP)
  END SUBROUTINE HELCRDFROMFILE

END MODULE HELIXUTILS
