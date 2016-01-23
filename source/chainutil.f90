MODULE CHAINUTILS
  ! utilities for defining chain objects
  ! and basic manipulations of chain configurations

  USE QUATUTILS
  USE GENUTILS

  IMPLICIT NONE
  LOGICAL :: CHAINTEST

  TYPE BEND
     ! a bend in the DNA treated as a rigid body
     ! specifically each one of these represents a nucleosome

     ! ----------
     ! structural parameters that define the nature of the bend
     ! ----------
     ! orientations of +/- tangents relative to the coordinate system of the overall bend
     TYPE(QUATERNION) :: TP,TM
     ! position of +/- endpoints, relative to bend center
     DOUBLE PRECISION :: POSP(3), POSM(3)
     ! parameters defining shape for steric repulsion
     ! each bend object is assumed cylindrical in shape     
     DOUBLE PRECISION :: RADIUS,HEIGHT,MAXRAD
     ! filler DNA beads that are held rigid with the bend
     ! these are only used for visual output of the configuration
     INTEGER :: NFILL
     DOUBLE PRECISION, POINTER :: FILLBEADS(:,:)
     TYPE(QUATERNION), POINTER :: FILLQ(:)
    
     ! -------------- 
     ! information for DisCO internucleosomal potential
     ! --------------
     ! charge positions and magnitudes 
     DOUBLE PRECISION, pointer :: CHARGEPOS(:,:), CHARGEMAG(:)
     INTEGER :: NCHARGE
     ! flexible tails
     ! FLEXTAILS: array of coordinates for the tail beads
     ! for each tail, lists: 3 coordinates for each bead relative to 
     ! the bend coordinate system
     DOUBLE PRECISION, POINTER :: FLEXTAILS(:)
     INTEGER :: NTAIL ! number of tails
     INTEGER, POINTER :: TAILnBEAD(:) ! number of beads in each tail
     INTEGER :: TOTTAILBEAD ! total number of tail beads
     ! equilibrium length and stretch force constant for each tail segment
     ! i-th value is for the segment between the i-th and i+1st bead
     DOUBLE PRECISION, POINTER :: TAILL0(:), TAILKL(:)
     ! same thing for the equilibrium angle and bend force
     DOUBLE PRECISION, pointer :: TAILANG0(:), TAILKA(:)
     ! attachment position for each tail
     DOUBLE PRECISION, POINTER :: TAILPOS(:,:)
     ! charges for all the tails, listed consecutively
     DOUBLE PRECISION, POINTER :: TAILCHARGES(:)       

     ! -----------
     ! current data for the bend (changes during simulation)
     ! -----------
     ! center position
     DOUBLE PRECISION :: CENT(3)
     ! orientation of the entire bend, relative to canonical axes
     TYPE(QUATERNION) :: ORIENT

     LOGICAL :: ARRAYSET ! have all the arrays been set for this bend?

     ! ------------
     ! extra marked points to output
     ! ------------
     INTEGER :: NMARKPOINT
     DOUBLE PRECISION, POINTER :: MARKPOINTS(:,:)
     CHARACTER*3, POINTER :: MARKPOINTNAME(:)

  END type BEND

  TYPE CHAIN
     ! object defining a regular (periodic) helical chain
     
     ! ------------
     ! basic geometric properties
     ! -----------     
     ! list of bend objects, one for each nucleosome
     TYPE(BEND), POINTER :: BENDS(:)
     ! list of nucleosomal positions
     ! nucleosome #i is located between linker segment i-1 and segment i
     INTEGER, POINTER :: BENDIND(:)
     ! segment length for linker DNA
     DOUBLE PRECISION :: LS

     ! NSEG: number of linker segments between nucleosomes
     ! NCRD: number of coordinates in the single vector representation of the chain
     ! NCRDDNA: number of coordinates for the DNA and nucleosomes 
     ! (as opposed to flexible tails if using DisCO potential)
     INTEGER :: NSEG, NCRD, NCRDDNA

     ! number of nucleosomes and number of linker DNA "beads"
     ! for the full chain representation 
     ! (not for the periodic chain used for the optimization)
     INTEGER :: NBEND, NBEAD

     ! ------------
     ! energetic parameters for the DNA
     ! ------------
     ! LP = bending persistence length of the chain
     ! LTW = twist persistence length
     ! LSTRETCH = stretch persistence length; 
     ! TWIST = natural twist of the DNA
     ! TENSION = externally applied tension
     DOUBLE PRECISION :: LP, LTW,LSTRETCH, TWIST, TENSION
     ! magnitude of steric energy
     DOUBLE PRECISION :: ESTERIC
     ! turn on sterics between bends? between DNA segments? between the two?
     LOGICAL :: BENDSTERICS,SEGSTERICS,BENDSEGSTERICS
     ! steric radius of DNA
     DOUBLE PRECISION :: DNARAD
     ! how many bends to go out for interactions
     INTEGER :: NINTERACT

     ! -------------
     ! current data about the chain (changes during simulation)
     ! -------------
     ! list of the current positions and orientations
     ! of all the beads defining the linker DNA
     ! used for input and output of chain configuration, NOT for minimization
     DOUBLE PRECISION, POINTER :: BEADS(:,:)
     TYPE(QUATERNION), POINTER :: BEADQ(:)
     ! single vector representation of chain configuration
     ! this lists, in order:
     ! 6 superhelix parameters 
     ! (height per nucleosome, angle per nucleosome, radius of fiber, 
     ! three euler angles for nucleosome orientation relative to fiber)
     ! for each linker segment, the alpha+gamma euler angle followed by 
     ! coordinates for the position of the bead representing the endpoint 
     ! of this segment (only the euler angle is given for the last 
     ! linker segment)     
     DOUBLE PRECISION, POINTER :: VEC(:)
     ! current energy and gradient
     ! ENERGYPARTS contains bend, twist, stretch, nuc-nuc, seg-seg, nuc-seg sterics if no DISCO or LANGOWSKI potentials are used
     DOUBLE PRECISION :: ENERGY, ENERGYPARTS(6)
     DOUBLE PRECISION, POINTER :: GRAD(:)
     ! overall diameter of the chain 
     ! (out to the furthest point of the nucleosomes from the axis)
     DOUBLE PRECISION :: DIAM
     
     ! -------------
     ! parameters for internucleosome interactions
     ! ------------
     ! for DiSCO lennard jones interactions (see Arya, BPS 2006)
     ! excluded volume distance for tail-tail, tail-core, core-core,
     ! tail-linker, and core-linker interactions
     DOUBLE PRECISION :: SIGTT,SIGTC,SIGCC,SIGTL,SIGCL
     ! excluded volume strengths for tail-tail (EPSEVT) and everything else (EPSEV)
     DOUBLE PRECISION :: EPSEV,EPSEVT
     ! inverse debye length
     DOUBLE PRECISION :: LDEBYEINV

      ! parameters for anisotropic LJ potential (Wedemann & Langowski 2004)
     DOUBLE PRECISION :: LJX, LJXP ! anisotropies in potential width and depth
     DOUBLE PRECISION :: EPS0, SIG0 ! potential depth and width scales   
     
     ! -----------
     ! bookkeeping
     ! -----------
     LOGICAL :: ARRAYSET=.false. ! have the various arrays been allocated?

  END type CHAIN


CONTAINS

  SUBROUTINE SAMECONFIG(CHAIN1,CHAIN2,ECUT,DISTCUT,SAME)
    ! check if two chain configurations are the same
    ! two chains are considered the same if all of the following holds:
    ! (1) the energies are within ECUT of each other
    ! (2) With the 0th nucleosomes aligned, the distance between the 1st nucs on the two chains is greater than rmscut
    ! (3) the angle of rotation to get from 0th to 1st nucleosome orientation does not differ by more than RMSCUT

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAIN1, CHAIN2
    DOUBLE PRECISION, INTENT(IN) :: ECUT, DISTCUT
    LOGICAL, INTENT(OUT) :: SAME
    DOUBLE PRECISION :: POS1(3), POS2(3), ANG, DIST,DIFF(3)
    TYPE(QUATERNION) :: Q1,Q2,QREL
    INTEGER :: B 
    LOGICAL :: SINGLENUC

    SAME = .TRUE.

    SINGLENUC = CHAIN1%VEC(1).EQ.0 .AND. CHAIN1%VEC(2).EQ.0 & 
         & .AND. CHAIN2%VEC(1).EQ.0 .AND. CHAIN2%VEC(2).EQ.0

    ! check energies
    IF (ABS(CHAIN1%ENERGY-CHAIN2%ENERGY).GT.ECUT) THEN
       SAME = .FALSE.
    ELSE
       ! check nucleosome distances
       POS1 = (/CHAIN1%VEC(3)*COS(CHAIN1%VEC(2)), &
            & CHAIN1%VEC(3)*SIN(CHAIN1%VEC(2)),CHAIN1%VEC(1)/)
       POS2 = (/CHAIN2%VEC(3)*COS(CHAIN2%VEC(2)), &
            & CHAIN2%VEC(3)*SIN(CHAIN2%VEC(2)),CHAIN2%VEC(1)/)
       DIST = NORM(POS1-POS2)

       !print*, 'testx1:', chain1%energy, chain2%energy, dist, distcut
       
       IF (DIST.GT.DISTCUT*CHAIN1%LS) THEN
          SAME = .FALSE.
       ELSEIF (.NOT.SINGLENUC) THEN
          ! check nucleosome orientations (unless working with single nucleosome)
          CALL EULER2QUAT(CHAIN1%VEC(4:6),Q1)
          CALL EULER2QUAT(CHAIN2%VEC(4:6),Q2)
          QREL = INVQUAT(Q1)*Q2
          ! overall angle to rotate from one to the other
          ANG = 2*ACOS(QREL%W)
          ANG = ANGLE0(ANG)
          !print*, 'testx2:', chain1%energy, chain2%energy, ang, distcut

          SAME = (abs(ANG).LT.DISTCUT)
       ENDIF

       ! also check positions of linker beads       
       IF (SAME) THEN
          !IF (SINGLENUC) THEN
             ! look at positions relative to nucleosome orientation
             DIST = 0D0
             ! orientation of nucleosomes
             CALL EULER2QUAT(CHAIN1%VEC(4:6), Q1)
             Q1 = INVQUAT(Q1)
             CALL EULER2QUAT(CHAIN2%VEC(4:6), Q2)
             Q2 = INVQUAT(Q2)

             DO B = 1,CHAIN1%NSEG-1
                POS1 = CHAIN1%VEC(4*B+4:4*B+6) - (/CHAIN1%VEC(3),0D0,0D0/)
                POS1 = QUAT2PT(Q1*PTQUAT(POS1)/Q1)
                POS2 = CHAIN2%VEC(4*B+4:4*B+6) - (/CHAIN2%VEC(3),0D0,0D0/)
                POS2 = QUAT2PT(Q2*PTQUAT(POS2)/Q2)

                DIFF = POS1 - POS2
                DIST = DIST + DOT_PRODUCT(DIFF,DIFF)/CHAIN1%LS**2         
                !print*, 'testx3:', sqrt(DOT_PRODUCT(DIFF,DIFF)/CHAIN1%LS**2)
             ENDDO
             DIST = SQRT(DIST/(CHAIN1%NSEG-1))
             SAME = DIST.LT.DISTCUT
          !ELSE
             ! look at absolute positions
          !   DIST = 0D0
          !   DO B = 1,CHAIN1%NSEG-1
          !      DIFF = CHAIN1%VEC(4*B+4:4*B+6) - CHAIN2%VEC(4*B+4:4*B+6)
          !      DIST = DIST +DOT_PRODUCT(DIFF,DIFF)/CHAIN1%LS**2
          !     print*, 'testx3:', b, sqrt(DOT_PRODUCT(DIFF,DIFF)/CHAIN1%LS**2)
           !  ENDDO
           !  DIST = SQRT(DIST/(CHAIN1%NSEG-1))
           !  print*, 'testx4:', dist
            ! SAME = DIST.LT.DISTCUT
          !ENDIF
       ENDIF


    ENDIF
  END SUBROUTINE SAMECONFIG

  SUBROUTINE PLACETAILBEADS(CHAINP,RANDOMIZE,RANDRANGE)
    ! Place straight tails at either end of the fiber structure
    ! used for input/output 
    ! (not relevant for the regular helix coords used in the actual minimization)
    ! if RANDOMIZE is present and true then randomize the tail beads a little (not perfectly straight)
    ! RANDRANGE sets how much they are randomized

    USE KEYS, ONLY : BPLEN

    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    LOGICAL, OPTIONAL :: RANDOMIZE
    DOUBLE PRECISION, OPTIONAL :: RANDRANGE
    TYPE(QUATERNION) :: QP, QZ
    TYPE(BEND), POINTER :: BP
    INTEGER :: I
    DOUBLE PRECISION :: AX(3), RANDAX(3)

    IF (PRESENT(RANDOMIZE)) THEN
       IF (RANDOMIZE.AND..NOT.PRESENT(RANDRANGE)) THEN
          PRINT*, 'ERROR IN PLACETAILBEADS: if randomizing tail beads, need to supply the randomization range'
          STOP 1
       ENDIF
    ENDIF

    ! place the beads in the first tail (keeping the tail straight and properly twisted)
    BP => CHAINP%BENDS(1)
    QZ = (BP%ORIENT*BP%TM)*PTQUAT((/0D0,0D0,1D0/))/(BP%ORIENT*BP%TM)
    AX = QUAT2PT(QZ)

    QP = BP%ORIENT*PTQUAT(BP%POSM)/BP%ORIENT
    CHAINP%BEADS(CHAINP%BENDIND(1),:) = BP%CENT+(/QP%X,QP%Y,QP%Z/)-AX*BPLEN/2
    CHAINP%BEADQ(CHAINP%BENDIND(1)) = BP%ORIENT*BP%TM

    DO I = CHAINP%BENDIND(1)-1,1,- 1       
       IF (PRESENT(RANDOMIZE)) THEN
          IF (RANDOMIZE) THEN
             CALL RANDOMAXIS(AX,RANDRANGE,RANDAX)
             AX = RANDAX      
          ENDIF
       ENDIF
       CHAINP%BEADS(I,:) = CHAINP%BEADS(I+1,:) -CHAINP%LS*AX
       ! bead orientation, properly twisted
       IF (I.EQ.CHAINP%BENDIND(1)-1) THEN
          CHAINP%BEADQ(I) = CHAINP%BEADQ(I+1)*ROTQUAT(-CHAINP%TWIST*(CHAINP%LS/2+BPLEN/2),(/0D0,0D0,1D0/))
       ELSE
          CHAINP%BEADQ(I) = CHAINP%BEADQ(I+1)*ROTQUAT(-CHAINP%TWIST*CHAINP%LS,(/0D0,0D0,1D0/))
       ENDIF
    ENDDO

    ! place the beads in the second tail (keeping the tail straight and properly twisted)
    BP => CHAINP%BENDS(CHAINP%NBEND)
    QZ = PTQUAT((/0D0,0D0,1D0/))
    QZ = (BP%ORIENT*BP%TP)*QZ/(BP%ORIENT*BP%TP)
    AX = QUAT2PT(QZ)

    QP = BP%ORIENT*PTQUAT(BP%POSP)/BP%ORIENT
    CHAINP%BEADS(CHAINP%BENDIND(CHAINP%NBEND)+1,:) = BP%CENT+(/QP%X,QP%Y,QP%Z/)+AX*BPLEN/2
    CHAINP%BEADQ(CHAINP%BENDIND(CHAINP%NBEND)+1) = BP%ORIENT*BP%TP*ROTQUAT(CHAINP%TWIST*(CHAINP%LS/2+BPLEN/2),(/0D0,0D0,1D0/))
  
    DO I = CHAINP%BENDIND(CHAINP%NBEND)+2,CHAINP%NBEAD       
       IF (PRESENT(RANDOMIZE)) THEN
          IF (RANDOMIZE) THEN
             CALL RANDOMAXIS(AX,RANDRANGE,RANDAX)
             AX = RANDAX
          ENDIF
       ENDIF

       CHAINP%BEADS(I,:) = CHAINP%BEADS(I-1,:) + CHAINP%LS*AX       
       ! bead orientation, properly twisted       
       CHAINP%BEADQ(I) = CHAINP%BEADQ(I-1)*ROTQUAT(CHAINP%TWIST*CHAINP%LS,(/0D0,0D0,1D0/))
    ENDDO
  END SUBROUTINE PLACETAILBEADS

 SUBROUTINE OUTPUTCHAIN(CHAINP,OUTFILE)
   ! output the structure of a chain for visualizing
   use keyS, ONLY : FILLBEADLIST
   IMPLICIT NONE
   TYPE(CHAIN), POINTER :: CHAINP
   CHARACTER(LEN=*) :: OUTFILE
   TYPE(BEND), POINTER :: BP
   INTEGER :: I, J, IND, COUNT, TMPI, MAXIND, B
   TYPE(QUATERNION) :: QAX
   DOUBLE PRECISION :: BEAD(3), BRANCH(3), POS(3)

!   print*, 'Outputting chain structure to file: ', OUTFILE

   OPEN(UNIT=99,FILE=OUTFILE,STATUS='UNKNOWN')
   DO J = 0,CHAINP%NBEND
      IF (J.EQ.0) THEN
         IND = 1
      ELSE
         IND = CHAINP%BENDIND(J)+1

         ! output the bound "filler beads" on the J-th bend
         BP=>CHAINP%BENDS(J)
         IF(BP%NFILL.EQ.0) THEN
            ! put in a bead half-way between edges          
            ! with an arbitrary twist
            CALL CROSS_PRODUCT(CHAINP%BEADS(IND,:)-CHAINP%BEADS(IND-1,:),(/1D0,0D0,0D0/), BRANCH)
            CALL NORMALIZE(BRANCH)
            WRITE(99,'(A1,1x,6F25.15)') 'F',&
                 & (CHAINP%BEADS(IND,:)+CHAINP%BEADS(IND-1,:))/2, &
                 & (CHAINP%BEADS(IND,:)+CHAINP%BEADS(IND-1,:))/2 + BRANCH
         ELSE
            DO I = 1,BP%NFILL
               BEAD = BP%CENT + QUAT2PT(BP%ORIENT*PTQUAT(BP%FILLBEADS(I,:))/BP%ORIENT)
               BRANCH = QUAT2PT(BP%ORIENT*BP%FILLQ(I)*PTQUAT((/1D0,0D0,0D0/))&
                    & /(BP%ORIENT*BP%FILLQ(I)))
               WRITE(99,'(A1,1x,6F25.15)') 'F', BEAD, BEAD+BRANCH                
            ENDDO
         ENDIF
      ENDIF

      ! place the linker beads
      IF (J.EQ.CHAINP%NBEND) THEN
         MAXIND = CHAINP%NBEAD
      ELSE
         MAXIND = CHAINP%BENDIND(J+1)
      ENDIF

      DO I = IND,MAXIND
         ! place branch according to quaternion
         QAX = PTQUAT((/1D0,0D0,0D0/))
         QAX = CHAINP%BEADQ(I)*QAX/CHAINP%BEADQ(I)                  
         BRANCH = CHAINP%BEADS(I,:) + QUAT2PT(QAX)

         WRITE(99,'(A1,1x,6F25.15)') 'A',CHAINP%BEADS(I,:), BRANCH
      ENDDo
   ENDDO

   ! write out the centers and orientations (as a quaternion) for each bend    
   ! also write out flexible tails with the last letter indicating whether
   ! the bead should be connected to the previous one
   DO B = 1,CHAINP%NBEND
      BP => CHAINP%BENDS(B)
      WRITE(99,'(A1,1x,7F25.15)'), 'B',BP%CENT, BP%ORIENT%W, BP%ORIENT%X, BP%ORIENT%Y, BP%ORIENT%Z

      ! write out marked points
      DO I = 1,BP%NMARKPOINT
         POS = BP%CENT +  QUAT2PT(BP%ORIENT*PTQUAT(BP%MARKPOINTS(I,:))/BP%ORIENT)   
         WRITE(99,'(A,1x,3F25.15)'), 'M'//BP%MARKPOINTNAME(I), POS
      ENDDO

      ! write out disco charges
      IF (BP%NCHARGE.GT.0) THEN
         DO I = 1,BP%NCHARGE
            POS = BP%CENT +  QUAT2PT(BP%ORIENT*PTQUAT(BP%CHARGEPOS(I,:))/BP%ORIENT)   
            WRITE(99,'(A1,1x,3F25.15)'), 'C', POS
         ENDDO
      ENDIF

      ! write out tail positions
      IF (BP%NTAIL.GT.0) THEN
         COUNT = 0
         DO I = 1,BP%NTAIL
            DO J = 1,BP%TAILNBEAD(I)
               COUNT = COUNT + 1
               POS = BP%FLEXTAILS(3*(COUNT-1)+1:3*COUNT)
               POS = BP%CENT +  QUAT2PT(BP%ORIENT*PTQUAT(POS)/BP%ORIENT)   
               IF (J.EQ.1) THEN
                  TMPI = 0
               ELSE
                  TMPI = 1
               ENDIF
               WRITE(99,'(A1,1x,3F25.15,I3)'), 'T',POS,TMPI
            ENDDO
         ENDDO
      ENDIF
   ENDDO

   CLOSE(99)

 END SUBROUTINE OUTPUTCHAIN

 SUBROUTINE CREATECHAIN(CHAINP,NBENDIN)
   ! create a chain object pointed to by the pointer CHAINP; 
   ! using parameters from keywords
   ! allocate the various arrays and set the parameters
   ! but DO NOT set the current geometry
   ! if the number of bends is not supplied explicitly, it is assumed to be 2
   USE KEYS
   IMPLICIT NONE
   TYPE(CHAIN), POINTER :: CHAINP
   INTEGER, INTENT(IN), OPTIONAL :: NBENDIN
   INTEGER :: NBEAD, B, J, I, COUNT
   DOUBLE PRECISION :: ROTMAT(3,3)
   TYPE(BEND), POINTER :: BP
   INTEGER :: NBEND

   IF (.NOT.ASSOCIATED(CHAINP)) THEN
      print*, 'ERROR in CREATECHAIN: CHAINP doesnt point to anything'
      STOP 1
   ENDIF

   IF (PRESENT(NBENDIN)) THEN
      NBEND = NBENDIN
   ELSE
      NBEND = 2
   ENDIF

   ! basic geometry parameters
   CHAINP%LP = LP
   CHAINP%LS = BENDSEP/NSEGPERLINK
   CHAINP%LTW = LTW
   CHAINP%LSTRETCH = LSTRETCH
   CHAINP%TWIST = TWIST
   CHAINP%TENSION = TENSION

   ! elastic and steric energetics
   CHAINP%BENDSTERICS = BENDSTERICS
   CHAINP%SEGSTERICS = SEGSTERICS
   CHAINP%BENDSEGSTERICS = BENDSEGSTERICS
   CHAINP%ESTERIC = ESTERIC
   CHAINP%DNARAD = DNARAD
   CHAINP%NINTERACT = MAXBENDINTER

   ! langowski potential stuff
   CHAINP%SIG0 = SIG0
   CHAINP%EPS0 = EPS0
   CHAINP%LJX = LJX
   CHAINP%LJXP = LJXP

   ! disco potential stuff
   CHAINP%EPSEV = EPSEV
   CHAINP%EPSEVT = EPSEVT
   CHAINP%SIGCC = SIGCC
   CHAINP%SIGCL = SIGCL
   CHAINP%SIGTC = SIGTC
   CHAINP%SIGTL = SIGTL
   CHAINP%SIGTT = SIGTT
   CHAINP%LDEBYEINV = LDEBYEINV

   IF (USEFIXDIAM) CHAINP%DIAM = FIXDIAMETER

   ! overall number of beads for outputting the chain
   NBEAD = NTAIL1SEG + NTAIL2SEG + 2 + (NBEND-1)*(NSEGPERLINK + 1)

   ! allocate arrays (if not already done)
   IF (.NOT.CHAINP%ARRAYSET) THEN
      IF (DISCOPOTENTIAL) THEN   
         IF (FILLINBEADS.GT.0) THEN
            CALL MAKEBARECHAIN(CHAINP,NBEAD,NBEND,NSEGPERLINK,FILLINBEADS,&
                 & NCHARGE,NTAIL,TAILNBEAD,NMARKPOINT)
         ELSE
            CALL MAKEBARECHAIN(CHAINP,NBEAD,NBEND,NSEGPERLINK,NCHARGE=NCHARGE,&
                 & NTAIL=NTAIL,TAILNBEAD=TAILNBEAD,NMARK=NMARKPOINT)
         ENDIF
      ELSE
         IF (FILLINBEADS.GT.0) THEN
            CALL MAKEBARECHAIN(CHAINP,NBEAD,NBEND,NSEGPERLINK,FILLINBEADS,NCHARGE,&
                 & NMARK=NMARKPOINT)
         ELSE
            CALL MAKEBARECHAIN(CHAINP,NBEAD,NBEND,NSEGPERLINK,NCHARGE=NCHARGE,&
                 & NMARK=NMARKPOINT)
         ENDIF
      ENDIF
   ENDIF

   ! Set up parameters for the bends (nucleosomes)
   DO B = 1,CHAINP%NBEND
      BP=>CHAINP%BENDS(B)
      BP%RADIUS = STERRADIUS
      BP%HEIGHT = STERHEIGHT
      BP%MAXRAD = SQRT(STERRADIUS**2 + (STERHEIGHT/2)**2)

      IF (SPIRALBENDS) THEN
         ! filler beads and edge conditions based on an idealized spiral
         CALL SETSPIRALBEND(BP,SPLEN,SPRADIUS,SPHEIGHT,SPTWIST,SPTWIST0,SPHAND) 
      ELSE
         ROTMAT(:,3) = BENDTANP
         ROTMAT(:,1) = XAXP
         CALL CROSS_PRODUCT(BENDTANP,XAXP,ROTMAT(:,2))
         BP%TP = ROTMAT2QUAT(ROTMAT)

         ROTMAT(:,3) = BENDTANM
         ROTMAT(:,1) = XAXM
         CALL CROSS_PRODUCT(BENDTANM,XAXM,ROTMAT(:,2))       
         BP%TM = ROTMAT2QUAT(ROTMAT)       

         BP%POSP = POSP
         BP%POSM = POSM

         ! filler beads directly from parameter file
         BP%NFILL = FILLINBEADS
         IF (B.EQ.1) THEN
            DO I = 1,FILLINBEADS
               BP%FILLBEADS(I,:) = FILLBEADLIST(I,1:3)                          
               ROTMAT(:,1) = FILLBEADLIST(I,4:6) - FILLBEADLIST(I,1:3)
               CALL NORMALIZE(ROTMAT(:,1))
               IF (I.EQ.FILLINBEADS) THEN
                  IF (ABS(ROTMAT(3,1)).GE.1D0-EPSILON(1D0)) THEN
                     CALL CROSS_PRODUCT(ROTMAT(:,1),(/0D0,1D0,0D0/),ROTMAT(:,3))
                  ELSE
                     CALL CROSS_PRODUCT(ROTMAT(:,1),(/0D0,0D0,1D0/),ROTMAT(:,3))
                  ENDIF
               ELSE
                  ROTMAT(:,3) = FILLBEADLIST(I+1,1:3)-FILLBEADLIST(I,1:3)
               ENDIF
               ! orthogonolize x and z axes
               ROTMAT(:,3) = ROTMAT(:,3) - DOT_PRODUCT(ROTMAT(:,3),ROTMAT(:,1))*ROTMAT(:,1)
               CALL NORMALIZE(ROTMAT(:,3))
               CALL CROSS_PRODUCT(ROTMAT(:,3),ROTMAT(:,1),ROTMAT(:,2))
               BP%FILLQ(I) = ROTMAT2QUAT(ROTMAT)
            ENDDO
         ELSE
            ! relative filler bead positions same for all bends
            BP%FILLBEADS = CHAINP%BENDS(1)%FILLBEADS
            BP%FILLQ = CHAINP%BENDS(1)%FILLQ
         END IF
      ENDIF

      ! DISCO potential stuff
      BP%NCHARGE = NCHARGE
      IF (NCHARGE.GT.0) THEN
         BP%CHARGEPOS = CHARGEPOS(1:NCHARGE,:)
         BP%CHARGEMAG = CHARGEMAG(1:NCHARGE)     
      ELSE
         NULLIFY(BP%CHARGEPOS); NULLIFY(BP%CHARGEMAG)
      ENDIF
      
      ! set up flexible tail stuff
      IF (USEFLEXTAILS) THEN          
         BP%FLEXTAILS = 0D0
         COUNT = 0
         DO I = 1,NTAIL
            BP%TAILPOS(I,:) = TAILPOS(I,:)             
            DO J = 1,TAILNBEAD(I)
               COUNT = COUNT + 1
               BP%TAILCHARGES(COUNT) = TAILBEADINFO(I,J,1)
               BP%TAILL0(COUNT) = TAILBEADINFO(I,J,2)
               BP%TAILKL(COUNT) = TAILBEADINFO(I,J,3)
               BP%TAILANG0(COUNT) = TAILBEADINFO(I,J,4)
               BP%TAILKA(COUNT) = TAILBEADINFO(I,J,5)
            ENDDO
         ENDDO
      ELSE
         NULLIFY(BP%TAILPOS); NULLIFY(BP%TAILCHARGES)
         NULLIFY(BP%TAILL0); NULLIFY(BP%TAILKL)
         NULLIFY(BP%TAILANG0); NULLIFY(BP%TAILKA)
      ENDIF

      ! mark points on bends
      DO I = 1,NMARKPOINT
         BP%MARKPOINTS(I,:) = MARKPOINTS(I,:)
         BP%MARKPOINTNAME(I) = MARKPOINTNAME(I)
      ENDDO
   ENDDO

   ! location of bends on the chain
   CHAINP%BENDIND(1) = NTAIL1SEG + 1
   DO B = 2,NBEND
      CHAINP%BENDIND(B) = CHAINP%BENDIND(B-1) + NSEGPERLINK+1
   ENDDO

 END SUBROUTINE CREATECHAIN


 SUBROUTINE MAKEBARECHAIN(CHAINP,NBEAD,NBEND,NSEG,NFILL,NCHARGE,NTAIL,TAILNBEAD,NMARK)
   ! create a chain object, allocating the appropriate arrays
   ! and setting the minimal number of parameters
   ! NOTE: this does not fill in any numerical parameters
   ! if NFILL is provided, then also allocate space for NFILL filler beads in each bend
   ! if NCHARGE is provided, allocate charge arrays for disco
   ! if NTAIL and TAILNBEAD are provided, allocate arrays for flexible tails

   TYPE(CHAIN), POINTER :: CHAINP
   INTEGER, INTENT(IN) :: NBEAD, NBEND, NSEG
   INTEGER, INTENT(IN), OPTIONAL :: NFILL,NCHARGE,NTAIL, TAILNBEAD(:),NMARK
   INTEGER :: B
   TYPE(BEND), POINTER :: BP    
   LOGICAL :: DOMARK

   IF (.NOT.ASSOCIATED(CHAINP)) THEN
      print*, 'ERROR in MAKEBARECHAIN: CHAINP doesnt point to anything'
      STOP 1
   ENDIF

   CHAINP%NSEG = NSEG
   CHAINP%NBEAD = NBEAD; CHAINP%NBEND = NBEND    

   ALLOCATE(CHAINP%BEADS(NBEAD,3),CHAINP%BEADQ(NBEAD))
   ALLOCATE(CHAINP%BENDS(NBEND),CHAINP%BENDIND(NBEND+1))    

   ! number of coords defining DNA configuration
   ! 6 coordinates for the helix (hpernuc,theta,rad,alpha,beta,gamma)
   ! 4 coordinates for each linker segment (alpha+gamma and bead position) 
   ! except last linker segment has alpha+gamma coordinates only
   CHAINP%NCRDDNA = 4*NSEG+3
   ! possibly extra coords for flexible tails
   IF (PRESENT(NTAIL).AND.PRESENT(TAILNBEAD)) THEN
      CHAINP%NCRD = CHAINP%NCRDDNA + SUM(TAILNBEAD(1:NTAIL))*3
   ELSE
      CHAINP%NCRD = CHAINP%NCRDDNA
   ENDIF

   ALLOCATE(CHAINP%VEC(CHAINP%NCRD),CHAINP%GRAD(CHAINP%NCRD))

   IF (PRESENT(NFILL)) THEN
      !  allocate locations of filler beads
      DO B = 1,CHAINP%NBEND
         CHAINP%BENDS(B)%NFILL = NFILL
         ALLOCATE(CHAINP%BENDS(B)%FILLBEADS(NFILL,3),CHAINP%BENDS(B)%FILLQ(NFILL))          
      ENDDO
   ELSE
      DO B = 1,CHAINP%NBEND
         CHAINP%BENDS(B)%NFILL = -1
         NULLIFY(CHAINP%BENDS(B)%FILLBEADS); NULLIFY(CHAINP%BENDS(B)%FILLQ)
      ENDDO
   ENDIF

   ! Marked points
   IF (PRESENT(NMARK)) THEN
      DOMARK = NMARK.GT.0
   ELSE
      DOMARK= .FALSE.
   ENDIF
   IF (DOMARK) THEN
      DO B = 1,CHAINP%NBEND
         CHAINP%BENDS(B)%NMARKPOINT = NMARK
         ALLOCATE(CHAINP%BENDS(B)%MARKPOINTS(NMARK,3),CHAINP%BENDS(B)%MARKPOINTNAME(NMARK))
      ENDDO
   ELSE
      DO B = 1,CHAINP%NBEND
         CHAINP%BENDS(B)%NMARKPOINT = 0
         NULLIFY(CHAINP%BENDS(B)%MARKPOINTS, CHAINP%BENDS(B)%MARKPOINTNAME)
      ENDDO
   ENDIF

    ! nucleosome charges for disco potential
   IF (PRESENT(NCHARGE)) THEN
      IF (NCHARGE.GT.0) THEN
         DO B = 1,CHAINP%NBEND
            CHAINP%BENDS(B)%NCHARGE = NCHARGE    
            ALLOCATE(CHAINP%BENDS(B)%CHARGEPOS(NCHARGE,3),&
                 & CHAINP%BENDS(B)%CHARGEMAG(NCHARGE))
         ENDDO
      ELSE
         DO B = 1,CHAINP%NBEND
            CHAINP%BENDS(B)%NCHARGE = 0
            NULLIFY(CHAINP%BENDS(B)%CHARGEPOS)
            NULLIFY(CHAINP%BENDS(B)%CHARGEMAG)
         ENDDO
      ENDIF
   ELSE
      DO B = 1,CHAINP%NBEND
         CHAINP%BENDS(B)%NCHARGE = 0
         NULLIFY(CHAINP%BENDS(B)%CHARGEPOS)
         NULLIFY(CHAINP%BENDS(B)%CHARGEMAG)
      ENDDO
   ENDIF

   IF (PRESENT(NTAIL).AND.PRESENT(TAILNBEAD)) THEN
      IF (NTAIL.GT.0) THEN
         ! set up arrays for flexible tails
         DO B = 1,CHAINP%NBEND
            BP=> CHAINP%BENDS(B)

            BP%TOTTAILBEAD = SUM(TAILNBEAD(1:NTAIL))
            ALLOCATE(BP%TAILNBEAD(NTAIL))
            ALLOCATE(BP%TAILL0(BP%TOTTAILBEAD),BP%TAILKL(BP%TOTTAILBEAD),BP%TAILANG0(BP%TOTTAILBEAD),&
                 & BP%TAILKA(BP%TOTTAILBEAD),BP%TAILPOS(NTAIL,3),BP%TAILCHARGES(BP%TOTTAILBEAD))
            BP%TAILNBEAD = TAILNBEAD(1:NTAIL)
            BP%NTAIL = NTAIL
            BP%FLEXTAILS => CHAINP%VEC(CHAINP%NCRDDNA+1:CHAINP%NCRD)          
         ENDDO
      ELSE
         DO B = 1,CHAINP%NBEND
            CHAINP%BENDS(B)%NTAIL = 0
            CHAINP%BENDS(B)%TOTTAILBEAD = 0
         ENDDO
      ENDIF
   ELSE
      DO B = 1,CHAINP%NBEND
         CHAINP%BENDS(B)%NTAIL = 0
         CHAINP%BENDS(B)%TOTTAILBEAD = 0
      ENDDO
   ENDIF

   CHAINP%ARRAYSET = .TRUE.    
   DO B = 1,CHAINP%NBEND
      CHAINP%BENDS(B)%ARRAYSET=.TRUE.
   ENDDO
 END SUBROUTINE MAKEBARECHAIN

 SUBROUTINE CLEANUPCHAIN(CHAINP)
   ! clean up allocated arrays for a chain object
   IMPLICIT NONE
   TYPE(CHAIN), POINTER :: CHAINP
   INTEGER :: B
   TYPE(BEND), POINTER :: BP

   DO B = 1,CHAINP%NBEND
      BP=>CHAINP%BENDS(B)
      IF (BP%NFILL.GT. 0) THEN
         DEALLOCATE(BP%FILLBEADS,BP%FILLQ)
      ENDIF
      IF (BP%NCHARGE.GT.0) DEALLOCATE(BP%CHARGEPOS,BP%CHARGEMAG)
      IF (BP%NTAIL.GT.0) THEN
         DEALLOCATE(BP%FLEXTAILS, BP%TAILNBEAD, BP%TAILL0, BP%TAILKL, &
              & BP%TAILANG0, BP%TAILKA, BP%TAILPOS, BP%TAILCHARGES, BP%MARKPOINTS,&
              & BP%MARKPOINTNAME)
      ENDIF
   ENDDO

   DEALLOCATE(CHAINP%BEADS,CHAINP%BENDIND,CHAINP%BENDS,CHAINP%BEADQ,CHAINP%VEC, CHAINP%GRAD)    

   CHAINP%ARRAYSET = .FALSE.

 END SUBROUTINE CLEANUPCHAIN

 SUBROUTINE COPYCHAIN(CHAIN1, CHAIN2, COPYALL)
   ! copy data from one chain object to another (from 1 to 2)
   ! if COPYALL is false, only copy the current state of the chain, 
   ! not the various parameters
   ! WARNING: this assumes the appropriate arrays have already been allocated
   ! and that the chains are the same size (same NCRD, NBEAD, NBEND, NSEG)
   TYPE(CHAIN), POINTER :: CHAIN1, CHAIN2
   LOGICAl, INTENT(IN) :: COPYALL
   INTEGER :: B
   TYPE(BEND), POINTER :: BP1, BP2   

   IF (CHAIN1%NBEAD.NE.CHAIN2%NBEAD&
        & .OR.CHAIN1%NBEND.NE.CHAIN2%NBEND &
        & .OR. CHAIN1%NCRD.NE.CHAIN2%NCRD&
        & .OR. CHAIN1%NSEG.NE.CHAIN2%NSEG) THEN
      print*, 'ERROR in COPYCHAIN: chains must have the same numbers of bends, beads, coords, segments'
      PRINT*, CHAIN1%NBEAD, CHAIN1%NBEND, CHAIN1%NCRD, CHAIN1%NSEG 
      PRINT*, CHAIN2%NBEAD, CHAIN2%NBEND, CHAIN2%NCRD, CHAIN1%NSEG
      STOP 1
   ELSEIF (.NOT.CHAIN2%ARRAYSET.OR..NOT.CHAIN1%ARRAYSET) THEN
      PRINT*, 'ERROR IN COPYCHAIN: both chains must have their arrays set to be properly copied'
      PRINT*, CHAIN1%ARRAYSET, CHAIN2%ARRAYSET
      STOP 1
   ENDIF

   CHAIN2%BEADS = CHAIN1%BEADS
   CHAIN2%BEADQ = CHAIN1%BEADQ    
   CHAIN2%ENERGY = CHAIN1%ENERGY
   CHAIN2%GRAD = CHAIN1%GRAD
   CHAIN2%VEC = CHAIN1%VEC
   CHAIN2%DIAM = CHAIN1%DIAM
   CHAIN2%ENERGYPARTS = CHAIN1%ENERGYPARTS

   DO B = 1,CHAIN1%NBEND
      BP1 => CHAIN1%BENDS(B); BP2 => CHAIN2%BENDS(B)
      CALL COPYBEND(BP1, BP2, COPYALL)
   ENDDO


   IF (COPYALL) THEN
      CHAIN2%BENDIND = CHAIN1%BENDIND
      CHAIN2%LP = CHAIN1%LP
      CHAIN2%LS = CHAIN1%LS
      CHAIN2%LTW = CHAIN1%LTW
      CHAIN2%LSTRETCH = CHAIN1%LSTRETCH
      CHAIN2%TWIST = CHAIN1%TWIST
      CHAIN2%ESTERIC = CHAIN1%ESTERIC
      CHAIN2%BENDSTERICS = CHAIN1%BENDSTERICS
      CHAIN2%SEGSTERICS = CHAIN1%SEGSTERICS
      CHAIN2%BENDSEGSTERICS = CHAIN1%BENDSEGSTERICS
      CHAIN2%TENSION = CHAIN1%TENSION
      CHAIN2%SIGTT = CHAIN1%SIGTT
      CHAIN2%SIGTC = CHAIN1%SIGTC
      CHAIN2%SIGCC = CHAIN1%SIGCC
      CHAIN2%SIGTL = CHAIN1%SIGTL
      CHAIN2%SIGCL = CHAIN1%SIGCL
      CHAIN2%EPSEV = CHAIN1%EPSEV
      CHAIN2%EPSEVT = CHAIN1%EPSEVT
      CHAIN2%LJX = CHAIN1%LJX
      CHAIN2%LJXP = CHAIN1%LJXP
      CHAIN2%EPS0 = CHAIN1%EPS0
      CHAIN2%SIG0 = CHAIN1%SIG0       
   ENDIF
 END SUBROUTINE COPYCHAIN

 SUBROUTINE COPYBEND(BEND1, BEND2, COPYALL)
   ! copy data from one bend object to another (from 1 to 2)
   ! if COPYALL is false, only copy the current position / orientation
   ! if COPYALL is true, copy everything
   ! Assumes the two bends have arrays all of the same size

   TYPE(BEND), POINTER :: BEND1, BEND2
   LOGICAL, INTENT(IN) :: COPYALL

   IF (.NOT.BEND2%ARRAYSET.OR..NOT.BEND1%ARRAYSET) THEN       
      PRINT*, 'ERROR IN COPYBEND: both bends must have their arrays set to be properly copied'
      PRINT*, BEND1%ARRAYSET, BEND2%ARRAYSET
      STOP 1
   ELSEIF (BEND2%NFILL.NE.BEND1%NFILL.OR.BEND2%NTAIL.NE.BEND1%NTAIL&
        & .OR.BEND1%NCHARGE.NE.BEND2%NCHARGE) THEN
      PRINT*, 'ERROR IN COPYBEND: bends do not have same size arrays'
      PRINT*, BEND1%NFILL, BEND1%NTAIL, BEND1%NCHARGE
      PRINT*, BEND2%NFILL, BEND2%NTAIL, BEND2%NCHARGE
      STOP 1
   ENDIF

   BEND2%ORIENT = BEND1%ORIENT
   BEND2%CENT = BEND1%CENT

   IF (BEND2%NFILL.GT.0) THEN
      BEND2%FILLBEADS = BEND1%FILLBEADS
   ENDIF
   IF (BEND2%NTAIL.GT.0) THEN
      BEND2%TAILPOS = BEND1%TAILPOS
   ENDIF

   IF (COPYALL) THEN
      BEND2%TP = BEND1%TP; BEND2%TM = BEND1%TM
      BEND2%POSP = BEND1%POSP; BEND2%POSM = BEND1%POSM
      BEND2%RADIUS = BEND1%RADIUS
      BEND2%HEIGHT = BEND1%HEIGHT
      BEND2%MAXRAD = BEND1%MAXRAD
      IF (BEND2%NFILL.GT.0) BEND2%FILLQ = BEND1%FILLQ

      ! disco stuff
      IF (BEND1%NCHARGE.GT.0) THEN
         BEND2%CHARGEPOS = BEND1%CHARGEPOS
         BEND2%CHARGEMAG = BEND1%CHARGEMAG
      ENDIF

      IF (BEND2%NTAIL.GT.0) THEN
         BEND2%FLEXTAILS = BEND1%FLEXTAILS
         BEND2%TAILNBEAD = BEND1%TAILNBEAD
         BEND2%TOTTAILBEAD = BEND1%TOTTAILBEAD
         BEND2%TAILL0 = BEND1%TAILL0
         BEND2%TAILKL = BEND1%TAILKL
         BEND2%TAILANG0 = BEND1%TAILKA 
         BEND2%TAILCHARGES = BEND1%TAILCHARGES
      ENDIF

      IF (BEND2%NMARKPOINT.GT.0) THEN
         IF (BEND2%NMARKPOINT.NE.BEND1%NMARKPOINT) THEN
            PRINT*, 'ERROR IN COPYBEND: the two bends have different numbers of marked points'
            STOP 1
         ENDIF
         BEND2%MARKPOINTS = BEND1%MARKPOINTS
         BEND2%MARKPOINTNAME = BEND1%MARKPOINTNAME
      ENDIF
   ENDIF
 END SUBROUTINE COPYBEND

 SUBROUTINE SETSPIRALBEND(BENDP,LB,RADIUS,HEIGHT,TWIST,TWIST0,SPHAND)
   ! fill up the FILLBEADS and FILLQ arrays of a bend object
   ! according to an idealized spiral
   ! also set the end conditions (POSP,POSM,BENDTANP, BENDTANM, XAXP, XAXM according to the spiral)
   ! BENDP%NFILL is the number of filler beads to use
   ! RADIUS, HEIGHT are the radius and height per turn of the spiral
   ! TWIST is the twist of the DNA within the spiral
   ! also sets the NFILL, LB parameters of the bend
   ! TWIST0 is the twist position half-way along the spiral
   ! SPHAND is positive for a RH helix and negative for a LH helix
   IMPLICIT NONE
   TYPE(BEND), POINTER :: BENDP
   INTEGER, INTENT(IN) :: SPHAND
   DOUBLE PRECISION, INTENT(IN) :: LB,RADIUS, HEIGHT, TWIST, TWIST0
   DOUBLE PRECISION :: LT, SP, CTS, CSL,STS,SSL, HSL
   INTEGER :: I, NB
   DOUBLE PRECISION :: RS(3), US(3), XS(3), ROTMAT(3,3)
   TYPE(QUATERNION) :: DNAQ
   
   NB = BENDP%NFILL

   ! length per turn
   LT = sqrt((2*pi*RADIUS)**2+HEIGHT**2)    
   ! space per bead
   SP =LB/(NB-1)

   DO I = 0,NB-1
      CTS = COS(TWIST0 - TWIST*LB/2 + TWIST*I*SP); STS = SIN(TWIST0 - TWIST*LB/2 + TWIST*I*SP)
      CSL = COS(2*PI/LT*(SP*I-LB/2)); SSL = SIN(2*PI/LT*(SP*I-LB/2)); HSL = HEIGHT*(SP*I-LB/2)/LT

      ! bead position
      RS = (/RADIUS*CSL,RADIUS*SSL*SPHAND, HSL/)
      ! bead tangent
      US = (/-2*PI*RADIUS/LT*SSL,2*PI*RADIUS/LT*CSL*SPHAND,HEIGHT/LT/)
      ! bead x-axis
      XS = (/-CSL*CTS + HEIGHT/LT*SSL*STS*SPHAND,&
           & -SSL*CTS*SPHAND-HEIGHT/LT*CSL*STS,2*PI*RADIUS/LT*STS*SPHAND/)
      ! quaternion for DNA axis orientation
      ROTMAT(:,3) = US; ROTMAT(:,1) = XS
      CALL CROSS_PRODUCT(ROTMAT(:,3),ROTMAT(:,1),ROTMAT(:,2))
      DNAQ = ROTMAT2QUAT(ROTMAT)

      BENDP%FILLBEADS(I+1,:) = RS
      BENDP%FILLQ(I+1) = DNAQ
   ENDDO

   BENDP%POSM = BENDP%FILLBEADS(1,:)
   BENDP%POSP = BENDP%FILLBEADS(NB,:)
   BENDP%TM = BENDP%FILLQ(1)
   BENDP%TP = BENDP%FILLQ(NB)

 END SUBROUTINE SETSPIRALBEND
END MODULE CHAINUTILS
